"""
Protein Hydrophobicity Profiler v2.0
=====================================
Kyte-Doolittle 疎水性スケールおよび Hopp-Woods 抗原性スケールを用いて
FASTAファイルからタンパク質の膜貫通領域・抗原性サイトを2D/3Dで可視化する
Streamlit アプリケーション。

ESMFold API による構造予測 → py3Dmol による3D可視化に対応。

Requirements:
    pip install -r requirements.txt
Run:
    streamlit run hydrophobicity_profiler.py
"""

import io
import re
import numpy as np
import requests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO
import py3Dmol


# ============================================================
# 1. Amino Acid Property Scales
# ============================================================

# Kyte-Doolittle (1982) — 疎水性スケール
# 正の値 = 疎水性（膜内部に埋もれやすい）
# 負の値 = 親水性（表面に露出しやすい）
KD_SCALE: dict[str, float] = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}

# Hopp-Woods (1981) — 抗原性（親水性）スケール
# 正の値 = 親水性が高い（表面露出 → 抗原性サイト候補）
# 負の値 = 疎水性が高い（内部に埋没）
# ※ Kyte-Doolittle とは正負の向きが逆
HW_SCALE: dict[str, float] = {
    "A": -0.5, "R":  3.0, "N":  0.2, "D":  3.0, "C": -1.0,
    "Q":  0.2, "E":  3.0, "G":  0.0, "H": -0.5, "I": -1.8,
    "L": -1.8, "K":  3.0, "M": -1.3, "F": -2.5, "P":  0.0,
    "S":  0.3, "T": -0.4, "W": -3.4, "Y": -2.3, "V": -1.5,
}


# ============================================================
# 2. Core Analysis Functions
# ============================================================

def calc_profile(
    sequence: str,
    scale: dict[str, float],
    window: int = 9,
) -> tuple[np.ndarray, np.ndarray]:
    """スライディングウィンドウでプロパティスコアを計算する。

    Args:
        sequence: アミノ酸配列（1文字表記）
        scale: アミノ酸→スコアの辞書
        window: ウィンドウサイズ（奇数推奨）

    Returns:
        scores: 各ウィンドウ中心位置の平均スコア (ndarray)
        positions: 残基番号 1-indexed (ndarray)
    """
    half = window // 2
    seq_len = len(sequence)

    if seq_len < window:
        return np.array([]), np.array([])

    # 全残基のスコアを一括取得してから畳み込み（np.convolve で高速化）
    raw_scores = np.array([scale.get(aa, 0.0) for aa in sequence])
    kernel = np.ones(window) / window
    convolved = np.convolve(raw_scores, kernel, mode="valid")

    positions = np.arange(half + 1, seq_len - half + 1)  # 1-indexed
    return convolved, positions


def detect_regions(
    scores: np.ndarray,
    positions: np.ndarray,
    threshold: float,
    min_length: int,
    direction: str = "high",
) -> list[dict]:
    """閾値を超える（または下回る）連続領域を検出する。

    Args:
        scores: スコア配列
        positions: 残基番号配列
        threshold: 判定閾値
        min_length: 最小連続残基数
        direction: "high" = 閾値以上を検出 / "low" = 閾値以下を検出

    Returns:
        検出領域のリスト [{"start", "end", "length", "avg_score"}, ...]
    """
    if direction == "high":
        mask = scores >= threshold
    else:
        mask = scores <= threshold

    regions: list[dict] = []
    in_region = False
    start_idx = 0

    for i, hit in enumerate(mask):
        if hit and not in_region:
            in_region = True
            start_idx = i
        elif not hit and in_region:
            in_region = False
            length = i - start_idx
            if length >= min_length:
                regions.append({
                    "start": int(positions[start_idx]),
                    "end": int(positions[i - 1]),
                    "length": length,
                    "avg_score": float(np.mean(scores[start_idx:i])),
                })

    # 末端処理
    if in_region:
        length = len(scores) - start_idx
        if length >= min_length:
            regions.append({
                "start": int(positions[start_idx]),
                "end": int(positions[-1]),
                "length": length,
                "avg_score": float(np.mean(scores[start_idx:])),
            })

    return regions


def per_residue_scores(
    sequence: str,
    scale: dict[str, float],
    window: int = 9,
) -> np.ndarray:
    """全残基に対してスコアを返す（3D構造マッピング用）。

    ウィンドウ端の残基はスコアが割り当てられないため、
    端は最近傍のスコアで埋める（パディング）。

    Returns:
        全残基分のスコア配列 (length == len(sequence))
    """
    scores, _ = calc_profile(sequence, scale, window)
    if len(scores) == 0:
        return np.zeros(len(sequence))

    half = window // 2
    full = np.zeros(len(sequence))
    full[half : half + len(scores)] = scores
    full[:half] = scores[0]
    full[half + len(scores):] = scores[-1]

    return full


# ============================================================
# 3. Structure Prediction API Integration
# ============================================================

# ESM Atlas API（Meta公式エンドポイント）
ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_ESMFOLD_LENGTH = 400  # ESMFold API の推奨上限

# AlphaFold DB（UniProt ID から既存の予測構造を取得）
ALPHAFOLD_DB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"


def predict_structure_esmfold(sequence: str, timeout: int = 120) -> str | None:
    """ESMFold API を呼び出して PDB 文字列を取得する。

    SSL証明書の問題が既知のため、以下の順にフォールバックする:
      1. ESM Atlas API（通常接続）
      2. ESM Atlas API（SSL検証スキップ）

    Args:
        sequence: アミノ酸配列（最大 400 残基推奨）
        timeout: タイムアウト秒数

    Returns:
        PDB形式の文字列。失敗時は None。
    """
    headers = {"Content-Type": "application/x-www-form-urlencoded"}

    # --- 1st try: 通常のSSL検証あり ---
    try:
        response = requests.post(
            ESMFOLD_URL,
            data=sequence,
            headers=headers,
            timeout=timeout,
            verify=True,
        )
        response.raise_for_status()
        if "ATOM" in response.text:
            return response.text
    except requests.exceptions.SSLError:
        # SSL証明書エラー → verify=False でリトライ
        st.warning(
            "ESM Atlas APIのSSL証明書に問題があるため、"
            "SSL検証をスキップしてリトライします。"
        )
    except requests.exceptions.RequestException as e:
        st.warning(f"ESMFold API (1st attempt): {e}")

    # --- 2nd try: SSL検証スキップ ---
    try:
        import urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

        response = requests.post(
            ESMFOLD_URL,
            data=sequence,
            headers=headers,
            timeout=timeout,
            verify=False,
        )
        response.raise_for_status()
        if "ATOM" in response.text:
            st.success("ESMFold API (SSL skip) で構造を取得しました。")
            return response.text
    except requests.exceptions.RequestException as e:
        st.warning(f"ESMFold API (SSL skip): {e}")

    return None


def fetch_alphafold_structure(uniprot_id: str, timeout: int = 30) -> str | None:
    """AlphaFold Protein Structure Database から既存の予測構造を取得する。

    ESMFold APIが利用できない場合のフォールバック。
    UniProt IDが必要（例: P07550）。

    Args:
        uniprot_id: UniProt accession（例: "P07550"）
        timeout: タイムアウト秒数

    Returns:
        PDB形式の文字列。失敗時は None。
    """
    url = ALPHAFOLD_DB_URL.format(uniprot_id=uniprot_id)
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        if "ATOM" in response.text:
            return response.text
        return None
    except requests.exceptions.RequestException as e:
        st.warning(f"AlphaFold DB error: {e}")
        return None


def inject_bfactor(pdb_string: str, residue_scores: np.ndarray) -> str:
    """PDB の B-factor カラムをスコア値で置換する。

    py3Dmol で spectrum カラーリングする際に B-factor を利用して
    残基ごとのスコアをグラデーション表示する。

    Args:
        pdb_string: PDB形式テキスト
        residue_scores: 残基ごとのスコア配列

    Returns:
        B-factor 置換済みの PDB 文字列
    """
    s_min, s_max = residue_scores.min(), residue_scores.max()
    if s_max - s_min < 1e-6:
        scaled = np.full_like(residue_scores, 50.0)
    else:
        scaled = (residue_scores - s_min) / (s_max - s_min) * 100.0

    lines = []
    for line in pdb_string.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                res_seq = int(line[22:26].strip()) - 1  # 0-indexed
                if 0 <= res_seq < len(scaled):
                    bfactor = f"{scaled[res_seq]:6.2f}"
                    line = line[:60] + bfactor + line[66:]
            except (ValueError, IndexError):
                pass
        lines.append(line)

    return "\n".join(lines)


# ============================================================
# 4. 3D Visualization (py3Dmol)
# ============================================================

def show_structure_3d(
    pdb_string: str,
    residue_scores: np.ndarray,
    tm_regions: list[dict],
    antigenic_sites: list[dict],
    color_mode: str = "score",
    width: int = 700,
    height: int = 500,
):
    """py3Dmol で 3D 構造を描画する。

    Args:
        pdb_string: PDB 形式テキスト
        residue_scores: 残基ごとのスコア配列
        tm_regions: 膜貫通領域リスト
        antigenic_sites: 抗原性サイトリスト
        color_mode: "score" = スコアグラデーション / "regions" = 領域ハイライト
    """
    modified_pdb = inject_bfactor(pdb_string, residue_scores)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(modified_pdb, "pdb")

    if color_mode == "score":
        # B-factor ベースのグラデーション (blue=low → white → red=high)
        view.setStyle(
            {},
            {"cartoon": {
                "colorscheme": {
                    "prop": "b",
                    "gradient": "rwb",
                    "min": 0,
                    "max": 100,
                },
            }},
        )
    elif color_mode == "regions":
        # ベースカラー: ライトグレー
        view.setStyle({}, {"cartoon": {"color": "#D3D3D3"}})

        # 膜貫通領域: 赤
        for region in tm_regions:
            view.setStyle(
                {"resi": list(range(region["start"], region["end"] + 1))},
                {"cartoon": {"color": "#E24B4A"}},
            )

        # 抗原性サイト: 青
        for site in antigenic_sites:
            view.setStyle(
                {"resi": list(range(site["start"], site["end"] + 1))},
                {"cartoon": {"color": "#378ADD"}},
            )

    view.zoomTo()
    view.spin(False)

    # stmol を使わず streamlit.components.v1.html で直接レンダリング
    html_str = view._make_html()
    components.html(html_str, height=height, width=width, scrolling=False)


# ============================================================
# 5. 2D Visualization (Matplotlib)
# ============================================================

def plot_dual_profile(
    positions: np.ndarray,
    kd_scores: np.ndarray,
    hw_scores: np.ndarray,
    tm_regions: list[dict],
    antigenic_sites: list[dict],
    protein_name: str = "Protein",
    kd_tm_threshold: float = 1.6,
    hw_antigenic_threshold: float = 1.0,
) -> plt.Figure:
    """Kyte-Doolittle と Hopp-Woods を上下2段で可視化する。"""

    fig, (ax_kd, ax_hw) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # --- Upper panel: Kyte-Doolittle ---
    ax_kd.plot(positions, kd_scores, color="black", linewidth=0.8)
    ax_kd.fill_between(positions, kd_scores, 0,
                        where=kd_scores > 0, color="#FFCCCC", alpha=0.3)
    ax_kd.fill_between(positions, kd_scores, 0,
                        where=kd_scores <= 0, color="#CCE5FF", alpha=0.3)

    for region in tm_regions:
        ax_kd.axvspan(region["start"], region["end"],
                       color="#E24B4A", alpha=0.25)
        mid = (region["start"] + region["end"]) / 2
        ax_kd.annotate("TM", xy=(mid, np.max(kd_scores) * 0.85),
                         ha="center", fontsize=8, color="#A32D2D",
                         fontweight="bold")

    ax_kd.axhline(y=kd_tm_threshold, color="#E24B4A", linestyle="--",
                   linewidth=0.8, alpha=0.6)
    ax_kd.axhline(y=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.4)
    ax_kd.set_ylabel("Kyte-Doolittle score", fontsize=11)
    ax_kd.set_title(
        f"Hydrophobicity / Antigenicity profile — {protein_name}",
        fontsize=13, fontweight="bold",
    )

    tm_patch = mpatches.Patch(color="#E24B4A", alpha=0.25, label="TM candidate")
    kd_line = plt.Line2D([0], [0], color="#E24B4A", linestyle="--",
                          label=f"TM threshold ({kd_tm_threshold})")
    ax_kd.legend(handles=[tm_patch, kd_line], loc="upper right", fontsize=8)

    # --- Lower panel: Hopp-Woods ---
    ax_hw.plot(positions, hw_scores, color="#185FA5", linewidth=0.8)
    ax_hw.fill_between(positions, hw_scores, 0,
                        where=hw_scores > 0, color="#CCE5FF", alpha=0.3)
    ax_hw.fill_between(positions, hw_scores, 0,
                        where=hw_scores <= 0, color="#FFCCCC", alpha=0.3)

    for site in antigenic_sites:
        ax_hw.axvspan(site["start"], site["end"],
                       color="#378ADD", alpha=0.20)
        mid = (site["start"] + site["end"]) / 2
        ax_hw.annotate("Ag", xy=(mid, np.max(hw_scores) * 0.85),
                         ha="center", fontsize=8, color="#185FA5",
                         fontweight="bold")

    ax_hw.axhline(y=hw_antigenic_threshold, color="#378ADD", linestyle="--",
                   linewidth=0.8, alpha=0.6)
    ax_hw.axhline(y=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.4)
    ax_hw.set_xlabel("Residue position", fontsize=11)
    ax_hw.set_ylabel("Hopp-Woods score", fontsize=11)

    ag_patch = mpatches.Patch(color="#378ADD", alpha=0.20, label="Antigenic candidate")
    hw_line = plt.Line2D([0], [0], color="#378ADD", linestyle="--",
                          label=f"Antigenic threshold ({hw_antigenic_threshold})")
    ax_hw.legend(handles=[ag_patch, hw_line], loc="upper right", fontsize=8)

    fig.tight_layout()
    return fig


def plot_comparison_scatter(
    kd_full: np.ndarray,
    hw_full: np.ndarray,
    tm_regions: list[dict],
    antigenic_sites: list[dict],
    protein_name: str = "Protein",
) -> plt.Figure:
    """KD vs HW の散布図。TM / Antigenic 領域を色分けする。"""

    fig, ax = plt.subplots(figsize=(7, 6))

    n = len(kd_full)
    colors = np.full(n, "#AAAAAA")

    for region in tm_regions:
        for i in range(region["start"] - 1, min(region["end"], n)):
            colors[i] = "#E24B4A"

    for site in antigenic_sites:
        for i in range(site["start"] - 1, min(site["end"], n)):
            colors[i] = "#378ADD"

    ax.scatter(kd_full, hw_full, c=colors, s=8, alpha=0.6, edgecolors="none")

    # 凡例用ダミー
    ax.scatter([], [], c="#E24B4A", s=30, label="TM candidate")
    ax.scatter([], [], c="#378ADD", s=30, label="Antigenic candidate")
    ax.scatter([], [], c="#AAAAAA", s=30, label="Other")

    ax.set_xlabel("Kyte-Doolittle (hydrophobicity)", fontsize=11)
    ax.set_ylabel("Hopp-Woods (antigenicity)", fontsize=11)
    ax.set_title(f"KD vs HW scatter — {protein_name}", fontsize=13,
                  fontweight="bold")
    ax.legend(fontsize=9)
    ax.axhline(y=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.4)
    ax.axvline(x=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.4)

    fig.tight_layout()
    return fig


# ============================================================
# 6. Streamlit UI
# ============================================================

SAMPLE_FASTA = (
    ">sp|P07550|ADRB2_HUMAN Beta-2 adrenergic receptor (partial)\n"
    "MGQPGNGSAFLLAPNRSHAPDHDVTQQRDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAK"
    "FERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTAS"
    "IETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEA"
    "INCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHV"
    "QNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDN"
    "LIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGE"
    "QSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDNIDSQGRNCSTNDSLL\n"
)


def render_sidebar() -> dict:
    """サイドバーのパラメータ設定UI。設定値を辞書で返す。"""
    with st.sidebar:
        st.header("⚙️ Parameters")

        st.subheader("Sliding window")
        window_size = st.slider(
            "Window size", 5, 25, 9, step=2,
            help="スコア平均化ウィンドウ（奇数推奨）",
        )

        st.subheader("Kyte-Doolittle (TM detection)")
        kd_tm_thresh = st.slider(
            "KD — TM threshold", 0.0, 3.0, 1.6, step=0.1,
            help="KDスコアがこの値以上の連続領域を膜貫通候補とする",
        )
        kd_tm_minlen = st.slider(
            "KD — TM min length", 10, 30, 19, step=1,
        )

        st.subheader("Hopp-Woods (Antigenic detection)")
        hw_ag_thresh = st.slider(
            "HW — Antigenic threshold", 0.0, 3.0, 1.0, step=0.1,
            help="HWスコアがこの値以上の連続領域を抗原性サイト候補とする",
        )
        hw_ag_minlen = st.slider(
            "HW — Antigenic min length", 3, 15, 6, step=1,
        )

        st.divider()

        st.subheader("3D Structure")
        enable_3d = st.checkbox("Enable 3D structure view", value=False)
        color_3d = st.radio(
            "3D color mode",
            ["score", "regions"],
            format_func=lambda x: {
                "score": "Score gradient (blue → red)",
                "regions": "TM / Antigenic highlight",
            }[x],
            disabled=not enable_3d,
        )
        score_type_3d = st.radio(
            "Score for 3D coloring",
            ["Kyte-Doolittle", "Hopp-Woods"],
            disabled=not enable_3d,
        )

        st.divider()
        st.caption(
            "**References**\n\n"
            "Kyte & Doolittle (1982) *J. Mol. Biol.* 157:105-132\n\n"
            "Hopp & Woods (1981) *PNAS* 78:3824-3828\n\n"
            "Lin et al. (2023) ESMFold — *Science* 379:1123-1130\n\n"
            "Jumper et al. (2021) AlphaFold — *Nature* 596:583-589"
        )

    return {
        "window": window_size,
        "kd_tm_thresh": kd_tm_thresh,
        "kd_tm_minlen": kd_tm_minlen,
        "hw_ag_thresh": hw_ag_thresh,
        "hw_ag_minlen": hw_ag_minlen,
        "enable_3d": enable_3d,
        "color_3d": color_3d,
        "score_type_3d": score_type_3d,
    }


def render_analysis(record, params: dict):
    """1つの配列に対する解析・可視化をレンダリングする。"""

    seq_str = str(record.seq).upper()
    protein_name = record.description or record.id
    st.subheader(f"📊 {protein_name}")

    # --- 2D Profile Calculation ---
    kd_scores, positions = calc_profile(seq_str, KD_SCALE, params["window"])
    hw_scores, _ = calc_profile(seq_str, HW_SCALE, params["window"])

    if len(kd_scores) == 0:
        st.warning(
            f"配列が短すぎます（{len(seq_str)} aa < window {params['window']}）"
        )
        return

    # 領域検出
    # KD: 高スコア → 膜貫通候補
    tm_regions = detect_regions(
        kd_scores, positions,
        params["kd_tm_thresh"], params["kd_tm_minlen"],
        direction="high",
    )
    # HW: 高スコア → 抗原性サイト候補
    antigenic_sites = detect_regions(
        hw_scores, positions,
        params["hw_ag_thresh"], params["hw_ag_minlen"],
        direction="high",
    )

    # --- Tabs ---
    tab_2d, tab_scatter, tab_3d, tab_table = st.tabs([
        "📈 2D Profile", "🔄 KD vs HW", "🧊 3D Structure", "📋 Results",
    ])

    # ========== Tab 1: 2D Profile ==========
    with tab_2d:
        fig_2d = plot_dual_profile(
            positions, kd_scores, hw_scores,
            tm_regions, antigenic_sites,
            protein_name, params["kd_tm_thresh"], params["hw_ag_thresh"],
        )
        st.pyplot(fig_2d)
        plt.close(fig_2d)

    # ========== Tab 2: KD vs HW Scatter ==========
    with tab_scatter:
        kd_full = per_residue_scores(seq_str, KD_SCALE, params["window"])
        hw_full = per_residue_scores(seq_str, HW_SCALE, params["window"])
        fig_sc = plot_comparison_scatter(
            kd_full, hw_full, tm_regions, antigenic_sites, protein_name,
        )
        st.pyplot(fig_sc)
        plt.close(fig_sc)

        st.caption(
            "KD（疎水性）が高く HW（親水性）が低い残基 → 膜内部に位置する傾向。"
            "逆のパターン → 表面露出した抗原性サイト候補。"
        )

    # ========== Tab 3: 3D Structure ==========
    with tab_3d:
        if not params["enable_3d"]:
            st.info(
                "3D構造表示を有効にするには、サイドバーの "
                "**Enable 3D structure view** をオンにしてください。"
            )
        else:
            # 配列長チェック
            if len(seq_str) > MAX_ESMFOLD_LENGTH:
                st.warning(
                    f"配列長 {len(seq_str)} aa は ESMFold API の推奨上限 "
                    f"({MAX_ESMFOLD_LENGTH} aa) を超えています。"
                    f"先頭 {MAX_ESMFOLD_LENGTH} 残基で予測します。"
                )
                seq_for_fold = seq_str[:MAX_ESMFOLD_LENGTH]
            else:
                seq_for_fold = seq_str

            # PDB を session_state にキャッシュ（同一配列は再予測しない）
            cache_key = f"pdb_{hash(seq_for_fold)}"

            # --- 構造取得: 3段階フォールバック ---
            if cache_key not in st.session_state:

                # Stage 1: ESMFold API
                with st.spinner("ESMFold API で構造予測中..."):
                    pdb_str = predict_structure_esmfold(seq_for_fold)

                if pdb_str:
                    st.session_state[cache_key] = pdb_str
                    st.success("ESMFold で構造を取得しました。")
                else:
                    st.warning(
                        "ESMFold API から構造を取得できませんでした "
                        "（SSL証明書の問題が既知です）。\n\n"
                        "**代替手段:** 以下から構造を取得できます。"
                    )

            # Stage 2: AlphaFold DB（UniProt ID 入力）
            if cache_key not in st.session_state:
                st.markdown("---")
                st.markdown("**🔄 代替 1: AlphaFold DB から取得**")
                st.caption(
                    "FASTAヘッダにUniProt IDが含まれる場合は自動検出します。"
                    "手動入力も可能です（例: P07550）。"
                )

                # UniProt ID を FASTA ヘッダから自動抽出を試みる
                auto_uniprot = ""
                header = record.description or record.id or ""
                # sp|P07550|ADRB2_HUMAN のようなパターン
                match = re.search(
                    r"(?:sp|tr)\|([A-Z0-9]+)\|", header
                )
                if match:
                    auto_uniprot = match.group(1)

                uniprot_id = st.text_input(
                    "UniProt ID",
                    value=auto_uniprot,
                    placeholder="P07550",
                    key=f"uniprot_{protein_name}",
                )

                if st.button(
                    "📥 AlphaFold DB から取得",
                    key=f"af_btn_{protein_name}",
                ):
                    if uniprot_id.strip():
                        with st.spinner(
                            f"AlphaFold DB から {uniprot_id} を取得中..."
                        ):
                            pdb_str = fetch_alphafold_structure(
                                uniprot_id.strip()
                            )
                        if pdb_str:
                            st.session_state[cache_key] = pdb_str
                            st.success(
                                f"AlphaFold DB から構造を取得しました "
                                f"(UniProt: {uniprot_id})。"
                            )
                            st.rerun()
                        else:
                            st.error(
                                f"UniProt ID '{uniprot_id}' の構造が "
                                f"AlphaFold DB に見つかりませんでした。"
                            )
                    else:
                        st.error("UniProt ID を入力してください。")

            # Stage 3: 手動 PDB アップロード
            if cache_key not in st.session_state:
                st.markdown("---")
                st.markdown("**📁 代替 2: PDB ファイルを手動アップロード**")

            uploaded_pdb = st.file_uploader(
                "PDB file (manual upload — optional)",
                type=["pdb"],
                key=f"pdb_upload_{protein_name}",
            )
            if uploaded_pdb:
                st.session_state[cache_key] = (
                    uploaded_pdb.getvalue().decode("utf-8")
                )

            # --- 3D 描画 ---
            if cache_key in st.session_state:
                pdb_data = st.session_state[cache_key]

                # カラーリング用スコア
                scale_for_3d = (
                    KD_SCALE if params["score_type_3d"] == "Kyte-Doolittle"
                    else HW_SCALE
                )
                scores_3d = per_residue_scores(
                    seq_for_fold, scale_for_3d, params["window"]
                )

                mode_label = (
                    "gradient (blue=low → red=high)"
                    if params["color_3d"] == "score"
                    else "TM=red / Ag=blue / Other=gray"
                )
                st.markdown(
                    f"**Color:** {params['score_type_3d']} — {mode_label}"
                )

                show_structure_3d(
                    pdb_data, scores_3d,
                    tm_regions, antigenic_sites,
                    color_mode=params["color_3d"],
                )

    # ========== Tab 4: Results table ==========
    with tab_table:
        col1, col2 = st.columns(2)

        with col1:
            st.markdown("**🔴 Transmembrane region candidates (KD)**")
            if tm_regions:
                for i, r in enumerate(tm_regions, 1):
                    st.markdown(
                        f"- **TM{i}**: Pos {r['start']}–{r['end']} "
                        f"({r['length']} aa, avg KD: {r['avg_score']:.2f})"
                    )
            else:
                st.caption("膜貫通候補は検出されませんでした。")

        with col2:
            st.markdown("**🔵 Antigenic site candidates (HW)**")
            if antigenic_sites:
                for i, s in enumerate(antigenic_sites, 1):
                    st.markdown(
                        f"- **Ag{i}**: Pos {s['start']}–{s['end']} "
                        f"({s['length']} aa, avg HW: {s['avg_score']:.2f})"
                    )
            else:
                st.caption("抗原性サイト候補は検出されませんでした。")

        # 統計サマリー
        with st.expander("📈 Sequence statistics"):
            c1, c2, c3, c4, c5 = st.columns(5)
            c1.metric("Length", f"{len(seq_str)} aa")
            c2.metric("Mean KD", f"{np.mean(kd_scores):.2f}")
            c3.metric("Mean HW", f"{np.mean(hw_scores):.2f}")
            c4.metric("TM regions", len(tm_regions))
            c5.metric("Ag sites", len(antigenic_sites))

        # CSV エクスポート
        with st.expander("💾 Export scores (CSV)"):
            kd_full = per_residue_scores(seq_str, KD_SCALE, params["window"])
            hw_full = per_residue_scores(seq_str, HW_SCALE, params["window"])

            csv_lines = ["position,residue,kd_score,hw_score"]
            for i, aa in enumerate(seq_str):
                csv_lines.append(
                    f"{i+1},{aa},{kd_full[i]:.4f},{hw_full[i]:.4f}"
                )

            csv_text = "\n".join(csv_lines)
            st.download_button(
                label="📥 Download CSV",
                data=csv_text,
                file_name=f"{record.id}_hydrophobicity.csv",
                mime="text/csv",
            )

    st.divider()


def main():
    st.set_page_config(
        page_title="Protein Hydrophobicity Profiler v2",
        page_icon="🧬",
        layout="wide",
    )

    st.title("🧬 Protein Hydrophobicity Profiler v2")
    st.markdown(
        "**Kyte-Doolittle 疎水性スケール** と **Hopp-Woods 抗原性スケール** の "
        "デュアル解析。ESMFold による 3D 構造予測・スコアマッピングに対応。"
    )

    params = render_sidebar()

    # --- Input area ---
    upload_tab, sample_tab, paste_tab = st.tabs([
        "📁 Upload FASTA", "🧪 Sample (GPCR)", "📝 Paste sequence",
    ])

    with upload_tab:
        uploaded_file = st.file_uploader(
            "FASTAファイルをアップロード",
            type=["fasta", "fa", "faa", "txt"],
        )

    with sample_tab:
        st.markdown(
            "β2アドレナリン受容体（7回膜貫通型GPCR）のサンプル配列を使用します。"
        )
        if st.button("▶ Run with sample", type="primary"):
            st.session_state["use_sample"] = True

    with paste_tab:
        pasted = st.text_area(
            "アミノ酸配列を貼り付け（FASTA形式 or 生配列）",
            height=120,
            placeholder=">my_protein\nMKWVTFISLLFLFSSAYSRGV...",
        )
        if st.button("▶ Run with pasted sequence"):
            if pasted.strip():
                st.session_state["pasted_seq"] = pasted.strip()

    # --- Determine input source ---
    fasta_handle = None
    if uploaded_file is not None:
        fasta_handle = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    elif st.session_state.get("use_sample"):
        fasta_handle = io.StringIO(SAMPLE_FASTA)
    elif st.session_state.get("pasted_seq"):
        text = st.session_state["pasted_seq"]
        if not text.startswith(">"):
            text = ">pasted_sequence\n" + text
        fasta_handle = io.StringIO(text)

    if fasta_handle is None:
        st.info(
            "FASTAファイルをアップロード、サンプル配列を使用、"
            "または配列を貼り付けてください。"
        )
        return

    # --- Parse & analyze each record ---
    records = list(SeqIO.parse(fasta_handle, "fasta"))
    if not records:
        st.error("有効な配列が見つかりませんでした。FASTA形式を確認してください。")
        return

    for record in records:
        render_analysis(record, params)


if __name__ == "__main__":
    main()
