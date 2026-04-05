# 🧬 Protein Hydrophobicity Profiler v2

Kyte-Doolittle 疎水性スケールと Hopp-Woods 抗原性スケールを用いて、タンパク質配列の**膜貫通領域**と**抗原性サイト**を予測・可視化する Streamlit Web アプリケーション。

ESMFold API による3D構造予測と、残基レベルのスコアマッピングに対応。

## Features

- **Dual-scale analysis** — Kyte-Doolittle (疎水性) + Hopp-Woods (抗原性) を同時解析
- **2D profiling** — 上下2段のプロファイルプロットで膜貫通領域・抗原性サイトをハイライト
- **KD vs HW scatter** — 2スケールの相関を散布図で可視化し、領域を色分け
- **ESMFold 3D mapping** — API経由で構造予測 → スコアグラデーション or 領域ハイライトで3D表示
- **Interactive parameters** — ウィンドウサイズ・閾値をサイドバーでリアルタイム調整
- **CSV export** — 全残基のスコアをダウンロード可能
- **Multiple input** — FASTA upload / サンプル配列 / テキスト貼り付けに対応

## Quick Start

```bash
pip install -r requirements.txt
streamlit run hydrophobicity_profiler.py
```

## Screenshot

*(アプリ起動後にスクリーンショットを追加)*

## Algorithm

| Analysis | Method | Default Parameters |
|---|---|---|
| Hydrophobicity | Kyte-Doolittle sliding window | Window = 9 |
| Antigenicity | Hopp-Woods sliding window | Window = 9 |
| TM prediction | KD threshold-based continuous region | Score ≥ 1.6, Length ≥ 19 aa |
| Antigenic sites | HW threshold-based continuous region | Score ≥ 1.0, Length ≥ 6 aa |
| 3D structure | ESMFold API (≤ 400 aa) | B-factor replacement for coloring |

## Architecture

```
FASTA input
  ├── BioPython SeqIO ──→ Sequence parsing
  │
  ├── Kyte-Doolittle ──→ Hydrophobicity profile ──→ TM region detection
  │                                                       │
  ├── Hopp-Woods ──────→ Antigenicity profile ───→ Antigenic site detection
  │                                                       │
  ├── Matplotlib ──────→ 2D dual-panel plot              │
  │                   └→ KD vs HW scatter plot            │
  │                                                       │
  └── ESMFold API ─────→ PDB structure ─→ B-factor injection
                                        └→ py3Dmol 3D viewer (score gradient / region highlight)
```

## References

- Kyte, J. & Doolittle, R.F. (1982) "A simple method for displaying the hydropathic character of a protein." *J. Mol. Biol.* 157:105-132.
- Hopp, T.P. & Woods, K.R. (1981) "Prediction of protein antigenic determinants from amino acid sequences." *Proc. Natl. Acad. Sci. USA* 78:3824-3828.
- Lin, Z. et al. (2023) "Evolutionary-scale prediction of atomic-level protein structure with a language model." *Science* 379:1123-1130.

## Tech Stack

- Python 3.10+
- Streamlit — Web UI framework
- BioPython — FASTA parsing (SeqIO)
- NumPy — Vectorized sliding window (np.convolve)
- Matplotlib — 2D visualization
- py3Dmol / stmol — 3D molecular visualization
- ESMFold API — AI-based structure prediction
