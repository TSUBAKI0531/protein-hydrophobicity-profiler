# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Activate environment
source venv/bin/activate

# Run the app
streamlit run hydrophobicity_profiler.py

# Install dependencies
venv/bin/python -m pip install -r requirements.txt
```

No test suite or linter is configured.

## Architecture

The entire application lives in a single file: `hydrophobicity_profiler.py` (~900 lines), structured in 6 sections:

1. **Amino Acid Property Scales** — `KD_SCALE` and `HW_SCALE` dictionaries (Kyte-Doolittle and Hopp-Woods)
2. **Core Analysis Functions** — `calc_profile()` uses `np.convolve` for sliding-window scoring; `detect_regions()` finds contiguous stretches above/below threshold; `per_residue_scores()` pads scores to full sequence length for 3D mapping
3. **Structure Prediction API Integration** — 3-stage fallback:
   - Stage 1/2: `predict_structure_esmfold()` — POST to ESM Atlas API, retries with `verify=False` on SSL error
   - Stage 3 (UI-driven): `fetch_alphafold_structure()` — GET from AlphaFold EBI DB by UniProt ID
   - Manual: PDB file upload via `st.file_uploader`
4. **3D Visualization** — `inject_bfactor()` replaces PDB B-factor column with normalized scores (0–100); `show_structure_3d()` renders via py3Dmol using `streamlit.components.v1.html` directly (no stmol dependency)
5. **2D Visualization** — `plot_dual_profile()` renders upper (KD) / lower (HW) panels; `plot_comparison_scatter()` cross-plots both scales with region color-coding
6. **Streamlit UI** — `render_sidebar()` returns parameter dict; `render_analysis()` drives per-sequence tabs (2D Profile / KD vs HW / 3D Structure / Results); `main()` handles input routing (file upload / sample FASTA / paste)

## Key Design Decisions

- **PDB session cache**: predicted structures are stored in `st.session_state` keyed by `hash(sequence)` to avoid redundant API calls on re-render
- **KD vs HW polarity**: KD is positive-hydrophobic; HW is positive-hydrophilic. Both use `direction="high"` in `detect_regions()` for their respective feature types (TM regions for KD, antigenic sites for HW)
- **ESMFold limit**: sequences >400 aa are truncated to first 400 residues before API submission; the full sequence is still used for 2D profiling
- **Sample sequence**: β2-adrenergic receptor (UniProt P07550), a 7-TM GPCR, is hardcoded as `SAMPLE_FASTA` for demonstration
