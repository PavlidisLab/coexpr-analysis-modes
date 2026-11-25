# Original analysis pipeline (R)

The code that produced the coexpression networks, the results and most figures of the submitted
manuscript (written 2023-2024).

| Script | What it does |
|---|---|
| `041_SETUP.R` | workspace setup and shared functions (expression QC, coexpression, rank normalisation) |
| `041_MAIN_01_exprmats.R` | loads the raw single-cell and bulk data; QC and normalisation |
| `041_MAIN_02_coexmats.R` | xCell, xSubject and xBulk coexpression matrices |
| `041_MAIN_00_misc.R` | gene annotation |
| `041_MAIN_03_go_integrity.R` | GO AUROCs of the networks (EGAD) |
| `041_ANALYSIS_01` ... `07` | one script per results section (05 clusters the xSubject networks and draws no figure); the others write their panels to `figures/pipeline_panels/` |
| `041_supplementary_data_generation.R` | supplementary data files (coexpression matrices and tables) |
| `041_SANDBOX.R` | exploratory code; also draws `figures_sandbox/schematic_01.png` |
| `misc/` | annotation inputs: GO annotations, curated regulatory interactions (DTRI), marker genes, ribosomal protein genes, housekeeping genes, Lau cell-type labels |

The scripts were run on the lab server. They expect a shared utility library (`GlobalSession.R`
and `utils/`) that is not included here, and they read and write cached R objects under
server paths set in `041_SETUP.R`.

`code/revision_analyses/lib/eric_engine.py` re-implements the xCell and xSubject definitions
from `041_SETUP.R` in Python, and reproduces the saved networks and every Fisher odds ratio.
