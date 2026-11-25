# Annotation inputs

Included: curated direct transcriptional regulatory interactions (`dtri_exp.tsv`, `dtri_plos.tsv`),
mouse marker genes (`mancarci_mouseMarkerGenesNCBI.rds`), ribosomal protein genes, Eisenberg
housekeeping genes, and the IHC cell-type marker lists in `patrick_ihc/`.

Not included (large or public; download them into this folder to run the pipeline):
- `Generic_human_ensemblIds_bioProcess.an.txt`, `Generic_human_ncbiIds_bioProcess.an.txt`:
  Gemma annotation files for the generic human platform (Lim et al. 2021, Database baab006;
  https://gemma.msl.ubc.ca), generated 2022-06-30.
- `pro_GO.csv`: GO biological-process annotations of human genes.
- `lau_celltypes.tsv`: cell-type labels per cell barcode for the Lau et al. (2020) dataset.
