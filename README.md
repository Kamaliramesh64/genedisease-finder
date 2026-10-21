# GeneDisease Finder

Identification of Disease Associations from Selected Gene Features Through Structural and Functional Annotation and Pathway Mapping.

Paste or upload a gene list (text, CSV, TXT or Excel) and get:

1. **Gene features**: structure (location, exons, protein length, domains, AlphaFold link, ClinVar variant counts) and function (GO terms, top tissues, druggability)
2. **Pathways**: enriched Reactome pathways, shared features, interaction network
3. **Diseases**: ranked disease associations, evidence-score filter, related diseases, drug targets
4. **Compare**: two gene sets side by side, charts, CSV downloads, PDF report

## Data sources
MyGene.info, UniProt, AlphaFold DB, Enrichr (Reactome), STRING, Open Targets, MyVariant.info (ClinVar). All free; an internet connection is required.

## Run locally
    pip install -r requirements.txt
    python app.py
Open http://127.0.0.1:5000 and click "Load example genes".

Department of Biomedical Engineering, SIMATS Engineering, Saveetha University, Chennai.
