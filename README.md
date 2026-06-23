# Genomic representation gaps in Colombian bat diversity (Chiroptera)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20804719.svg)](https://doi.org/10.5281/zenodo.20804719)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0002--0450--5586-brightgreen)](https://orcid.org/0009-0002-0450-5586)

**Veritas Lab | Bogotá, Colombia**

Bioinformatic assessment of genomic sequence availability for Colombian bats in NCBI GenBank and BOLD Systems.

---

## Overview

Colombia harbors 224 bat species — approximately 15% of global Chiroptera richness. Despite this exceptional diversity, only 26% of Colombian bat species have any genomic representation in public databases. This repository contains all scripts used to quantify this gap, characterize available sequences by marker and taxonomy, and reconstruct maximum likelihood phylogenies from publicly available data.

**Key findings:**
- 445 sequences from NCBI GenBank representing 57 species
- 177 records from BOLD Systems representing 29 species
- 69 species total with genomic representation (25,7% of 222)
- 165 species with no sequence data in any public database
- COI barcode available for only 16 species in NCBI

---

## Repository structure
chiroptera-colombia-genomics/

├── 01_descarga_ncbi.py        # NCBI GenBank sequence retrieval via Biopython/Entrez

├── 02_descarga_bold.py        # BOLD Systems v5 API retrieval

├── 03_analisis_vacios.py      # Taxonomic gap analysis and visualization

├── 04_filogenia.py            # MAFFT alignment + IQ-TREE phylogenetic inference

├── 05_visualizacion.R         # Basic phylogenetic tree visualization (ggtree)

└── 06_arbol_publicable.R      # Publication-ready phylogenetic figures

---

## Dependencies

**Python 3.11+**

biopython>=1.84

pandas

matplotlib

seaborn

requests

Install with:
```bash
pip install biopython pandas matplotlib seaborn requests
```

**R 4.4+**
ggtree

ape

ggplot2

dplyr

treeio

Install with:
```r
BiocManager::install(c("ggtree", "treeio"))
install.packages(c("ape", "ggplot2", "dplyr"))
```

**External tools**
- [MAFFT v7.520](https://mafft.cbrc.jp/alignment/software/)
- [IQ-TREE 2](http://www.iqtree.org/)

Install on macOS:
```bash
brew install mafft iqtree
```

---

## Usage

### 1. Retrieve sequences from NCBI
```bash
python3 01_descarga_ncbi.py
```
Output: `data/raw/chiroptera_colombia_raw.gb` and `data/raw/chiroptera_colombia_metadata.csv`

### 2. Retrieve barcodes from BOLD
```bash
python3 02_descarga_bold.py
```
Output: `data/raw/bold_chiroptera_colombia_raw.tsv`

### 3. Gap analysis
```bash
python3 03_analisis_vacios.py
```
Output: summary tables in `results/tables/` and figures in `results/figures/`

### 4. Phylogenetic inference
```bash
python3 04_filogenia.py
```
Output: aligned FASTA files and `.treefile` in `results/`

### 5. Visualization in R
Open `06_arbol_publicable.R` in RStudio and run.

---

## Data availability

Raw metadata tables and alignment files are deposited at Zenodo:
> Melo Escobar, D.A. (2026). Genomic representation gaps in Colombian bat diversity — dataset. *Zenodo*. https://doi.org/10.5281/zenodo.20804719

Sequences analyzed in this study are publicly available at:
- NCBI GenBank: https://www.ncbi.nlm.nih.gov
- BOLD Systems: https://portal.boldsystems.org

---

## Citation

If you use these scripts or data, please cite:

> Melo Escobar, D.A. (2026). Only 31% of Colombian bat species have genomic representation: a bioinformatic assessment of Chiroptera sequences in public databases. *PeerJ*. [DOI pending]

---

## Author

**Diego Alfonso Melo Escobar**
Veritas Lab, Bogotá, Colombia
ORCID: [0009-0002-0450-5586](https://orcid.org/0009-0002-0450-5586)

---

*Veritas Lab is an independent computational biology and science communication initiative based in Bogotá, Colombia.*
