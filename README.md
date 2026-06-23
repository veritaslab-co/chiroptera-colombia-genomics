# Genomic representation gaps in Colombian bat diversity (Chiroptera)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20804719.svg)](https://doi.org/10.5281/zenodo.20804719)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0002--0450--5586-brightgreen)](https://orcid.org/0009-0002-0450-5586)

**Veritas Lab | Bogotá, Colombia**

Bioinformatic assessment of genomic sequence availability for Colombian bats in NCBI GenBank and BOLD Systems, with comparative regional analysis across eight Neotropical countries.

---

## Overview

Colombia harbors 222 bat species — approximately 16% of global Chiroptera richness — including nine endemic species. Despite this exceptional diversity, only 25% of Colombian bat species have any genomic representation in public databases, and none of the nine endemic species has a single publicly available nucleotide sequence. This repository contains all scripts used to quantify this gap, characterize available sequences by marker, taxonomy, and geographic provenance, and reconstruct maximum likelihood phylogenies from publicly available data.

**Key findings:**

- 445 sequences from NCBI GenBank representing **56 valid species**
- 177 records from BOLD Systems representing 29 species
- **56 species** total with genomic representation (**25.2% of 222**)
- **166 species** with no sequence data in any public database
- COI barcode available for only **16 species** in NCBI (COI coverage index: **0.54** sequences/species — second-lowest in the Neotropics)
- **9 endemic bat species**: none has any publicly available sequence data
- **4 families entirely absent** from both databases: Emballonuridae, Mormoopidae, Noctilionidae, Furipteridae
- 42% of NCBI records retrieved carry non-Colombian collection localities (provenance bias)
- 39% of confirmed Colombian records originate from a single department (Córdoba)

---

## Repository structure

```
chiroptera-colombia-genomics/
├── scripts/
│   ├── 01_descarga_ncbi.py        # NCBI GenBank sequence retrieval via Biopython/Entrez
│   ├── 02_descarga_bold.py        # BOLD Systems v5 API retrieval
│   ├── 03_analisis_vacios.py      # Taxonomic gap analysis and visualization (Figs 1-3, 7)
│   ├── 04_filogenia.py            # MAFFT alignment + IQ-TREE phylogenetic inference
│   ├── 05_visualizacion.R         # Basic phylogenetic tree visualization (ggtree)
│   └── 06_arbol_publicable.R      # Publication-ready phylogenetic figures (Figs 4-5)
├── results/
│   ├── figures/                   # All manuscript figures (PNG, 300 dpi)
│   └── tables/                    # Coverage matrices, gap tables, geographic analysis
└── data/
    └── raw/                       # Raw downloads from NCBI and BOLD (see Zenodo)
```

---

## Dependencies

**Python 3.11+**

```
biopython>=1.84
pandas
matplotlib
seaborn
requests
```

Install with:

```bash
pip install biopython pandas matplotlib seaborn requests
```

**R 4.4+**

```
ggtree
ape
ggplot2
dplyr
treeio
```

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
python3 scripts/01_descarga_ncbi.py
```

Output: `data/raw/chiroptera_colombia_raw.gb` and `data/raw/chiroptera_colombia_metadata.csv`

### 2. Retrieve barcodes from BOLD

```bash
python3 scripts/02_descarga_bold.py
```

Output: `data/raw/bold_chiroptera_colombia_raw.tsv`

### 3. Gap analysis and figures

```bash
python3 scripts/03_analisis_vacios.py
```

Output: summary tables in `results/tables/` and figures in `results/figures/`

### 4. Phylogenetic inference

```bash
python3 scripts/04_filogenia.py
```

Output: aligned FASTA files and `.treefile` in `results/`

### 5. Publication-ready phylogenetic figures

Open `scripts/06_arbol_publicable.R` in RStudio and run section by section.

---

## Methods summary

- **Species checklist**: Ramírez-Chaves et al. (2024), Mamíferos de Colombia v1.14 (SCMas/SiB Colombia) — 222 species
- **Sequence retrieval**: NCBI GenBank (`Chiroptera[Organism] AND Colombia[Country]`); BOLD API v5
- **Quality filtering**: excluded taxonomic contaminants (*Plecotus austriacus* OX031246.1; *Rhinolophus ferrumequinum*); morphospecies and subspecies excluded from species-level counts
- **Geographic provenance**: `geo_loc_name` field extracted and parsed to department level
- **Alignment**: MAFFT v7.520 (`--auto`)
- **Phylogenetic inference**: IQ-TREE 2 with ModelFinder (BIC); best-fit models GTR+F+R4 (cytb) and TIM2+F+I+G4 (COI)
- **Outgroup**: *Pteropus vampyrus* (Yinpterochiroptera hypothesis; Teeling et al. 2005)
- **Visualization**: ggtree v3.10 in R 4.4

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

> Melo Escobar, D.A. (2026). Only 25% of Colombian bat species have genomic representation: a bioinformatic assessment of Chiroptera sequences in public databases. *PeerJ*. [DOI pending]

---

## Author

**Diego Alfonso Melo Escobar**
Independent Researcher | Bogotá, Colombia
ORCID: [0009-0002-0450-5586](https://orcid.org/0009-0002-0450-5586)
Email: veritaslab.latam@gmail.com

---

## License

MIT License — see LICENSE file for details.
