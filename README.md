# Dosage sensitivity is a major determinant of human copy number variant pathogenicity

## Alan M. Rice and Aoife McLysaght

This repository contains code and data used for analyses in the manuscript.

### Code tested on

* Mac OS X - *10.11*
* Python - *v2.7.11*
* Bedtools - *v2.23.0*
* Bedops - *v2.4.12*
* R - *v3.2.3*

Python dependancies:
* pydash

### To run code and reproduce:

1. Install all required dependencies (Bedtools, Bedops, R)

2. Download code

  ```sh
  git clone git@github.com:alanrice/paper-dosage-sensitivity-copy-number-variation.git
  cd paper-dosage-sensitivity-copy-number-variation
  ```

3. Download archived datasets and copy to ```paper-dosage-sensitivity-copy-number-variation/analysis/datasets``` folder

  * [dbVar GRCh37 human remapped CNVs dated 31st Oct 2013](https://drive.google.com/open?id=0B856ApoNIDDPb09DMi14cnBtSkE)
  * [dbVar GRCh37 human submitted CNVs dated 31st Oct 2013](https://drive.google.com/open?id=0B856ApoNIDDPUF9HN1poa1hZbTQ)
  * [Ensembl GRCh37 homology data in JSON format](https://drive.google.com/open?id=0B856ApoNIDDPNFhoX01YU2N1d1E)
  * [Ensembl GRCh37 gene trees](https://drive.google.com/open?id=0B856ApoNIDDPeE4yVHc4c3NXQWs)
  * [Ensembl GRCh37 gene gain/loss trees](https://drive.google.com/open?id=0B856ApoNIDDPLTN3a2x0RDFZUnM)

4. Run code

  ```sh
  cd analysis
  chmod +x main.sh
  ./main.sh
  ```

### Important files to note

```
.
└── analysis
    ├── copyNumberAnalysis.py - Script for counting gene duplication/loss events in mammalian genomes
    ├── dbVarFilter.py - Script to filter dbVar database files
    ├── main.sh - Main bash script for running analysis and generating output files in files/
    └── files/
        └── speciesCopyNumberAnalysis.sorted.txt - Mammalian gene duplication/loss counts (Num unchanged - duplicated - no ortholog)
```

### Licence

Code - MIT

### Citation

Rice, A. M. & McLysaght, A. Dosage sensitivity is a major determinant of human copy number variant pathogenicity. Nat. Commun. 8, 14366 doi: 10.1038/ncomms14366 (2017).
