![R Version](https://img.shields.io/badge/R-%3E%3D%204.3-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Reproducibility](https://img.shields.io/badge/reproducibility-renv-orange)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.18677425-blue.svg)](https://doi.org/10.5281/zenodo.18677425)

# Netris LapNet
### Transcriptomic Analysis of PDAC Tissues (LapNet Study)

----

This repository contains the R code and analysis pipeline for the manuscript: 

*Roth et al. (2026). Netrin-1 blockade alleviates resistance to first line chemotherapy in locally advanced pancreatic cancer. Nature.*

The analysis covers the preprocessing of microBulk RNA-seq from FFPE samples, molecular subtyping, differential expression analysis (Pre vs. Post therapy), and integration with public datasets (GSE253260, GSE131050).

----

### 1. Project Overview
- Preprocessing: Alignment with Rsubread and quantification for FFPE-specific considerations.

- Classification: Molecular subtyping using PDACMOC, PurIST, and ESTIMATE.

- Differential Expression: edgeR pipeline for paired and unpaired analysis.

- Validation: Integration and GSEA comparison with Nicolle et al. and Linehan et al. datasets.


### 2. Getting Started
#### Prerequisites:

You will need R (>= 4.0) and the following core Bioconductor packages:
edgeR, Rsubread, clusterProfiler, org.Hs.eg.db, and singscore.


#### Installation:

Clone this repository:

```bash
git clone https://github.com/Genomics-Consulting/Netris_LapNet.git
```

Open `Netris_LapNet.Rproj` in RStudio.

Restore the R environment using `renv`:

```r
install.packages("renv")
renv::restore() 
```


#### Molecular Subtyping Setup:

1. First, create the base environment containing the R and Python dependencies:

```bash
# bash
conda env create -f env/PDACMOC_env.yml
conda activate PDACMOC
```
2. The PDACMOC and ADVOCATE packages must be installed from source after activating the environment. 
From within the directory where you have cloned the PDACMOC repository, run:

```r
# Inside an R session
# Install the ADVOCATE dependency required for stroma classification
install.packages("inst/packages/ADVOCATE_0.1.0.tar.gz", repos = NULL, type = "source")

# Install the main PDACMOC package
install.packages(".", repos = NULL, type = "source")
```
3. The PDACMOC package uses Python-based machine learning models. 
Script 2_molecular_classification.R is configured to use the Python interpreter provided by the Conda environment:

```r
reticulate::use_condaenv("PDACMOC")
```


#### External Data:

While data/pdata_LapNet.csv is included in this repo, the raw expression data and large objects should be downloaded from GEO:

- LapNet Data: GEO Accession [GSE319924](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319924)

- Public Data: GSE253260 (Nicolle et al.), GSE131050 (Linehan et al.), and GSE225691 (Cassier et al.).


### 3. Usage
The analysis is divided into sequential scripts. Please run them in order:

1. `scripts/1_preprocessing.R`: Alignment and initial QC.

2. `scripts/2_molecular_classification.R`: Application of PDAC-specific classifiers.

3. `scripts/3_differential_Pro.vs.Pre.R`: Identification of DEGs and GSEA.

4. `scripts/4_public_datasets.R`: Processing and integration of validation cohorts.

5. `scripts/5_LapNet_pdata.R`: Consolidation of metadata for final figures.

6. `scripts/6_Figures.R`: Generation of all main and extended figures.


### 4. Citation
If you use this code or the LapNet dataset, please cite:

*Roth et al. (2026). Netrin-1 blockade alleviates resistance to first line chemotherapy in locally advanced pancreatic cancer. Nature.*

----

> #### 📧 Contact & Collaboration
> For questions or collaboration opportunities, please contact us at **contact@genomicsconsulting.eu** or visit our [website](https://www.genomicsconsulting.eu/).
>
> **Genomics Consulting SARL** – 999 261 738 R.C.S. Lyon  
> Registered Office: 36B rue de la Batterie, 69500 Bron, France

----