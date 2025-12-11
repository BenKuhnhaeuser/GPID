# GeneParliamentID: A pipeline for multi-gene DNA barcoding [UNDER CONSTRUCTION]

**GeneParliamentID** by Benedikt Kuhnhäuser, Royal Botanic Gardens Kew  
Current version: 1.0 (December 2025)

## Purpose
GeneParliamentID (GPID) is a pipeline for sample identification using hundreds or thousands of genes, such as those generated using targeted sequence capture.  

GeneParliamentID integrates species identifications inferred from individual genes to provide an overall identification that reflects the relative support for each alternative identification. We conceptualise this process as a “Gene Parliament” in which each gene represents one part of the genomic identity of an individual, and where the overall species identity is established through consideration of the number of genes supporting each different identification. This approach allows explicit assessment of congruence and discordance among multiple genes in species identification, comparable to the established concept of gene tree discordance in phylogenomics.  

The pipeline incorporates a sequence of several filters to increase the overall accuracy of identification by reducing noise introduced by low-confidence identifications, and is accompanied by a separate script for identifying the optimal filtering thresholds to facilitate implementation in any given dataset.

## Wiki
Our [wiki](https://github.com/BenKuhnhaeuser/GPID/wiki) contains detailed instructions on the following topics:
- Installation
- GeneParliamentID pipeline
  - Pipeline parameters
  - Input file formats
- Building a reference database
- Method calibration
  - Test data
  - Identify best-performing genes for identification
  - Select optimal filtering thresholds
  - Estimate confidence
- Interpretation of the Gene Parliament

## Setup
### Clone GeneParliamentID repository
Clone the GeneParliamentID repository from GitHub:  
`git clone https://github.com/BenKuhnhaeuser/GPID.git`

### Install dependencies
We recommend installation of all dependencies using [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) with a new environment:  
`conda create --name gpid`

Activate conda environment:  
`conda activate gpid`

Install BLAST and R:  
`conda install blast=2.17.0 r=4.5.0`  

Install R packages:  
`conda install r-optparse=1.7.5 r-dplyr=1.1.4 r-tidyr=1.3.1 r-strex=2.0.1 r-withr=3.0.2 r-ggplot2=4.0.0 r-rcolorbrewer=1.1_3 r-ggtext=0.1.2`


## Run GeneParliamentID pipeline
You can run the GeneParliamentID pipeline like this:
`Rscript gpid.r --input <file> [--performance <file>] [--thresholds <file>] [--confidence <file>] [--grouping <file>]`

**Required arguments:**  
`-i` / `--input`: File containing results of multi-gene BLAST search for unknown sample. [Example].  

**Optional arguments:**  
`-p` / `--performance`: Optional calibration file containing gene performances. [Example].  
`-t` / `--thresholds`: Optional calibration file containing filtering thresholds. If omitted, default thresholds are applied. [Example].  
`-c` / `--confidence`:Optional calibration file providing confidence estimates depending on gene support. If omitted, no confidence estimate is added. [Example].  
`-g` / `--grouping`: Optional file assigning species to species-groups. If omitted, no grouping information is implemented. [Example].  

All files have specific formatting requirements. Please see the [wiki](https://github.com/BenKuhnhaeuser/GPID/wiki) for instructions and check the example files linked above.

## Outputs
The pipeline produces three outputs:
1) Gene Parliament table listing proportion of genes for for each identification: `Sample1_gpid.csv`  
2) Gene Parliament graph visualisting proportion of genes for for each identification: `Sample1_gpid.pdf`  
3) Top identification table detailing top identification and confidence in the calibration: `Sample1_topid.csv`

For a detailed description of these outputs and their interpretation, please see the [wiki](https://github.com/BenKuhnhaeuser/GPID/wiki).  


## Citation
Kuhnhäuser, B.G., Quintero-Berns, C., Schley, R., Stevenson, J., Ndiade Bourobou, D., Cziba, L., Deklerck, V., Gallego, B., Lisingo, J., Baker, W.J. & Bellot, S. (2025). **GeneParliamentID: A pipeline for multi-gene DNA barcoding.** Submitted to Molecular Ecology Resources.
