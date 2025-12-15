# GeneParliamentID:<br>A pipeline for multi-gene DNA barcoding [UNDER CONSTRUCTION]


**GeneParliamentID** by Benedikt Kuhnhäuser, Royal Botanic Gardens, Kew  
Current version: 1.0 (December 2025)

## Overview
GeneParliamentID (GPID) is a pipeline for sample identification using hundreds or thousands of genes, such as those generated using targeted sequence capture.  

GeneParliamentID integrates species identifications inferred from individual genes to provide an overall identification that reflects the relative support for each alternative identification. We conceptualise this process as a **“Gene Parliament”** in which each gene represents one part of the genomic identity of an individual, and where the overall species identity is established through consideration of the number of genes supporting each different identification. This approach allows explicit assessment of congruence and discordance among multiple genes in species identification, comparable to the established concept of gene tree discordance in phylogenomics.  

The pipeline incorporates a sequence of several filtering steps to increase the accuracy of identification. It is accompanied by a calibration script that allows identifying the optimal filtering thresholds in any given dataset.

For detailed instructions on how to use, calibrate and interpret GeneParliamentID, please visit our [wiki](https://github.com/BenKuhnhaeuser/GPID/wiki). 

## Quick start
### Installation
Clone the repository:  
`git clone https://github.com/BenKuhnhaeuser/GPID.git`  

Install dependencies in new [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) environment:  
`conda create --name gpid`  
`conda activate gpid`  
`conda install blast=2.17.0 r=4.5.0 r-optparse=1.7.5 r-dplyr=1.1.4 r-tidyr=1.3.1 r-strex=2.0.1 r-withr=3.0.2 r-ggplot2=4.0.0 r-rcolorbrewer=1.1_3 r-ggtext=0.1.2`  

For more detailed instructions, see [Setup](https://github.com/BenKuhnhaeuser/GPID/wiki/Setup/).

### Run pipeline
Activate conda environment:  
`conda activate gpid`  

Run GeneParliamentID pipeline:  
`bash gpid.sh --sample <directory> --reference <directory> --performance <file> --thresholds <file> --confidence <file> [--grouping <file>]`

When you are done, deactivate the environment again using:  
`conda deactivate`  

## Pipeline inputs
**Required arguments:**  
Test sample and reference database:  
`-s` / `--sample`: Directory containing multipe genes for a sample of unknown identity  
`-r` / `--reference`: Directory with corresponding reference databases of each gene for lineage of interest  

Calibration files to set pipeline parameters:  
`-p` / `--performance`: Gene performance (percentage of correctly identified test samples for each gene)  
`-t` / `--thresholds`: Filtering thresholds  
`-c` / `--confidence`: Confidence estimates depending on gene support  

**Optional argument:**  
`-g` / `--groups`: User-defined groups of closely related species  
  
All files have specific formatting requirements. Please see the [Pipeline inputs](https://github.com/BenKuhnhaeuser/GPID/wiki/Pipeline-inputs) for detailed instructions and example files.

## Pipeline outputs
The pipeline produces three outputs:
1) Gene Parliament table listing proportion of genes for each identification: `<Sample1>_gpid.csv`  
2) Gene Parliament graph visualisting proportion of genes for each identification: `<Sample1>_gpid.pdf`  
3) Top identification, including confidence estimate: `<Sample1>_topid.csv`

For a detailed description of these outputs and their interpretation, please see the [wiki](https://github.com/BenKuhnhaeuser/GPID/wiki#outputs).


## Citation
Kuhnhäuser, B.G., Quintero-Berns, C., Schley, R., Stevenson, J., Ndiade Bourobou, D., Cziba, L., Deklerck, V., Gallego, B., Lisingo, J., Baker, W.J. & Bellot, S. (2025). **GeneParliamentID: A pipeline for multi-gene DNA barcoding.** Submitted to Molecular Ecology Resources.
