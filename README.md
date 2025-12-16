# GeneParliamentID:<br>A pipeline for multi-gene DNA barcoding [UNDER CONSTRUCTION]


**GeneParliamentID** by Benedikt Kuhnhäuser, Royal Botanic Gardens, Kew  
Current version: 1.0 (December 2025)

## Overview
GeneParliamentID (GPID) is a pipeline for sample identification using hundreds or thousands of genes, such as those generated using targeted sequence capture.  

GeneParliamentID integrates species identifications inferred from individual genes to provide an overall identification that reflects the relative support for each alternative identification. We conceptualise this process as a **“Gene Parliament”** in which each gene represents one part of the genomic identity of an individual, and where the overall species identity is established through consideration of the number of genes supporting each different identification. This approach allows explicit assessment of congruence and discordance among multiple genes in species identification, comparable to the established concept of gene tree discordance in phylogenomics.  

The pipeline incorporates a sequence of several filtering steps to increase the accuracy of identification. It is accompanied by a calibration script that allows identifying the optimal filtering thresholds in any given dataset.

For detailed instructions on how to use, calibrate and interpret GeneParliamentID, please visit our [wiki](https://github.com/BenKuhnhaeuser/GPID/wiki). 

## Quick start
The fundamental commands for running GeneParliamentID are given here. For a worked example using test data, see the [Tutorial](https://github.com/BenKuhnhaeuser/GPID/wiki/Tutorial).

### Installation
Clone the repository:  
`git clone https://github.com/BenKuhnhaeuser/GPID.git`  

Install dependencies in new [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) environment:  
`conda create --name gpid`  
`conda activate gpid`  
`conda install blast=2.17.0 r=4.5.0 r-optparse=1.7.5 r-dplyr=1.1.4 r-tidyr=1.3.1 r-strex=2.0.1 r-withr=3.0.2 r-ggplot2=4.0.0 r-rcolorbrewer=1.1_3 r-ggtext=0.1.2`  

For more detailed instructions, see [Setup](https://github.com/BenKuhnhaeuser/GPID/wiki/Setup/).

### Export path to directory containing scripts  
`export PATH=$PATH:</path/to/gpid/scripts>`  
This allows to execute the scripts from any location.  
  
Notes:  
- Change `</path/to/gpid/scripts>` to the actual path on your machine
- This path will only be valid for the length of your session
- Remember to omit the angle brackets, e.g. $PATH:/path/to/gpid/scripts
- Do not include the script itself in the path

### Run pipeline
Activate conda environment:  
`conda activate gpid`  

Run GeneParliamentID pipeline:  
`bash gpid.sh --sample <directory> --reference <directory> --performance <file> --thresholds <file> --confidence <file> [--groups <file>]`

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
  
See [Pipeline parameters](https://github.com/BenKuhnhaeuser/GPID/wiki/Pipeline-parameters) for detailed instructions on the requirements for each argument.

## Pipeline outputs
The GeneParliamentID pipeline summarises all individual gene identifications in a Gene Parliament, which represents the percentage of genes supporting all competing identifications.  
  
The Gene Parliament is presented both as a table `<Sample>_gpid.csv` and as a figure `<Sample>_gpid.pdf`.
  
For a detailed description of these outputs and their interpretation, see [Interpretation](https://github.com/BenKuhnhaeuser/GPID/wiki/Interpretation).

## Method calibration
When running the GeneParliamentID pipeline for the first time for a new lineage, method calibration using a test dataset is required to identify the optimal pipeline parameters for this lineage. This is not necessary if optimal pipeline parameters have already been established for the lineage of interest. 

For details, see [Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/Method-calibration).

## Citation
Kuhnhäuser, B.G., Quintero-Berns, C., Schley, R., Stevenson, J., Ndiade Bourobou, D., Cziba, L., Deklerck, V., Gallego, B., Lisingo, J., Baker, W.J. & Bellot, S. (2025). **GeneParliamentID: A pipeline for multi-gene DNA barcoding.** Under review.
