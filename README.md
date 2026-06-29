
# GeneParliamentID:<br>A pipeline for multi-gene DNA barcoding

**GeneParliamentID** by Benedikt Kuhnhäuser, Royal Botanic Gardens, Kew  
Current version: 1.1.2 (June 2026)

**Citation**  
Kuhnhäuser, B.G., Quintero-Berns, C., Schley, R., Stevenson, J., Ndiade Bourobou, D., Cziba, L., Deklerck, V., Gallego, B., Lisingo, J., Baker, W.J. & Bellot, S. **GeneParliamentID: A pipeline for multi-gene DNA barcoding.** (In review)

## Overview
GeneParliamentID (GPID) is a pipeline for sample identification using hundreds or thousands of genes, such as those generated using targeted sequence capture.  

GPID integrates species identifications inferred from individual genes to provide an overall identification that reflects the relative support for each alternative identification. We conceptualise this process as a **“Gene Parliament”** in which each gene represents one part of the genomic identity of an individual, and where the overall species identity is established through consideration of the number of genes supporting each different identification. This approach allows explicit assessment of congruence and discordance among multiple genes in species identification. Besides sample identification, the pipeline includes reference directory preparation, method calibration and method validation.    

The pipeline is structured into four main commands that are explained in detail in the [Wiki](https://github.com/BenKuhnhaeuser/GPID/wiki):  
1. `gpid reference`: Prepare a reference directory. See [1. Reference construction](https://github.com/BenKuhnhaeuser/GPID/wiki/1.-Reference-construction).  
2. `gpid calibrate`: Run the calibration workflow to identify the optimal pipeline settings. See [2. Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration).  
3. `gpid validate`: Run validation analyses on samples with known identity to test the accuracy of identification. See [3. Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation).  
4. `gpid identify`: Run the identification workflow for sample identification using optimal pipeline settings. See [4. Sample identification](https://github.com/BenKuhnhaeuser/GPID/wiki/4.-Sample-identification).

The Wiki also contains guidance on the [Interpretation](https://github.com/BenKuhnhaeuser/GPID/wiki/5.-Interpretation) of the identification results and a hands-on [Tutorial](https://github.com/BenKuhnhaeuser/GPID/wiki/6.-Tutorial) using example data.


## Pipeline summary
The key steps of the GeneParliamentID pipeline are:
1. Select genes with performance of species identification above user-defined threshold
2. **Align each gene** against reference, select top match based on highest Bit-score
3. Remove low-confidence matches that don't meet the user-defined alignment filtering thresholds
4. Summarise the number of genes supporting each identification to the **Gene Parliament**
5. Flag the sample as data-deficient if the number of genes (parliament size) is below the user-defined threshold
6. Select identification with most support as the **top identification**
7. Evaluate confidence in the top identification based on the percentage of genes supporting this identification  
<img width="600" alt="GeneParliamentID pipeline" src="https://github.com/user-attachments/assets/2fe0f86a-e1a3-4603-b3f5-0e6623bd176b" />

## Setup

### System requirements
GeneParliamentID is a command-line tool written for Unix-operating systems such as Linux.  

The minimum requirements for running GPID are 1 CPU and 1 GB memory. Depending on the size of the dataset analysed, more processing power or memory may be needed.

### Install GPID
We recommend installation of GPID including all dependencies using [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) with a new environment:  
`conda create --name gpid gpid`

### Activate GPID environment
To activate the GPID conda environment, use:  
`conda activate gpid`

### Confirm successful installation
To confirm that the installation has worked and show a help message on how to use GPID, simply run:  
`gpid`  

## Quick start
To give you a first taste of the capabilities of GPID, this is a minimal example that only covers sample identification.  

[Reference construction](https://github.com/BenKuhnhaeuser/GPID/wiki/1.-Reference-construction), [Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration) and [Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation) have already been performed. Note that **these steps only need to be conducted once** for a given lineage and set of genes.  

For a full worked example including reference construction to method calibration and validation, see the [Tutorial](https://github.com/BenKuhnhaeuser/GPID/wiki/6.-Tutorial).  


### Download data
To download the folder, run:  
`wget https://github.com/BenKuhnhaeuser/GPID/blob/main/quickstart.tar.gz`  

Then, extract the files in the folder using:  
`tar -zxvf example_data.tar.gz`  

The extracted directory contains the following files and folders:  
- `reference`: folder containing BLAST reference databases for all genes. See [1. Reference construction](https://github.com/BenKuhnhaeuser/GPID/wiki/1.-Reference-construction).
- `calibration_gene_performance.csv`: file listing performance of each gene, i.e. the percentage of samples correctly identified to species, estimated using method calibration. See [2. Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration).
- `calibration_filtering_thresholds.csv`: file with optimal filtering thresholds selected using method calibration. See [2. Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration).
- `validation_confidence_support.csv`: file listing the probability of the top identification being correct, close or wrong depending on the percentage of genes supporting the identification; produced during method validation. See [3. Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation).
- `species_groups.csv`: optional file specifying for each species a user-defined group of closely related species. See See [3. Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation).
- `samples`: folder containig a sub-folder for each sample to be identified. Each sub-folder contains one fasta file per gene for the sample to be identified, and each file contains a single corresponding gene sequence for the sample. See [4. Sample identification](https://github.com/BenKuhnhaeuser/GPID/wiki/4.-Sample-identification).   

### Sample identification
Sample identification is conducted using the following command:  
`gpid identify -i samples/CQL_2/ -r reference/ -g calibration_gene_performance.csv -t calibration_filtering_thresholds.csv -c calibration_confidence_support.csv -s species_groups.csv`  

The required input arguments are:  
`-i`: Sample directory containing one FASTA file per gene for the sample to identify  
`-r`:  Reference directory containing one FASTA file per gene and the corresponding BLAST databases  

For method calibration and validation, the standard file names and locations produced during method calibration and validation are used by default. This can be specified using:  
`-g`: Gene performance file  
`-t`: Filtering thresholds file  
`-c`: Confidence estimates file  

Optionally, manually defined groups of closely related species can be included in the results:  
`-s`: Species groups file
 
See [4. Sample identification](https://github.com/BenKuhnhaeuser/GPID/wiki/4.-Sample-identification) for a full list of arguments and detailed instructions on the requirements for each argument.


### Pipeline outputs
GPID summarises all individual gene identifications in a Gene Parliament, which represents the percentage of genes supporting all competing identifications. The Gene Parliament is presented both as a table and as a figure.  

For a detailed description of the Gene Parliament figure and table and their interpretation, see [Interpretation](https://github.com/BenKuhnhaeuser/GPID/wiki/5.-Interpretation).

#### Gene Parliament figure
The Gene Parliament figure gives a quick overview of the top 10 identifications that were retrieved and their relative support.  

<img width="800" alt="Gene Parliament CQL_2" src="https://github.com/user-attachments/assets/d3b1d9da-e18f-4507-bfb4-8ee523509d68" />

#### Gene Parliament table
To see all identifications that were retrieved, we can have a look at the Gene Parliament table. Importantly, the table contains for the top identification information on the probability of the identification being correct (correct to species level), close (correct to species group) or wrong (neither correct to species nor to species group).  

| Sample | Rank | Identification | Species_group | Support_pct | Support_count | Parliament_size | Data_checks | ID_correct_pct | ID_close_pct |ID_wrong_pct|  
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| CQL_2 | 1 | *Entandrophragma_angolense* | Entandrophragma | 45.08 | 55| 122 | PASSED | 82.14 | 17.86 | 0 |
| CQL_2 | 2 | *Entandrophragma_excelsum* | Entandrophragma | 18.85 | 23 |
| CQL_2 | 3 | *Entandrophragma_congoense* | Entandrophragma | 16.39 | 20 |
| ... |  |  |  |  |  |

### And the identification is...?
Inspection of the Gene Parliament figure and table shows that a clear majority of genes (45.1%) support the identification as *Entandrophragma angolense*, whilst *Entandrophragma excelsum* (18.85%) and *Entandrophragma congoense* (16.39%) also get sizeable support. Other species have almost negligible support, but all belong to the same genus *Entandrophragma*.  

The table indicates a probability of 82.14% that the identification is `correct` to species level, and a further 17.86% (thus totaling 100%) that the identification is `close`, i.e. correct to species group (in this case genus). The probability that the identification is `wrong`, i.e. neither correct nor close, is estimated to be 0. Overall, we can be fully confident that the sample belongs to the genus *Entandrophragma*, and have high confidence that it was taken from the species *Entandrophragma_angolense*.


## Next steps
To work through the full GPID workflow using example data, explore the [Tutorial](https://github.com/BenKuhnhaeuser/GPID/wiki/6.-Tutorial).
