
# GeneParliamentID:<br>A pipeline for multi-gene species identification

**GeneParliamentID** by Benedikt Kuhnhäuser, Royal Botanic Gardens, Kew  
Current version: 1.1.2 (June 2026)

**Citation**  
Kuhnhäuser, B.G., Quintero-Berns, C., Schley, R., Stevenson, J., Ndiade Bourobou, D., Cziba, L., Deklerck, V., Gallego, B., Lisingo, J., Baker, W.J. & Bellot, S. **GeneParliamentID: A pipeline for multi-gene species identification.** (In review)

## Overview
GeneParliamentID (GPID) is a pipeline for identification of biological samples to the species level using hundreds or thousands of genes, such as those generated using targeted sequence capture.  

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

The computational requirements for running GPID depend on the dataset. As a point of reference, running the tutorial with 1 CPU requires approx. 175 MB memory.

### Install GPID
We recommend installation of GPID including all dependencies using [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) with a new environment:  
`conda create --name gpid gpid`  

The dependencies installed are [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and [R](https://www.r-project.org/) with the packages [dplyr](https://dplyr.tidyverse.org), [ggplot2](https://ggplot2.tidyverse.org), [ggpubr](https://rpkgs.datanovia.com/ggpubr), [ggtext](https://github.com/wilkelab/ggtext), [stringr](https://stringr.tidyverse.org), [svglite](https://github.com/r-lib/svglite), [tidyr](https://tidyr.tidyverse.org) and [withr](https://withr.r-lib.org).  

### Activate GPID environment
To activate the GPID conda environment, use:  
`conda activate gpid`

### Confirm successful installation
To confirm that the installation has worked and show a help message on how to use GPID, simply run:  
`gpid`  

## Quick start
To give you a first taste of the capabilities of GPID, this is a minimal example for the identification of  that only covers sample identification.  
For this purpose, we identify a rattan palm of the genus *Calamus* using a small test dataset comprised of 96 genes and 91 reference taxa.  

[Reference construction](https://github.com/BenKuhnhaeuser/GPID/wiki/1.-Reference-construction), [Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration) and [Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation) have already been performed. Note that **these steps only need to be conducted once** for a given lineage and set of genes.  

For a full worked example including reference construction to method calibration and validation, see the [Tutorial](https://github.com/BenKuhnhaeuser/GPID/wiki/6.-Tutorial).  


### Download data
To download the folder, run:  
`wget https://github.com/BenKuhnhaeuser/GPID/blob/main/quickstart.tar.gz`  

Then, extract the files in the folder using:  
`tar -zxvf quickstart.tar.gz`  

The extracted directory contains the following files and folders:  
- `reference`: folder containing BLAST reference databases for all genes. See [1. Reference construction](https://github.com/BenKuhnhaeuser/GPID/wiki/1.-Reference-construction).
- `gene_performance.csv`: file listing performance of each gene, i.e. the percentage of samples correctly identified to species, estimated using method calibration. See [2. Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration).
- `thresholds.csv`: file with optimal filtering thresholds selected using method calibration. See [2. Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration).
- `confidence_support.csv`: file listing the probability of the top identification being correct, close or wrong depending on the percentage of genes supporting the identification; produced during method validation. See [3. Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation).
- `species_groups.csv`: optional file specifying for each species a user-defined group of closely related species. See See [3. Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation).
- `samples`: folder containing a sub-folder for each sample to be identified. Each sub-folder contains one fasta file per gene for the sample to be identified, and each file contains a single corresponding gene sequence for the sample. See [4. Sample identification](https://github.com/BenKuhnhaeuser/GPID/wiki/4.-Sample-identification).   

### Sample identification
Sample identification is conducted using the following command:  
`gpid identify -i samples/Calamus_sample_1/ -r reference/ -g gene_performance.csv -t thresholds.csv -c confidence_support.csv -s species_groups.csv`  

The required input arguments are:  
`-i`: Sample directory containing one `fasta` file per gene for the sample to identify  
`-r`:  Reference directory containing one `fasta` file per gene and the corresponding BLAST databases  

For method calibration and validation, the standard file names and locations produced during method calibration and validation are used by default. This can be specified using:  
`-g`: Gene performance file  
`-t`: Filtering thresholds file  
`-c`: Confidence estimates file  

Optionally, manually defined groups of closely related species can be included in the results:  
`-s`: Species groups file
 
See [4. Sample identification](https://github.com/BenKuhnhaeuser/GPID/wiki/4.-Sample-identification) for a full list of arguments and detailed instructions on the requirements for each argument.


### Pipeline outputs
GPID summarises all individual gene identifications in a Gene Parliament, which represents the percentage of genes supporting all competing identifications.  

The following output files are produced in the directory `identification/Calamus_sample_1`:
- Gene Parliament figure with top 10 identifications: `Calamus_sample_1_gpid.pdf`
- Gene Parliament table with all identifications: `Calamus_sample_1_gpid.csv`
- Individual BLAST identifications of each gene: `Calamus_sample_1_blast.tsv`
  
For a detailed description of the Gene Parliament figure and table and their interpretation, see [Interpretation](https://github.com/BenKuhnhaeuser/GPID/wiki/5.-Interpretation).

#### Gene Parliament figure
The Gene Parliament figure `Calamus_sample_1_gpid.pdf` gives a quick overview of the (up to) top 10 identifications that were retrieved and their relative support.  

<img width="1095" height="545" alt="image" src="https://github.com/user-attachments/assets/1fd4b993-5873-4477-9a32-4440ff866b73" />

In this case, a clear majority of genes support the identification as *Calamus melanochaetes*. There are several other species identifications, but these are only supported by one or two genes.  


#### Gene Parliament table
To see all identifications that were retrieved, we can have a look at the Gene Parliament table `Calamus_sample_1_gpid.csv`. This file contains all identifications in tabular format, which may be useful for further analysis of the results.  

Importantly, the table contains the results from [Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation), providing the percentage of validation samples that were correctly, closely (to species group) or wrongly identified at this level of support. In this case, we had decided to divide the validation samples into three bins of `0-33.33%`, `>33.33-66.66%` and `>66.66-100%`. Based on the percentage of genes supporting the identification of `Calamus_sample_1` as *Calamus melanochaetes*, the validation results from the third bin `>66.66-100%` are used. In this case, all validation samples at this level of support were correct identified. The variable `ID_correct_pct` thus contains the value `100`, whereas the variables `ID_close_pct` and `ID_wrong_pct` are `0`.  

| Sample | Rank | Identification | Species_group | Support_pct | Support_count | Parliament_size | Data_checks | ID_correct_pct | ID_close_pct |ID_wrong_pct|  
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| Calamus_sample_1 | 1 | Calamus_melanochaetes | Melanochaetes_group | 74.19 | 23 | 31 | PASSED | 100 | 0 | 0 |
| Calamus_sample_1 | 2 | Calamus_calicarpus | Melanochaetes_group | 6.45 | 2 |
| Calamus_sample_1 | 3 | Calamus_ater | Acamptostachys_group | 3.22 | 1 |
| Calamus_sample_1 | 3 | Calamus_mollispinus | Applanatus_group | 3.22 | 1 |
| Calamus_sample_1 | 3 | Calamus_ruber | Resiniferae_group | 3.22 | 1 |
| Calamus_sample_1 | 3 | Calamus_pseudoconcolor | Concolor_group | 3.22 | 1 |
| Calamus_sample_1 | 3 | Calamus_calapparius | Calapparius_group | 3.22 | 1 |
| Calamus_sample_1 | 3 | Calamus_subangulatus | Concolor_group | 3.22 | 1 |

### And the identification is...?
Based on these results, we can have high confidence that the sample is [*Calamus melanochaetes*](https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:665239-1), which is widespread across the Asian tropics.  

## Next steps
Have a go at identifying the other two samples, `Calamus_sample_2` and `Calamus_sample_3`. To check whether you are interpreting the results correctly, have a look at [Interpretation](https://github.com/BenKuhnhaeuser/GPID/wiki/5.-Interpretation).  

To work through the full GPID workflow including reference construction, method calibration and method validation, explore the [Tutorial](https://github.com/BenKuhnhaeuser/GPID/wiki/6.-Tutorial).
