# GeneParliamentID:<br>A pipeline for multi-gene DNA barcoding

**GeneParliamentID** by Benedikt Kuhnhäuser, Royal Botanic Gardens, Kew  
Current version: 1.1.2 (June 2026)

**Citation**  
Kuhnhäuser, B.G., Quintero-Berns, C., Schley, R., Stevenson, J., Ndiade Bourobou, D., Cziba, L., Deklerck, V., Gallego, B., Lisingo, J., Baker, W.J. & Bellot, S. **GeneParliamentID: A pipeline for multi-gene DNA barcoding.** (In review)

## Overview
GeneParliamentID (GPID) is a pipeline for sample identification using hundreds or thousands of genes, such as those generated using targeted sequence capture.  

GPID integrates species identifications inferred from individual genes to provide an overall identification that reflects the relative support for each alternative identification. We conceptualise this process as a **“Gene Parliament”** in which each gene represents one part of the genomic identity of an individual, and where the overall species identity is established through consideration of the number of genes supporting each different identification. This approach allows explicit assessment of congruence and discordance among multiple genes in species identification. Besides sample identification, the pipeline includes reference directory preparation, method calibration and method validation.    

The pipeline is structured into four main commands:  
1. `gpid reference`: Prepare a reference directory  
2. `gpid calibrate`: Run the calibration workflow to identify the optimal pipeline settings  
3. `gpid validate`: Run validation analyses on samples with known identity to test the accuracy of identification
4. `gpid identify`: Run the identification workflow for sample identification using optimal pipeline settings  

## Wiki
For detailed instructions on how to use and interpret GPID, please visit the [Wiki](https://github.com/BenKuhnhaeuser/GPID/wiki).  

The Wiki contains advice on the following topics, including a separate page for each main `gpid` command:  
1. [Reference construction](https://github.com/BenKuhnhaeuser/GPID/wiki/1.-Reference-construction): Preparing a reference directory using `gpid reference`.  
2. [Method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration): Optimising the pipeline settings using `gpid calibrate`.  
3. [Method validation](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation): Estimating the accuracy of identification using `gpid validate`.    
4. [Sample identification](https://github.com/BenKuhnhaeuser/GPID/wiki/4.-Sample-identification). Identifying a sample using `gpid identify`.  
5. [Interpretation](https://github.com/BenKuhnhaeuser/GPID/wiki/5.-Interpretation): Interpretation of the identification results.  
6. [Tutorial](https://github.com/BenKuhnhaeuser/GPID/wiki/6.-Tutorial): Hands-on tutorial using example data.    

## Pipeline summary
The key steps of the GeneParliamentID pipeline are:
1. Select genes with performance of species identification above user-defined threshold
2. **Align each gene** against reference, select top match based on highest Bit-score
3. Remove low-confidence matches that don't meet the user-defined alignment filtering thresholds
4. Summarise the number of genes supporting each identification to the **Gene Parliament**
5. Flag the sample as data-deficient if the number of genes (parliament size) is below the user-defined threshold
6. Select identification with most support as the **top identification**
7. Evaluate confidence in the top identification based on the percentage of genes supporting this identification  
<img width="600" alt="GeneParliamentID pipeline" src="https://github.com/user-attachments/assets/57cf5a9f-b1e0-46da-85ba-d020ff2c93a6" />  
  

## Setup

### Operating system
GeneParliamentID is a command-line tool written for Unix operating systems such as Linux.  

### System requirements
The system requirements for running GPID depend on the size of the reference dataset and the number of genes analysed. For the datasets with which GPID was tested for the publication, they were:  
- Rattans: XXXGB RAM, XXX CPU cores
- Mahoganies: XXXGB RAM, XXX CPU cores 

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
The fundamental commands for running GeneParliamentID are given here. For detailed instructions, see the [wiki](https://github.com/BenKuhnhaeuser/GPID/wiki).  

### 1. Reference construction
Run a few checks and prepare a [reference](https://github.com/BenKuhnhaeuser/GPID/wiki/1.-Reference-construction) directory with BLAST databases for all genes in the directory. There is one single command:  
`gpid reference -r <reference directory>`

### 2. Method calibration
When running the GeneParliamentID pipeline for a lineage for the first time, method calibration using a test dataset is highly recommended to identify the optimal pipeline parameters for this lineage, which will increase the accuracy of identification. This involves multiple steps that are explained in detail [here](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration). The base command is:  
`gpid calibrate`  

Method calibration is structured into five subcommands:  
`gpid calibrate prepare`: Prepare input files for calibration by matching each test sample against the reference dataset  
`gpid calibrate alignments`: Identify optimal alignment filtering thresholds  
`gpid calibrate genes`: Estimate gene performance and identify optimal gene threshold  
`gpid calibrate parliament`: Identify optimal minimum parliament size threshold  
`gpid calibrate combine`: Combine manually selected thresholds in a calibration file for subsequent use  
  
Note: It is possible to [bypass method calibration](https://github.com/BenKuhnhaeuser/GPID/wiki/2.-Method-calibration#bypassing-method-calibration). This is *not* recommended as it will most likely result in considerably reduced accuracy of identification compared to running the pipeline with optimal parameters. However, it may be justified e.g. for a first explorative analysis or if a test dataset for calibration is not available.

### 3. Method validation
This allows to assess the accuracy of identification using test samples of known identity. Method validation consists of three subcommands that are explained in detail [here](https://github.com/BenKuhnhaeuser/GPID/wiki/3.-Method-validation). The base command is:  
`gpid validate`  

Method validation is structured into three subcommands:  
`gpid validate prepare`: Prepare input files for calibration by matching each test sample against the reference dataset  
`gpid validate confidence`: Estimate validation confidence for different numbers of support bins  
`gpid validate bins`: Save validation confidence support probabilities for selected number of bins  

### 4. Sample identification
Once the reference has been constructed and method calibration and validation are completed, you can conduct sample identification without the need to repeat the above steps. Sample identification is conducted using a single command that takes the outputs prepared using `gpid reference`, `gpid calibrate` and `gpid validate` as inputs. There is a single command with multiple flags:  
`gpid identify -i <sample directory> -r <reference directory> [-g <gene performance file> -t <thresholds file> -c <confidence support file>] [-s <species groups file>]`  

The required input arguments are:  
`-i`: Sample directory containing one FASTA file per gene for the sample to identify  
`-r`:  Reference directory containing one FASTA file per gene and the corresponding BLAST databases  

For method calibration and validation, the standard file names and locations produced in steps 2-3 are used by default. This can be specified using:  
`-g`: Gene performance file  
`-t`: Filtering thresholds file  
`-c`: Confidence estimates file  

Optionally, manually defined groups of closely related species can be included in the results:  
`-s`: Species groups file
 
See [here](https://github.com/BenKuhnhaeuser/GPID/wiki/4.-Sample-identification) for a full list of arguments and detailed instructions on the requirements for each argument.


## Pipeline outputs
The GeneParliamentID pipeline summarises all individual gene identifications in a Gene Parliament, which represents the percentage of genes supporting all competing identifications.  
  
The Gene Parliament is presented both as a table and as a figure.

Example Gene Parliament figure:  
<img width="800" alt="Entandrophragma_angolense_CQL_EA9_geneparliament" src="https://github.com/user-attachments/assets/70cae691-0e4b-45fd-bbcb-15b02621c704" />
  
For a detailed description of the Gene Parliament figure and table and their interpretation, see [Interpretation](https://github.com/BenKuhnhaeuser/GPID/wiki/5.-Interpretation).


## Tutorial  
We provide a worked example of the GeneParliamentID pipeline using test data and calibration files.  
This allows you to quickly familiarise yourself with:  
- input file requirements and formats
- how to run the pipeline
- the output files produced

It might also provide a useful template for setting up your own analyses.

To access the worked example, see [Tutorial](https://github.com/BenKuhnhaeuser/GPID/wiki/6.-Tutorial).
