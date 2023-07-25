
# Data overview

- directories
  - *annotation*: gene annotations, chromosome lengths, sample and metastasis pair annotation
  - *bin*: covTest R package and regNet R package (installation see below)
  - *conda*: the required conda environment will be installed there (installation see below)
  - *data*: expression and methylation values of the cohort
  - *FiguresTables*: figures and tables of the paper will be stored there (see the analysis scripts below)
  - *regNet*: where the internal analysis data of the scripts are stored, including the validation cohort
  - *scripts*: all jupyter notebooks of the analysis explained below (*.r.ipynb)

 
# Analysis pipeline

- To just look at the jupyter notebook scripts without installing anything, starting an environment or jupyter notebooks: all scripts are stored as html files in scripts/html/

- alternatively you can run the R scripts saved in scripts/R/*.r in e.g. Rstudio instead of using jupyter notebooks. But the conda environment must be set up as explained below

## Prerequisites

### download additional data

- Additional data for the analysis is stored at Zenodo (doi: 10.5281/zenodo.8182913) 
- this includes
  - annotation.tar.gz - gene annotation, pathway definitions, R data with sample and group annotation
  - covTest_1.02.tar.gz - covTest R package needed (see below)
  - regNet-master.zip - regNet R package (see below)
  - data.tar.gz - MelBrainSys and TCGA expression and methylation data used
  
### Setting up conda environments

```
conda create --prefix conda -c conda-forge r-base=4 r-irkernel notebook r-ggplot2 r-scales r-gridextra r-circlize r-venndiagram r-ggrepel r-pheatmap r-reshape2 r-glmnet r-glmpath r-lars r-stringr r-xlsx r-readxl r-devtools r-pvclust
```

- activate and add 2 bioconductor packages:

```
conda activate ./conda
conda install bioconductor-org.hs.eg.db
conda install bioconductor-clusterprofiler
```

### Using jupyter notebooks
- start jupyter notebook to view and/or execute the *.r.ipynb scripts:

```
jupyter notebook
```
- if you run the code on a powerful remote server, use ssh and port forwarding:
  - on the server (where the above conda environment is running):

```
jupyter notebook --port=9000 --no-browser
```

  - on your local machine where you edit the notebooks in an internet browser:
  
```
ssh -N -L 9000:localhost:9000 userlogin@server.com
```


### Installing covtest and regNet

#### Get regNet
- download to bin directory

```
cd bin
git clone https://github.com/seifemi/regNet.git
``` 
- or download and unpack 
- create an empty regNet folder in your conda/lib/R/library/

#### Get covtest
- download the latest version https://cran.r-project.org/src/contrib/Archive/covTest/  ->  covTest_1.02.tar.gz
- or use the file here: bin/covTest_1.02.tar.gz (after download from Zenodo, see above)
```
cd bin
tar -xvf covTest_1.02.tar.gz
```

#### Then install both with jupyter notebook

- scripts/setup-covTest-and-regNet.r.ipynb
 
## Data preparation

- download additional data from Zenodo (doi: 10.5281/zenodo.8182913), see above
- unpack the gz files in regNet/Data/ and data/ and the file annotation/pathway-definitions.xls.gz

- scripts/create-training-test-data.r.ipynb
  - create test and training datasets from TCGA data
  - creates many files in regNet/Data/ - TestSet_ExpressionData_regNet_Run_1.txt, TestSet_MethylationData_regNet_Run_1.txt, TrainSet_ExpressionData_regNet_Run_1.txt, TrainSet_MethylationData_regNet_Run_1.txt ...

- scripts/prepare-pathways.r.ipynb
  - create an R object from sheets of the pathway definition table S4 (metabolic pathways, ...) and save for the scripts used later pathway-categories.rds (already done)


## Define metastasis pairs, subgroups, sample annotation colors

- scripts/define-metPairs-subgroups-sampleAnnotation.r.ipynb
  - define metastasis pairs of sample IDs found in expression and methylation data
  - define annotation colors for plots later on
  - saved in annotation/samplePairs-annotation-colors-clusters.rds

## Network inference

- scripts/network-training.r.ipynb
- for each of the 1-n networks: 
  - this will create the basic regNet folder structure for each network
  - and train a network with the corresponding training data (expression and methylation files from above)
  - the results are split up into single parts for each network, which are then combined
- Download NetworkModels.tar from Zenodo (doi: 10.5281/zenodo.8182913) to avoid training again, as this takes lots of RAM, CPU and time. The archive contains for each network (TrainNetwork-1/, .. TrainNetwork-25/) a subfolder NetworkModel/WholeNetwork/ with files inside. These can be used with the following scripts.

## Trained network evaluation

- networks-connectivity.r
- creates 
  - connectivity table for paper - Suppl-Table-5-connectivity-stats.xls
  - Dot plot of connections - SupplFigure-3-source-target-interactions.png
  - known melanoma driver gene ranking in networks - Suppl-Table-6-key-driver-ranks.xls

## Evaluate prediction of networks

- scripts/network-prediction-evaluation.r.ipynb
  - predicts expression of the test data with the trained and random networks
  - predicts expression of MelBrainSys data
  - creates correlation distribution plots of these predicted and original expression: Figure-2-gene-expr-prediction-correl.png
  
## Neighboring gene auto-correlation

- calculate-gene-gene-expr-correl.r.ipynb
- calculate auto-correlation of the expression of neighboring genes based on TCGA expression values, Figure S1

## Network propagation
- scripts/network-propagation.r.ipynb
  - computes network flow/impact matrix of each gene on every other gene for each sample

## Altered genes per metastases pair
- scripts/altered-gene-per-metastases-pair.r.ipynb
  - defines and prepares altered genes, i.e. differential promoter methylation + expression trend
  - saved in altered-genes-per-patient.rds
  - overlap of these genes between metastasis pairs
  - produces subfigures for Suppl. Figure S2 
  
## Calculate log impacts
- scripts/calculate-log-impacts.r.ipynb
- as basis for other scripts
- calculates for each metastases pair 
  - average impacts log-ratios (brain/extracranial) on pathways
  - average impacts log-ratios on genes
- from the altered genes of all met. pairs
- calculates for each sample and met. pair the average impact ratios (brain/extracranial) of the altered genes on all pathways
- saved for later: metPairs-impactRatios-onPathways.rds, metPairs-impactRatios-onAllTargetGenes.rds


## plot log impact distribution and perform cluster analysis

- scripts/impact-logratios-cluster-analysis-TableS7-S9.r.ipynb
- violin plots of impact log2-ratios, Suppl. Figure S4
- cluster heatmaps of log-ratio impacts on pathways, Figure 3
- Suppl. Table S7: log-impact and log-expression by pathway
- Suppl. Table S9: expression impact quadrants shared pathways

## Impact expression correlation

- scripts/impact-expression-correl-Fig4-SuppFigS5.r.ipynb
- calculates correlations between average log expression ratio and average log impact ratio of the altered genes on all pathways of each metastases pair
- creates SupplFigure 5 with scatter plots
- and Figure 4 with barplots of the correlations 

## Overrepresentation analysis
- scripts/overrepresentation-analyses.r.ipynb
- overrepresentation analyses of the top ranked genes of each subgroup
- determine exclusive shared overrepresented genes that are consistent within each group but unique for the group
- saved in Suppl. Table S8

## Target candidate genes
- scripts/target-candidates-circosPlots-Fig5-SupplFigS6-SupplTableS10.r.ipynb
- produces
  - Figure 5: circos plots of source and target genes
  - SupplTable S10 target genes mean impact and expr. ratios per group
  - SupplFigure S6 target gene candidates scatterplot impact vs. expression ratios

## Validation

- must copy all TCGA trained networks from discovery to validation cohort subdirectories
  - in regNet/ directory of discovery cohort:
  
  ```
  tar -cvf trained-NWs.tar.gz TrainNetwork-*/NetworkModel/WholeNetwork/*
  mv trained-NWs.tar.gz validation-cohort/
  cd validation-cohort/
  tar xvf trained-NWs.tar.gz 
  ```
- copy the methylation and expression data of the validation cohort (MelBrainSys_ExpressionData_2022_allNeededGenes_regNet.txt.gz, MelBrainSys_MethylationData_2022_allNeededGenes_regNet.txt.gz) from data (from Zenodo) to the directory regNet/validation-cohort/Data
- scripts/validation-cohort-network-flow-propagation.r.ipynb
  - generates folder structure (remaining subfolders in TrainNetwork-* )
  - computes network flow/impact matrix of each gene on each other gene for each sample
- scripts/validation-SupplTableS12-SupplFigureS7.r.ipynb
  - computes impact log-ratios of all source genes on all genes
  - clusters impact ratios of discovery and validation cohort together (Suppl. Figure S7)
  - determines validation cohort target gene candidates (Suppl. Table S10)
  - rankings of discovery cohort target genes in valid. cohort and impact log-ratios (Suppl. Table 12)


## TCGA classification and association

- correlate TCGA sample expression with MelBrainSys samples, assign MelBrainSys subgroups 
- calculate associations of assigned subgroups with clinical / molecular features, make Kaplan-Meier survival analysis
- download Suppl. Table S3 of "Genomic Classification of Cutaneous Melanoma", Cell. 2015 Jun 18; 161(7): 1681–1696, doi: 10.1016/j.cell.2015.05.044
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4580370/  -> NIHMS698912-supplement-3.xlsx
  - save it to annotation folder
- scripts/TCGA-associations.r.ipynb
- produces
  - Figure-6-TCGA-associations.png
  - SupplFigure-S8-TCGA-survival.png

