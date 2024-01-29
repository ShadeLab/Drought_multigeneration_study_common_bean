# Drought_multigeneration_study_common_bean
 Sequence processing and analysis for multigeneration drought experiment


## Github Repository for
# Belowground plant microbiome resistance and response after two generations of drought exposure in common bean
## Abby Sulesky-Grieb, Marie Simonin, A. Fina Bintarti, Brice Marolleau, Matthieu Barret, and Ashley Shade
<i>This work is not published but will soon be available on bioRxiv.</i>


### Data
The raw sequence data for this study are available in the NCBI SRA under bioproject [PRJNA1058980](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1058980/)



### To cite this work or code
Coming soon.


### Abstract
Coming soon.


### Contents

Code is split up into six directories. Files necessary to run R code are located in [R_analysis_files](https://github.com/ShadeLab/Drought_multigeneration_study_common_bean/tree/main/R_analysis_files).

#### Sequence processing
The root and rhizosphere samples from Michigan were process separately. Code used for sequence processing including read QC, ASV clustering, taxonomy assignment, can be found under  [Multigen_Root_processing.Rmd](https://github.com/ShadeLab/Drought_multigeneration_study_common_bean/blob/main/Multigen_Root_processing.Rmd) and [Multigen_Rhizosphere_processing.Rmd](https://github.com/ShadeLab/Drought_multigeneration_study_common_bean/blob/main/Multigen_Rhizosphere_processing.Rmd). Scripts were run in QIIME2 using SLURM on the MSU HPCC using slurm batch files with suffix .sb. 

#### Sequence decontamination and rarefaction
Code for sequence decontamination by extraction group was also completed separately for Michigan root and rhizosphere samples: [root_OTU_decontam.Rmd](https://github.com/ShadeLab/Drought_multigeneration_study_common_bean/blob/main/root_OTU_decontam.Rmd) and [rhizo_OTU_decontam.Rmd](https://github.com/ShadeLab/Drought_multigeneration_study_common_bean/blob/main/rhizo_OTU_decontam.Rmd). Output ASV table, taxonomy and metadata files from QIIME2 were used to create a Phyloseq object in R. 

After decontamination, the datasets were combine and rarefied together. Decontaminated and rarefied phyloseq object is labelled root_rhizo_control_drought_rarefiedphyloseq copy.rds in the R_analysis_files folder.  

#### Microbiome Analysis
Formal analysis of the Michigan and France datasets can be found under [combined_16s_analysis.Rmd](https://github.com/ShadeLab/Drought_multigeneration_study_common_bean/blob/main/combined_16s_analysis.Rmd). All analysis was run with R and code was run in Rmarkdown. Input files of R phyloseq oobjects for both datasets can be found in the R_analysis_files folder.

#### Plant Phenotypic trait analysis
Plant phenotypic measurements are located in the R_analysis_files folder. Analysis code can be found at [PlantHealth_analysis.Rmd](https://github.com/ShadeLab/Drought_multigeneration_study_common_bean/blob/main/PlantHealth_analysis.Rmd).

### Funding
This work was supported by the United States Department of Agriculture award 2019-67019-29305, and by the Michigan State University [Plant Resilience Institute](https://plantresilience.msu.edu). 

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)