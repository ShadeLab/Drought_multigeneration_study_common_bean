## Github Repository for
# Evaluating the legacy of drought exposure on root and rhizosphere bacterial microbiomes over two plant generations.  
## A Sulesky-Grieb†, AF Bintarti†, J Colovas, B Marolleau, T Boureau, M Simonin, M Barret and A Shade
<i>This work is submitted for peer review.</i>


### Data
The raw sequence data for this study are available in the NCBI SRA under bioproject [PRJNA1058980](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1058980/)


### To cite this work or code
Coming soon.


### Abstract
	Drought is a critical risk in developing countries for staple crops like common bean (Phaseolus vulgaris L.). We conducted an experiment to understand the legacy effects of repeated drought exposure across plant generations on the root and rhizosphere microbiome of the common bean, hypothesizing that a legacy of exposure improves overall plant microbiome resilience. We profiled the bacterial microbiome using marker gene amplicon sequencing over two plant generations in a complete factorial design for two common bean genotypes, Red Hawk and Flavert. We also repeated the experiment for Red Hawk in the two distinct production soils of Pays de la Loire, France, and Michigan, USA. Despite clear and relatively consistent drought effects on plant phenotype, there was neither response of the Red Hawk microbiomes to drought, nor a notable legacy of drought exposure. For Flavert, there was a legacy drought effect for the second generation in the rhizosphere microbiome beta diversity. These results suggest that a cross-generational legacy of drought on the belowground plant microbiome may not be consistent and can be affected by soil origin and host genotype. It also suggests a potential for resistance in the belowground plant microbiome to drought. 


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
This work was supported by the United States Department of Agriculture [grant number 2019-67019-29305](https://portal.nifa.usda.gov/web/crisprojectpages/1018838-inheritance-of-abiotic-stress-tolerance-through-seed-microbiome-modification.html) to AS and MB, and by the Michigan State University [Plant Resilience Institute](https://plantresilience.msu.edu/) to AS. This work was co-funded by the European Union [grant number ERC, MicroRescue, 101087042](https://cordis.europa.eu/project/id/101087042) to AS. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them. AS acknowledges support from the United States Department of Agriculture National Institute of Food and Agriculture and Michigan State University AgBioResearch, and the Centre National de la Recherche Scientifique (CNRS), France. The Angers Plant Phenotyping Platform PHENOTIC (DOI: 10.17180/ykbz-2v85) is acknowledged for the production and phenotyping of plants.

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)  
[Ecologie Microbienne Lyon](https://www.ecologiemicrobiennelyon.fr/eng)  
[IRHS](https://irhs.angers-nantes.hub.inrae.fr/)  
