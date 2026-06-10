
#title: "combined_16S_new_analysis"
#author: "Abby Sulesky & Ari Fina Bintarti"
#date: "2025-07-30"


install.packages("remotes")
install.packages("vegan")
library(remotes)
#remotes::install_github("david-barnett/microViz", force = TRUE)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(microViz)
library(corncob)
library(ggraph)
library(DT)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(decontam)
library(scales)
library(vegan)
library(microbiome)
library(ComplexHeatmap)
library(dplyr)
library(patchwork)
library(devtools)
library(pairwiseAdonis)
library(indicspecies)
library(VennDiagram)



##### Alpha Diversity Analysis #####


############## 1. MSU DATA ##############


# Un-rarefied data
MSU_phyloseq_unrare <- readRDS("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/root_rhizo_combined_phyloseq.rds")
MSU_phyloseq_unrare
MSU_meta_unrare <- as.data.frame(sample_data(MSU_phyloseq_unrare))
# subset only control and drought in Treatment
MSU_phyloseq_unrare_cont_dro <- MSU_phyloseq_unrare %>% 
  ps_filter(Treatment != "Nutrient", .keep_all_taxa = TRUE)
MSU_phyloseq_unrare_cont_dro
# subset only control and drought in G1_treatment
MSU_phyloseq_unrare_cont_dro2 <- MSU_phyloseq_unrare_cont_dro %>% 
  ps_filter(G1_treatment != "Nutrient", .keep_all_taxa = TRUE)
MSU_phyloseq_unrare_cont_dro2 #133 samples retained
# check sequencing depth
sort(colSums(otu_table(MSU_phyloseq_unrare_cont_dro2), na.rm = FALSE, dims = 1), decreasing = F) # not rarefied, lowest= 20976 reads
# MSU metadata
MSU_meta_ed <- as.data.frame(sample_data(MSU_phyloseq_unrare_cont_dro2))
MSU_meta_ed
#write_csv(MSU_meta_ed, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/MSU_meta_ed.csv")
#setwd("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/")
MSU_meta_ed <- read.csv("MSU_meta_ed.csv", check.names = F,row.names = 1)

# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc", force=T)
library(metagMisc)


#### Perform Multiple Rarefactions and Calculate Average Alpha Diversity (Richness) across Rarefactions ###

set.seed(13)
MSU_multrare_alpha <- phyloseq_mult_raref_div(MSU_phyloseq_unrare_cont_dro2,
                                              SampSize = min(sample_sums(MSU_phyloseq_unrare_cont_dro2)),
                                              iter = 100,
                                              divindex = "Observed",
                                              parallel = FALSE,verbose = TRUE)
MSU_multrare_alpha
# save the output in the local computer
# write.csv(MSU_multrare_alpha, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/MSU_multi_rare_alpha.csv",row.names=TRUE)

# read back the edited multi-rarefactions alpha diversity
MSU_multrare_alpha <- read.csv("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/MSU_multi_rare_alpha.csv",
                               row.names = 1)

# check if two dataframes have the same row names before adding alpha diversity data in the meta data
identical(rownames(MSU_meta_ed), rownames(MSU_multrare_alpha))
# adding observed/richness (average from multiple rarefactions) to the metadata
MSU_meta_ed$Observed <- MSU_multrare_alpha[rownames(MSU_meta_ed), "Estimate"]

# code for background theme for ggplot so I don't have to have it in every plot
my_theme <- theme(panel.background = element_rect(fill = "white", colour = "white"), 
                  panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', colour = "light gray"),
                  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', colour = "light gray"))



############## MSU: Generation 1 ##############


# subset G1 root data
MSU_G1_Root <- MSU_meta_ed %>%
  filter(Generation == "G1", Compartment == "root")
str(MSU_G1_Root)

# plot richness - G1 root
MSU_G1_Root_alpha <- ggplot(MSU_G1_Root, aes(x=Treatment, y=Observed, col=Treatment)) + 
  geom_boxplot(aes(color=Treatment)) + 
  geom_point(size=3, shape=21, fill="white") +
  labs(title="(USA)") +
  scale_y_continuous(limits = c(300, 1000))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), 
                     labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  my_theme +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), 
        plot.tag = element_text(face="bold")) + 
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Flavert"="Flavert", 
                                                                                 "Red_Hawk"="Red Hawk"))) 
MSU_G1_Root_alpha

# statistical test G1 root using Welch t-test (without assuming the homogeneity of the data) 

# check normality
library(ggpubr)
ggqqplot(MSU_G1_Root$Observed) # Okay
shapiro.test(MSU_G1_Root$Observed) # NS
# stats - G1 Root
MSU_G1_Root_stats <- t.test(Observed ~ G1_treatment, data=MSU_G1_Root)
MSU_G1_Root_stats # Treatment Not Significant

# subset G1 Rhizosphere data
MSU_G1_Rhizo <- MSU_meta_ed %>%
  filter(Generation == "G1", Compartment == "rhizosphere")
str(MSU_G1_Rhizo)

# plot richness - G1 Rhizosphere
MSU_G1_Rhizo_alpha <- ggplot(MSU_G1_Rhizo, aes(x=Treatment, y=Observed, col=Treatment)) + 
  geom_boxplot(aes(color=Treatment)) + 
  geom_point(size=3, shape=21, fill="white") +
  labs(title="(USA)") +
  scale_y_continuous(limits = c(800, 1800))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), 
                     labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  my_theme +
  theme(legend.position="none", 
        axis.title.y=element_blank(), 
        plot.tag = element_text(face="bold"),
        axis.title.x = element_blank()) + 
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Flavert"="Flavert", 
                                                                                 "Red_Hawk"="Red Hawk"))) 
MSU_G1_Rhizo_alpha

# statistical test G1 Rhizo using Welch t-test (without assuming the homogeneity of the data) 
# check normality
ggqqplot(MSU_G1_Rhizo$Observed) # a bit off
shapiro.test(MSU_G1_Rhizo$Observed) # NS
# stats - G1 Rhizo
MSU_G1_Rhizo_stats <- t.test(Observed ~ G1_treatment, data=MSU_G1_Rhizo)
MSU_G1_Rhizo_stats # Treatment Not Significant


############## MSU: Generation 2 #############


# subset G2 root data
MSU_G2_Root <- MSU_meta_ed %>%
  filter(Generation == "G2", Compartment == "root")
str(MSU_G2_Root)

# plot richness - G2 root
MSU_G2_Root_alpha <- ggplot(MSU_G2_Root, aes(x=G1_treatment, y=Observed, group=G1_G2, color = G1_G2)) + 
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  geom_point(position=position_dodge(width=0.8),aes(group=G1_G2, col=G1_G2),size=3, shape=21, fill="white")+
  labs(title="(USA)") + 
  scale_y_continuous(limits = c(300, 1650))+
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment", 
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), 
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) + ylab("ASVs per Root Sample") + 
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(), 
        plot.tag = element_text(face="bold")) +
  my_theme + facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Flavert"="Flavert", 
                                                                                            "Red_Hawk"="Red Hawk"))) 
MSU_G2_Root_alpha

# statistical test -G2 root with two-way ANOVA

# check normality
plot(MSU_G2_Root_stats, 2) #look okay
shapiro.test(residuals(MSU_G2_Root_stats)) #NS, good
# homogeneity
plot(MSU_G2_Root_stats,1)
MSU_G2_Root %>% 
  levene_test(Observed_bc ~ G1_treatment*G2_treatment) # its okay with transformed data
# predict the best transformation
set.seed(13)
test <- bestNormalize(MSU_G2_Root$Observed)
test # use box cox
# data transformation
set.seed(13)
bc_MSU_G2_Root <- boxcox(MSU_G2_Root$Observed)
MSU_G2_Root$Observed_bc <- predict(bc_MSU_G2_Root)
# apply ANOVA
MSU_G2_Root_stats <- aov(Observed_bc ~ G1_treatment*G2_treatment, data=MSU_G2_Root)
summary(MSU_G2_Root_stats) # NS
#                          Df  Sum Sq Mean Sq F value Pr(>F)  
#G1_treatment               1   0.15  0.1482   0.151  0.699
#G2_treatment               1   2.37  2.3707   2.422  0.127
#G1_treatment:G2_treatment  1   1.42  1.4210   1.452  0.235
#Residuals                 44  43.06  0.9786 
# post-hoc -G2 Root
MSU_G2_Root_tuk <- tukey_hsd(MSU_G2_Root, Observed ~ G1_treatment*G2_treatment)
MSU_G2_Root_tuk #NS


# subset the G2 Rhizosphere data
MSU_G2_Rhizo <- MSU_meta_ed %>%
  filter(Generation == "G2", Compartment == "rhizosphere")
str(MSU_G2_Rhizo)

# plot richness - G2 Rhizosphere
MSU_G2_Rhizo_alpha <- ggplot(MSU_G2_Rhizo, aes(x=G1_treatment, y=Observed, group=G1_G2, color = G1_G2)) + 
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  geom_point(position=position_dodge(width=0.8),aes(group=G1_G2, col=G1_G2),size=3, shape=21, fill="white")+
  labs(title="(USA)", x="G1 Treatment") + 
  scale_y_continuous(limits = c(500, 2600))+
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment", 
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), 
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) +
  ylab("ASVs per Rhizosphere Sample") + 
  theme(legend.position="none", axis.title.y=element_blank(), 
        plot.tag = element_text(face="bold"), axis.title.x = element_text(margin = margin(t = 15)))  +
  my_theme + facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Flavert"="Flavert", 
                                                                                            "Red_Hawk"="Red Hawk"))) 
MSU_G2_Rhizo_alpha

# statistical test -G2 Rhizosphere with two-way ANOVA

# check normality
plot(MSU_G2_Rhizo_stats, 2) #look okay
shapiro.test(residuals(MSU_G2_Rhizo_stats)) #NS, good
# homogeneity
plot(MSU_G2_Rhizo_stats,1)
MSU_G2_Rhizo %>% 
  levene_test(Observed ~ G1_treatment*G2_treatment) # NS, good
# apply ANOVA
MSU_G2_Rhizo_stats <- aov(Observed ~ G1_treatment*G2_treatment, data=MSU_G2_Rhizo)
summary(MSU_G2_Rhizo_stats) # significant interaction G1 x G2
#                           Df  Sum Sq Mean Sq F value Pr(>F)  
#G1_treatment               1   84582   84582   0.685 0.4124  
#G2_treatment               1   90090   90090   0.729 0.3977  
#G1_treatment:G2_treatment  1  512061  512061   4.145 0.0477 *
#Residuals                 45 5559670  123548                 

#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# post-hoc -G2 Rhizosphere
MSU_G2_Rhizo_tuk <- tukey_hsd(MSU_G2_Rhizo, Observed ~ G1_treatment*G2_treatment)
MSU_G2_Rhizo_tuk #NS

MSU_G2_Rhizo_pwc <- MSU_G2_Rhizo %>%
  dplyr::group_by(G1_treatment) %>%
  tukey_hsd(Observed ~ G2_treatment)
MSU_G2_Rhizo_pwc 
#G1_treatment term           group1  group2  null.value estimate conf.low conf.high  p.adj p.adj.signif
#1 Control      G2_treatment Control Drought          0    -115.  -404.        175. 0.421  ns          
#2 Drought      G2_treatment Control Drought          0     294.    -4.67      593. 0.0534 ns 



#############  2. INRAE DATA ############# 


# # Un-rarefied data
INRAE_root_rhizo_CD <- readRDS("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/INRAE_root_rhizo_CD.rds")
INRAE_root_rhizo_CD
# check sequencing depth
sort(colSums(otu_table(t(INRAE_root_rhizo_CD)), na.rm = FALSE, dims = 1), decreasing = T) # not rarefied 


#### Perform Multiple Rarefactions and Calculate Average Alpha Diversity (Richness) across Rarefactions ###

set.seed(13)
INR_multrare_alpha <- phyloseq_mult_raref_div(INRAE_root_rhizo_CD,
                                              SampSize = 2500,
                                              iter = 100,
                                              divindex = "Observed",
                                              parallel = FALSE,verbose = TRUE)
INR_multrare_alpha
# save the output in the local computer
# write.csv(INR_multrare_alpha, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/INR_multi_rare_alpha.csv",row.names=TRUE)

# INRAE metadata
INRAE_meta_ed <- as.data.frame(sample_data(INRAE_root_rhizo_CD))
INRAE_meta_ed
#write_csv(INRAE_meta_ed, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/INRAE_meta_ed.csv")
setwd("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/")
INRAE_meta_ed <- read.csv("INRAE_meta_ed.csv", check.names = F,row.names = 1) # still has 120 samples

# read back the edited multi-rarefactions alpha diversity
INR_multrare_alpha <- read.csv("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/INR_multi_rare_alpha.csv",
                               row.names = 1) # it has 110 samples because 10 samples are lost during rarefaction

# adding observed/richness (average from multiple rarefactions) to the metadata
# 1. Subset 'INRAE_meta_ed' to only rows that also exist in 'INR_multrare_alpha'
INRAE_meta_ed_filt <- INRAE_meta_ed[rownames(INRAE_meta_ed) %in% rownames(INR_multrare_alpha), ]

# 2. Add the new column from df1 to df2_filtered
INRAE_meta_ed_filt$Observed <- INR_multrare_alpha[rownames(INRAE_meta_ed_filt), "Estimate"]
# check which samples are dropped
setdiff(rownames(INRAE_meta_ed), rownames(INR_multrare_alpha))



############## INRAE: Generation 1 ##############


# Subset INRAE G1 Root and Rhizosphere
INRAE_G1 <- INRAE_meta_ed_filt %>%
  filter(Generation == "G1")

INRAE_G1_Root <- INRAE_G1 %>%
  filter(Compartment == "Root")

INRAE_G1_Rhizo <- INRAE_G1 %>%
  filter(Compartment == "Rhizosphere")

# plot richness - G1 Root
INRAE_G1_Root_alpha <- ggplot(INRAE_G1_Root, aes(x=Treatment, y=Observed, col=Treatment)) + 
  geom_boxplot(aes(color=Treatment)) + 
  geom_point(size=3, shape=21, fill="white") +
  labs(title="(USA)") +
  scale_y_continuous(limits = c(40, 200))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), 
                     labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  labs(title="G1 Root Alpha Diversity (FR)") + 
  ylab("ASVs per Root Sample") +
  my_theme +
  theme(legend.position="none", axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        plot.tag = element_text(face="bold")) +
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Flavert"="Flavert", 
                                                                                 "Red_Hawk"="Red Hawk")))
  #facet_wrap(vars(Genotype),nrow=1, scales="free")
INRAE_G1_Root_alpha

# statistical test G1 Root using two-way ANOVA 

# check normality
library(rstatix)
plot(INR_G1_Root_stats, 2)
shapiro.test(residuals(INR_G1_Root_stats)) # 
# checking outliers
INR_G1_Root_out <- INRAE_G1_Root %>%
  group_by(Genotype, G1_treatment) %>%
  identify_outliers(Observed) # just 1 extreme
INR_G1_Root_out
# homogeneity
plot(INR_G1_Root_stats,1)
INRAE_G1_Root %>% levene_test(Observed ~ Genotype*G1_treatment) # NS, so it is Good
#INRAE_G1_Root <- INRAE_G1_Root[, -120]

library(MASS)
MASS::boxcox(aov(Observed ~ Genotype*G1_treatment, data=INRAE_G1_Root)) # lambda 1 is within the 95% of CI, meaning that transformation is not needed

# stats - G1 Root
INR_G1_Root_stats <- aov(Observed ~ Genotype*G1_treatment, data=INRAE_G1_Root) 
summary(INR_G1_Root_stats) # Genotype is Significant & Treatment is Not Significant
#                       Df Sum Sq Mean Sq F value Pr(>F)  
#Genotype               1   7133    7133   7.651 0.0138 *
#G1_treatment           1   3079    3079   3.303 0.0879 .
#Genotype:G1_treatment  1   3295    3295   3.534 0.0785 .
#Residuals             16  14917     932           
              
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# post-hoc -G1 Root
TukeyHSD(INR_G1_Root_stats)
INR_G1_Root_tuk <- tukey_hsd(INRAE_G1_Root, Observed ~ Genotype*G1_treatment)
INR_G1_Root_tuk

#term                  group1           group2           null.value estimate conf.low conf.high  p.adj p.adj.signif
# Genotype              Flavert          Red_Hawk                  0  -37.8     -66.7      -8.82 0.0138 *    
# G1_treatment          Control          Drought                   0   24.8      -4.13     53.8  0.0879 ns
# Genotype:G1_treatment Flavert:Control  Flavert:Drought           0   50.5      -4.76    106.   0.0796 ns
# Genotype:G1_treatment Red_Hawk:Control Red_Hawk:Drought          0   -0.854   -56.1      54.4  1      ns
# Genotype:G1_treatment Flavert:Drought  Red_Hawk:Drought          0  -63.4    -119.       -8.19 0.0218 *

# checking pairwise comparison control vs drought within each genotype
INR_G1_Root_pwc <- INRAE_G1_Root %>%
  dplyr::group_by(Genotype) %>%
  tukey_hsd(Observed ~ G1_treatment)
INR_G1_Root_pwc

# check if include compartment
INR_G1_aov <- aov(Observed ~ Genotype*G1_treatment*Compartment, data=INRAE_G1) 
summary(INR_G1_aov) # Genotype is Significant & Treatment is Not Significant
# post-hoc -G1 Root
INR_G1_tuk <- tukey_hsd(INRAE_G1, Observed ~ Genotype*G1_treatment*Compartment)
INR_G1_tuk #NS


# plot richness - G1 Rhizosphere
INRAE_G1_Rhizo_alpha <- ggplot(INRAE_G1_Rhizo, aes(x=Treatment, y=Observed, col=Treatment)) + 
  geom_boxplot(aes(color=Treatment)) + 
  geom_point(size=3, shape=21, fill="white") +
  labs(title="(USA)") +
  scale_y_continuous(limits = c(100, 450))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), 
                     labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  labs(title="G1 Rhizosphere Alpha Diversity (FR)") +  
  ylab("ASVs per Rhizosphere Sample") +
  my_theme +
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Flavert"="Flavert", 
                                                                                 "Red_Hawk"="Red Hawk"))) +
  theme(legend.position="none", 
        plot.tag = element_text(face="bold"),
        axis.title.x = element_blank())
INRAE_G1_Rhizo_alpha

# statistical test G1 Rhizosphere using two-way ANOVA 

# check normality
plot(INR_G1_Rhizo_stats, 2)
shapiro.test(residuals(INR_G1_Rhizo_stats)) # Good
# homogeneity
plot(INR_G1_Rhizo_stats,1)
INRAE_G1_Rhizo %>% levene_test(Observed_bc ~ Genotype*G1_treatment) # Good with transformed data

INR_G1_Rhizo_bc <- MASS::boxcox(aov(Observed ~ Genotype*G1_treatment, data=INRAE_G1_Rhizo)) #The peak of the log-likelihood is now around λ ≈ 1
#The 95% confidence interval includes λ = 1
INR_G1_Rhizo_lambda <- INR_G1_Rhizo_bc$x[which.max(INR_G1_Rhizo_bc$y)]
INR_G1_Rhizo_lambda
# box-cox transformation
if (INR_G1_Rhizo_lambda == 0) {
  INRAE_G1_Rhizo$Observed_bc <- log(INRAE_G1_Rhizo$Observed)
} else {
  INRAE_G1_Rhizo$Observed_bc <- (INRAE_G1_Rhizo$Observed^INR_G1_Rhizo_lambda - 1) / INR_G1_Rhizo_lambda
}

# stats - G1 Rhizosphere
INR_G1_Rhizo_stats <- aov(Observed_bc ~ Genotype*G1_treatment, data=INRAE_G1_Rhizo) 
summary(INR_G1_Rhizo_stats) # Genotype is Not Significant & Treatment is Not Significant
#                       Df Sum Sq Mean Sq F value Pr(>F)  
#Genotype               1 5.604e-06 5.604e-06   3.395 0.0840 .
#G1_treatment           1 4.470e-07 4.470e-07   0.271 0.6100  
#Genotype:G1_treatment  1 5.393e-06 5.393e-06   3.267 0.0895 .
#Residuals             16 2.641e-05 1.651e-06                  

# post-hoc -G1 Rhizosphere
INR_G1_Rhizo_tuk <- tukey_hsd(INRAE_G1_Rhizo, Observed_bc ~ Genotype*G1_treatment)
INR_G1_Rhizo_tuk # Nothing is significant
# checking pairwise comparison control vs drought within each genotype
INR_G1_Rhizo_pwc <- INRAE_G1_Rhizo %>%
  dplyr::group_by(Genotype) %>%
  tukey_hsd(Observed_bc ~ G1_treatment)
INR_G1_Rhizo_pwc


############## INRAE: Generation 2 ##############


# Subset INRAE G1 Root and Rhizosphere
INRAE_G2 <- INRAE_meta_ed_filt %>%
  filter(Generation == "G2")

INRAE_G2_Root <- INRAE_G2 %>%
  filter(Compartment == "Root")

INRAE_G2_Rhizo <- INRAE_G2 %>%
  filter(Compartment == "Rhizosphere")

# Plot richness - INRAE G2 root
INRAE_G2_Root_alpha <- ggplot(data=INRAE_G2_Root, aes(x=G1_treatment, y=Observed, group=G1_G2, color = G1_G2)) + 
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  geom_point(position=position_dodge(width=0.8),aes(group=G1_G2, col=G1_G2),size=3, shape=21, fill="white")+
  labs(title="G2 Root Alpha Diversity (FR)") + 
  scale_y_continuous(limits = c(0, 45))+
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment", 
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), 
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) + ylab("ASVs per Root Sample") + 
  theme(legend.position="none", axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        plot.tag = element_text(face="bold")) +
  my_theme + facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Flavert"="Flavert", 
                                                                                            "Red_Hawk"="Red Hawk")))
INRAE_G2_Root_alpha

# statistical analysis for G2 - Root using three-way ANOVA

# check normality
plot(INR_G2_Root_stats, 2)
shapiro.test(residuals(INR_G2_Root_stats)) # Good with transformed data
# homogeneity
plot(INR_G2_Root_stats,1)
INRAE_G2_Root %>% levene_test(Observed_asin ~ Genotype*G1_treatment*G2_treatment) # NS, Good
# check the distribution
MASS::boxcox(aov(Observed ~ Genotype*G1_treatment*G2_treatment, data=INRAE_G2_Root)) # the optimum is between -1 and 0
# predict the best transformation
set.seed(13)
bN.INRAE_G2_Root <- bestNormalize(INRAE_G2_Root$Observed)
bN.INRAE_G2_Root # use Standardized asinh(x) Transformation
# transform the data using asinh(x) Transformation
set.seed(13)
asin.INRAE_G2_Root <- arcsinh_x(INRAE_G2_Root$Observed)
INRAE_G2_Root$Observed_asin <- asin.INRAE_G2_Root$x.t
# apply three-way ANOVA
INR_G2_Root_stats <- aov(Observed_asin ~ Genotype*G1_treatment*G2_treatment, data=INRAE_G2_Root) 
summary(INR_G2_Root_stats) # G2 is significant
#                                   Df Sum Sq Mean Sq F value Pr(>F)  
#Genotype                            1  2.098   2.098   2.490 0.1267  
#G1_treatment                        1  0.192   0.192   0.228 0.6372  
#G2_treatment                        1  5.917   5.917   7.023 0.0135 *
#Genotype:G1_treatment               1  0.547   0.547   0.650 0.4275  
#Genotype:G2_treatment               1  1.311   1.311   1.556 0.2234  
#G1_treatment:G2_treatment           1  1.001   1.001   1.188 0.2857  
#Genotype:G1_treatment:G2_treatment  1  0.028   0.028   0.034 0.8557  
#Residuals                          26 21.906   0.843                

# post-hoc -G2  Root
INR_G2_Root_tuk <- tukey_hsd(INRAE_G2_Root, Observed_asin ~ Genotype*G1_treatment*G2_treatment)
View(INR_G2_Root_tuk) # NS

# Separate between Genotype
IFlav_G2_Root <- INRAE_G2_Root %>%
filter(Genotype == "Flavert")

IRdHk_G2_Root <- INRAE_G2_Root %>%
  filter(Genotype == "Red_Hawk")

# check normality
plot(IFlav_G2_Root_aov, 2)
shapiro.test(residuals(IFlav_G2_Root_aov)) # NG
plot(IRdHk_G2_Root_aov, 2)
shapiro.test(residuals(IRdHk_G2_Root_aov)) # Good
# homogeneity
IFlav_G2_Root %>% levene_test(Observed ~ G1_treatment*G2_treatment) # Good
IRdHk_G2_Root %>% levene_test(Observed ~ G1_treatment*G2_treatment) # Good
# identify outliers
hist(IFlav_G2_Root$Observed)
out <- IFlav_G2_Root %>%
  group_by(G2_treatment) %>%
  identify_outliers (Observed) # there are outliers
# apply tWO-way ANOVA INRAE G2-Flavert
IFlav_G2_Root_aov <- aov(Observed ~ G1_treatment*G2_treatment, data=IFlav_G2_Root) 
summary(IFlav_G2_Root_aov) # NS
# apply tWO-way ANOVA INRAE G2-Red Hawk
IRdHk_G2_Root_aov <- aov(Observed ~ G1_treatment*G2_treatment, data=IRdHk_G2_Root) 
summary(IRdHk_G2_Root_aov) # G2_treatment               1  297.1  297.06   8.734  0.012 *
# post-hoc INRAE Red Hawk-G2  Root 
IRdHk_G2_Root_tuk <- tukey_hsd(IRdHk_G2_Root, Observed ~ G1_treatment*G2_treatment)
IRdHk_G2_Root_tuk #

IRdHk_G2_Root_pwc <- IRdHk_G2_Root %>%
  dplyr::group_by(G1_treatment) %>%
  tukey_hsd(Observed ~ G2_treatment)
IRdHk_G2_Root_pwc


# Plot Richness INRAE G2-Rhizosphere

#iG2rhizo_anno <- data.frame(x1 = 0.8, x2 = 1.75, 
                            #y1 = 315, y2 = 320, 
                            #xstar = 1.25, ystar = 322,
                            #lab = "*",
                            #Genotype = "Flavert")
iG2rhizo_anno <- data.frame(x1 = c(0.8, 0.8), x2 = c(1.75, 2.25), 
                            y1 = c(315, 330), y2 = c(320, 335), 
                            xstar = c(1.25, 1.5), ystar = c(322, 340),
                            lab = c("**", "**"),
                            Genotype = c("Flavert", "Flavert"))
iG2rhizo_anno


INRAE_G2_Rhizo_alpha <- ggplot(data=INRAE_G2_Rhizo, aes(x=G1_treatment, y=Observed, group=G1_G2, color = G1_G2)) + 
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  geom_point(position=position_dodge(width=0.8),aes(group=G1_G2, col=G1_G2),size=3, shape=21, fill="white")+
  ylab("ASVs per Rhizosphere Sample") +
  labs(title="G2 Rhizosphere Alpha Diversity (FR)", x="G1 Treatment") + 
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment", 
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), 
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) + 
  theme(legend.position="none", plot.tag = element_text(face="bold"),
        axis.title.x = element_text(margin = margin(t = 15))) +
  my_theme + facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Flavert"="Flavert", 
                                                                                            "Red_Hawk"="Red Hawk")))+
  scale_y_continuous(limits = c(0, 350))+
  geom_text(inherit.aes=FALSE, data = iG2rhizo_anno, aes(x = xstar,  y = ystar, label = lab), size=6) +
  geom_segment(inherit.aes=FALSE, data = iG2rhizo_anno, aes(x = x1, xend = x1, y = y1, yend = y2),colour = "black") +
  geom_segment(inherit.aes=FALSE, data = iG2rhizo_anno, aes(x = x2, xend = x2, y = y1, yend = y2), colour = "black") +
  geom_segment(inherit.aes=FALSE, data = iG2rhizo_anno, aes(x = x1, xend = x2, y = y2, yend = y2),colour = "black")

INRAE_G2_Rhizo_alpha

# statistical analysis for INRAE G2 - Rhizosphere using three-way ANOVA

# check normality
plot(INR_G2_Rhizo_stats, 2)
shapiro.test(residuals(INR_G2_Rhizo_stats)) # Good using transformed data
# homogeneity
plot(INR_G2_Rhizo_stats,1)
INRAE_G2_Rhizo %>% levene_test(Observed_ordN ~ Genotype*G1_treatment*G2_treatment) # 

# predict the best transformation
MASS::boxcox(aov(Observed ~ Genotype*G1_treatment*G2_treatment, data=INRAE_G2_Rhizo)) 
#INRAE_G2_Rhizo <- INRAE_G2_Rhizo[,c(-121,-122)]
set.seed(13)
bN.INRAE_G2_Rhizo <- bestNormalize(INRAE_G2_Rhizo$Observed)
bN.INRAE_G2_Rhizo #
# transform the data using orderNormalize Transformation
set.seed(13)
ordNorm.INRAE_G2_Rhizo <- orderNorm(INRAE_G2_Rhizo$Observed)
INRAE_G2_Rhizo$Observed_ordN <- ordNorm.INRAE_G2_Rhizo$x.t
# apply three-way ANOVA
INR_G2_Rhizo_stats <- aov(Observed_ordN ~ Genotype*G1_treatment*G2_treatment, data=INRAE_G2_Rhizo) 
summary(INR_G2_Rhizo_stats) 
#                                  Df Sum Sq Mean Sq F value  Pr(>F)   
#Genotype                            1  3.613   3.613   4.983 0.03378 * 
#G1_treatment                        1  6.621   6.621   9.132 0.00532 **
#G2_treatment                        1  0.193   0.193   0.266 0.61032   
#Genotype:G1_treatment               1  0.969   0.969   1.336 0.25748   
#Genotype:G2_treatment               1  0.416   0.416   0.574 0.45490   
#G1_treatment:G2_treatment           1  0.824   0.824   1.136 0.29555   
#Genotype:G1_treatment:G2_treatment  1  1.817   1.817   2.506 0.12463   
#Residuals                          28 20.301   0.725

# post-hoc INRAE-G2  Rhizosphere
INR_G2_Rhizo_tuk <- tukey_hsd(INRAE_G2_Rhizo, Observed ~ Genotype*G1_treatment*G2_treatment)
View(INR_G2_Rhizo_tuk)
# Fl vs. RH  ***
# G1_treatment  Control vs Drought   ****
# Fl G1 C vs.  Fl G1 D *
# Fl C_C vs. Fl D_C ****
# Fl C_C vs. Fl D_D ****
INR_G2_Rhizo_pwc <- INRAE_G2_Rhizo %>%
  group_by(Genotype) %>%
  pairwise_t_test(Observed ~ G1_G2, p.adjust.method = "BH")
INR_G2_Rhizo_pwc

# Separate between Genotype
IFlav_G2_Rhizo <- INRAE_G2_Rhizo %>%
  filter(Genotype == "Flavert")

IRdHk_G2_Rhizo <- INRAE_G2_Rhizo %>%
  filter(Genotype == "Red_Hawk")
# check normality
shapiro.test(residuals(IFlav_G2_Rhizo_aov)) # Good 
shapiro.test(residuals(IRdHk_G2_Rhizo_aov)) # Good
# homogeneity
IFlav_G2_Rhizo %>% levene_test(Observed ~ G1_treatment*G2_treatment) # robust, we have equal sample size
IRdHk_G2_Rhizo %>% levene_test(Observed ~ G1_treatment*G2_treatment) # Good
# order normalized transformation for INRAE G2 Flavert
set.seed(13)
#ordN.IFlav_G2_Rhizo <- orderNorm(IFlav_G2_Rhizo$Observed)
#IFlav_G2_Rhizo$Observed_ordN <- ordN.IFlav_G2_Rhizo$x.t

# apply tWO-way ANOVA-Flavert
IFlav_G2_Rhizo_aov <- aov(Observed ~ G1_treatment*G2_treatment, data=IFlav_G2_Rhizo) 
summary(IFlav_G2_Rhizo_aov) #G1_treatment               1  80096   80096  25.631 0.000173 ***
# apply tWO-way ANOVA-Red Hawk
IRdHk_G2_Rhizo_aov <- aov(Observed ~ G1_treatment*G2_treatment, data=IRdHk_G2_Rhizo) 
summary(IRdHk_G2_Rhizo_aov) #NS
# post-hoc INRAE Flavert-G2  Rhizosphere
IFlav_G2_Rhizo_tuk <- tukey_hsd(IFlav_G2_Rhizo, Observed ~ G1_treatment*G2_treatment)
IFlav_G2_Rhizo_tuk # C_C vs D_C **, C_C vs D_D **
# post-hoc INRAE Red Hawk-G2  Rhizosphere
IRH_G2_Rhizo_tuk <- tukey_hsd(IRdHk_G2_Rhizo, Observed ~ G1_treatment*G2_treatment)
IRH_G2_Rhizo_tuk 


### Full Figures (combined plots - Alpha diversity in Root and Rhizosphere at INRAE and MSU) ###

# G1 Alpha Diversity
G1_full_alpha <- (INRAE_G1_Root_alpha + MSU_G1_Root_alpha + plot_layout(ncol=2,widths=c(2,1))) / 
  (INRAE_G1_Rhizo_alpha + MSU_G1_Rhizo_alpha + plot_layout(ncol=2,widths=c(2,1))) + 
  plot_annotation(tag_levels = 'A')
G1_full_alpha

ggsave(plot=G1_full_alpha, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Figures/NewFigures31_07_25/Fig.4.tif",
       width=8, height=6,
       device = "tiff",
       units= "in", dpi = 400,
       compression="lzw", bg= "white")

# G2 Alpha Diversity
G2_full_alpha <- (INRAE_G2_Root_alpha + MSU_G2_Root_alpha + plot_layout(ncol=2,widths=c(2,1))) / 
  (INRAE_G2_Rhizo_alpha + MSU_G2_Rhizo_alpha + plot_layout(ncol=2,widths=c(2,1))) + plot_annotation(tag_levels = 'A')
G2_full_alpha

ggsave(plot=G2_full_alpha, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Figures/NewFigures31_07_25/Fig.5.tif",
       width=11, height=8,
       device = "tiff",
       units= "in", dpi = 400,
       compression="lzw", bg= "white")


##################################################################################################################################################

##### Beta Diversity Analysis #####


# 1. MSU Data

# Unrarefied dataset
MSU_phyloseq_unrare_cont_dro2

# Perform Multiple Rarefactions 
set.seed(13)
MSU_multrare <- phyloseq_mult_raref(MSU_phyloseq_unrare_cont_dro2,
                                   SampSize = min(sample_sums(MSU_phyloseq_unrare_cont_dro2)),
                                   replace = F,
                                   iter = 100)
sample_data(MSU_multrare$`1`)

# Save the result in the computer
#saveRDS(MSU_multrare, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/MSU_multi_rare.rds")
# load the rarefied data
MSU_multrare <- readRDS("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/MSU_multi_rare.rds")
MSU_multrare


############## MSU: Generation 1 ##############


# 1. Subset the phyloseq list to only include Generation 1 (G1)
MSU_G1_pseq <- lapply(MSU_multrare, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Generation == "G1"]
  ps_sub <- prune_samples(keep_samples,ps)
# Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
MSU_G1_pseq[[1]]

# 2. Subset the phyloseq list to only include Generation 1 (G1) Root
MSU_G1_Root_pseq <- lapply(MSU_G1_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Compartment == "root"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
MSU_G1_Root_pseq[[1]]
# Make the meta data
meta_MSU_G1_Root = data.frame(sample_data(MSU_G1_Root_pseq[[1]]))

# 3. Subset the phyloseq list to only include Generation 1 (G1) Rhizosphere
MSU_G1_Rhizo_pseq <- lapply(MSU_G1_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Compartment == "rhizosphere"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
MSU_G1_Rhizo_pseq[[1]]
# Make the meta data
meta_MSU_G1_Rhizo = data.frame(sample_data(MSU_G1_Rhizo_pseq[[1]]))

# 4. Calculate Bray-Curtis distances for each iteration

# MSU G1 Root
set.seed(13)
MSU_G1_Root_bray <- mult_dissim(MSU_G1_Root_pseq,
                                method = "bray",
                                average = F)
MSU_G1_Root_bray[[1]]
# MSU G1 Rhizo
set.seed(13)
MSU_G1_Rhizo_bray <- mult_dissim(MSU_G1_Rhizo_pseq,
                                method = "bray",
                                average = F)
MSU_G1_Rhizo_bray[[1]]

# 5. PERMANOVA test for each distance matrix

## MSU G1 Root
perm_mG1root <- how(nperm = 9999) 
set.seed(13)
MSU_G1_Root_permanova <- lapply(MSU_G1_Root_bray, function(D) {
  adonis2(D ~ G1_treatment, data = meta_MSU_G1_Root, permutations = perm_mG1root, by="term")
})
MSU_G1_Root_permanova[[100]] #NS

## MSU G1 Rhizo
perm_mG1rhizo <- how(nperm = 9999) 
set.seed(13)
MSU_G1_Rhizo_permanova <- lapply(MSU_G1_Rhizo_bray, function(D) {
  adonis2(D ~ G1_treatment, data = meta_MSU_G1_Rhizo, permutations = perm_mG1rhizo, by="term")
})
MSU_G1_Rhizo_permanova[[100]] #NS

# 6. Summarize the results

# 1) MSU G1 Root

# Extract F, R2, p-values 
MSU_G1_Root_perm_df <- do.call(rbind, lapply(MSU_G1_Root_permanova, function(x) {
  data.frame(F = x$F[1], R2 = x$R2[1], p = x$`Pr(>F)`[1])
}))

# Compute median and quantiles
MSU_G1_Root_perm.summ <- data.frame(
  F_median  = median(MSU_G1_Root_perm_df$F, na.rm = TRUE),
  R2_median = median(MSU_G1_Root_perm_df$R2, na.rm = TRUE),
  p_median  = median(MSU_G1_Root_perm_df$p, na.rm = TRUE),
  p_05      = as.numeric(quantile(MSU_G1_Root_perm_df$p, 0.05, na.rm = TRUE)),
  p_95      = as.numeric(quantile(MSU_G1_Root_perm_df$p, 0.95, na.rm = TRUE))
)

MSU_G1_Root_perm.summ # NS
#  F_median   R2_median   p_median    p_05      p_95
#  1.430279   0.1516688   0.1403      0.12939   0.15261

# 2) MSU G1 Rhizosphere

# Extract F, R2, p-values 
MSU_G1_Rhizo_perm_df <- do.call(rbind, lapply(MSU_G1_Rhizo_permanova, function(x) {
  data.frame(F = x$F[1], R2 = x$R2[1], p = x$`Pr(>F)`[1])
}))

# Compute median and quantiles
MSU_G1_Rhizo_perm.summ <- data.frame(
  F_median  = median(MSU_G1_Rhizo_perm_df$F, na.rm = TRUE),
  R2_median = median(MSU_G1_Rhizo_perm_df$R2, na.rm = TRUE),
  p_median  = median(MSU_G1_Rhizo_perm_df$p, na.rm = TRUE),
  p_05      = as.numeric(quantile(MSU_G1_Rhizo_perm_df$p, 0.05, na.rm = TRUE)),
  p_95      = as.numeric(quantile(MSU_G1_Rhizo_perm_df$p, 0.95, na.rm = TRUE))
)

MSU_G1_Rhizo_perm.summ # NS
#   F_median   R2_median   p_median   p_05       p_95
#   0.9167976  0.1028169   0.4311     0.393625   0.47629

###############################################################################

# Average Beta Diversity over multiple rarefaction iterations

# 1.) MSU G1 Root
set.seed(13)
MSU_G1_Root_average_bc <- mult_dist_average(MSU_G1_Root_bray)
MSU_G1_Root_average_bc
# 2.) MSU G1 Rhizo
set.seed(13)
MSU_G1_Rhizo_average_bc <- mult_dist_average(MSU_G1_Rhizo_bray)
MSU_G1_Rhizo_average_bc

# PERMANOVA test on the average BC distances

# 1.) MSU G1 Root
set.seed(13)
MSU_G1_Root_perm_average <- adonis2(MSU_G1_Root_average_bc ~ G1_treatment, data = meta_MSU_G1_Root, permutations = perm_mG1root, by="term")
MSU_G1_Root_perm_average # NS
#              Df SumOfSqs      R2      F Pr(>F)
#G1_treatment  1  0.13855 0.15159 1.4294 0.1398

# 2.) MSU G1 Rhizo
set.seed(13)
MSU_G1_Rhizo_perm_average <- adonis2(MSU_G1_Rhizo_average_bc ~ G1_treatment, 
                                     data = meta_MSU_G1_Rhizo, permutations = perm_mG1rhizo, by="term")
MSU_G1_Rhizo_perm_average # NS
#             Df SumOfSqs      R2      F Pr(>F)
#G1_treatment  1  0.13588 0.10282 0.9168 0.4316

# Figure PCoA Plot
 
# 1.) MSU G1 Root Ordination
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
MSU_G1_Root_cmd <- cmdscale(MSU_G1_Root_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.MSU_G1_Root <- MSU_G1_Root_cmd$points[,1]
ax2.scores.MSU_G1_Root <- MSU_G1_Root_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.MSU_G1_Root <- MSU_G1_Root_cmd$eig[1]/sum(MSU_G1_Root_cmd$eig)
ax2.MSU_G1_Root <- MSU_G1_Root_cmd$eig[2]/sum(MSU_G1_Root_cmd$eig)
meta_MSU_G1_Root_map <- cbind(meta_MSU_G1_Root,ax1.scores.MSU_G1_Root,ax2.scores.MSU_G1_Root)
# make the plot
meta_MSU_G1_Root_map$Genotype <- "Red Hawk"
MSU_G1_Root_pcoa <- ggplot(data = meta_MSU_G1_Root_map, 
                           aes(x = ax1.scores.MSU_G1_Root, y = ax2.scores.MSU_G1_Root))+
  geom_point(aes(x=ax1.scores.MSU_G1_Root, y=ax2.scores.MSU_G1_Root,color=Treatment, shape=Treatment),
             size=3, alpha=1)+
  labs(title = "(USA) G1 Root", x="PCoA1 [36.8%]", y="PCoA2 [23%]") +
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.MSU_G1_Root,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.MSU_G1_Root,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1) ) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + my_theme + 
  theme(plot.tag = element_text(face="bold")) + facet_wrap(~Genotype)
MSU_G1_Root_pcoa

# 2.) MSU G1 Rhizosphere Ordination
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
MSU_G1_Rhizo_cmd <- cmdscale(MSU_G1_Rhizo_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.MSU_G1_Rhizo <- MSU_G1_Rhizo_cmd$points[,1]
ax2.scores.MSU_G1_Rhizo <- MSU_G1_Rhizo_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.MSU_G1_Rhizo <- MSU_G1_Rhizo_cmd$eig[1]/sum(MSU_G1_Rhizo_cmd$eig)
ax2.MSU_G1_Rhizo <- MSU_G1_Rhizo_cmd$eig[2]/sum(MSU_G1_Rhizo_cmd$eig)
meta_MSU_G1_Rhizo_map <- cbind(meta_MSU_G1_Rhizo,ax1.scores.MSU_G1_Rhizo,ax2.scores.MSU_G1_Rhizo)
# make the plot
meta_MSU_G1_Rhizo_map$Genotype <- "Red Hawk"
MSU_G1_Rhizo_pcoa <- ggplot(data = meta_MSU_G1_Rhizo_map, 
                           aes(x = ax1.scores.MSU_G1_Rhizo, y = ax2.scores.MSU_G1_Rhizo))+
  geom_point(aes(x=ax1.scores.MSU_G1_Rhizo, y=ax2.scores.MSU_G1_Rhizo,color=Treatment, shape=Treatment),
             size=3, alpha=1)+
  labs(title = "(USA) G1 Rhizosphere", x="PCoA1 [32.1%]", y="PCoA2 [12%]") +
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.MSU_G1_Rhizo,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.MSU_G1_Rhizo,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1) ) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + my_theme + 
  theme(legend.position="none", plot.tag = element_text(face="bold")) + facet_wrap(~Genotype)
MSU_G1_Rhizo_pcoa


############## MSU: Generation 2 ##############


# 1. Subset the phyloseq list to only include Generation 2 (G2)
MSU_G2_pseq <- lapply(MSU_multrare, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Generation == "G2"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
MSU_G2_pseq[[1]]

# 2. Subset the phyloseq list to only include Generation 2 (G2) Root
MSU_G2_Root_pseq <- lapply(MSU_G2_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Compartment == "root"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
MSU_G2_Root_pseq[[1]]
# Make the meta data
meta_MSU_G2_Root = data.frame(sample_data(MSU_G2_Root_pseq[[1]]))

# 3. Subset the phyloseq list to only include Generation 2 (G2) Rhizosphere
MSU_G2_Rhizo_pseq <- lapply(MSU_G2_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Compartment == "rhizosphere"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
MSU_G2_Rhizo_pseq[[1]]
# Make the meta data
meta_MSU_G2_Rhizo = data.frame(sample_data(MSU_G2_Rhizo_pseq[[1]]))

# 4. Calculate Bray-Curtis distances for each iteration

# MSU G2 Root
set.seed(13)
MSU_G2_Root_bray <- mult_dissim(MSU_G2_Root_pseq,
                                method = "bray",
                                average = F)
MSU_G2_Root_bray[[1]]
# MSU G2 Rhizo
set.seed(13)
MSU_G2_Rhizo_bray <- mult_dissim(MSU_G2_Rhizo_pseq,
                                 method = "bray",
                                 average = F)
MSU_G2_Rhizo_bray[[1]]

# Average Beta Diversity over multiple rarefaction iterations

# 1.) MSU G2 Root
set.seed(13)
MSU_G2_Root_average_bc <- mult_dist_average(MSU_G2_Root_bray)
MSU_G2_Root_average_bc
# 2.) MSU G2 Rhizo
set.seed(13)
MSU_G2_Rhizo_average_bc <- mult_dist_average(MSU_G2_Rhizo_bray)
MSU_G2_Rhizo_average_bc

# PERMANOVA test on the average BC distances

# 1.) MSU G2 Root, block by planting group
perm_mG2root <- how(nperm = 9999) 
str(meta_MSU_G2_Root)
meta_MSU_G2_Root$planting_group <- as.factor(meta_MSU_G2_Root$planting_group)
setBlocks(perm_mG2root) <- with(meta_MSU_G2_Root, planting_group)

set.seed(13)
MSU_G2_Root_perm_average <- adonis2(MSU_G2_Root_average_bc ~ G2_treatment*G1_treatment, 
                                    data = meta_MSU_G2_Root, permutations = perm_mG2root, by="term")
MSU_G2_Root_perm_average # 
#                           Df SumOfSqs      R2      F Pr(>F)  
#G2_treatment               1   0.1941 0.02670 1.2517 0.0290 *
#G1_treatment               1   0.1237 0.01702 0.7980 0.2832  
#G2_treatment:G1_treatment  1   0.1284 0.01767 0.8283 0.2239  
#Residual                  44   6.8219 0.93861             

# 2.) MSU G2 Rhizo, block by planting group
perm_mG2rhizo <- how(nperm = 9999) 
str(meta_MSU_G2_Rhizo)
meta_MSU_G2_Rhizo$planting_group <- as.factor(meta_MSU_G2_Rhizo$planting_group)
setBlocks(perm_mG2rhizo) <- with(meta_MSU_G2_Rhizo, planting_group)

set.seed(13)
MSU_G2_Rhizo_perm_average <- adonis2(MSU_G2_Rhizo_average_bc ~ G2_treatment*G1_treatment, 
                                     data = meta_MSU_G2_Rhizo, permutations = perm_mG2rhizo, by="term")
MSU_G2_Rhizo_perm_average # NS

# Figure PCoA Plot

# 1.) MSU G2 Root Ordination 

# with constrained ordination (CAP = constrained analysis of principal coordinates)
# CAP with constraint and conditioning
MSU_G2_Root_cap <- capscale(MSU_G2_Root_average_bc ~ 1 + Condition(planting_group), data = meta_MSU_G2_Root)
# make the ordinates
MSU_G2_Root_cap.ord.df = data.frame(vegan::scores(MSU_G2_Root_cap, display="sites")) %>%
  tibble::rownames_to_column(var="SampleID") %>%
  dplyr::select(SampleID, MDS1, MDS2) %>%
  left_join(meta_MSU_G2_Root, by="SampleID") 

MSU_G2_Root_eigvec = vegan::eigenvals(MSU_G2_Root_cap)
MSU_G2_Root_fracvar = round(MSU_G2_Root_eigvec/sum(MSU_G2_Root_eigvec)*100, 2)
# make the plot
MSU_G2_Root_cap.ord.df$Genotype <- "Red Hawk"

MSU_G2_Root_pcoa <- ggplot(MSU_G2_Root_cap.ord.df, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=G1_G2, shape=G1_treatment), size=3, alpha=1) +
  labs(x=paste("PCoA1 [", MSU_G2_Root_fracvar[1], "%]", sep=""),
       y=paste("PCoA2 [", MSU_G2_Root_fracvar[2], "%]", sep=""), caption="G2 treatment *") +
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment",
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought"), guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  # add a title and delete the automatic caption
  labs(title = "(USA) G2 Root") + my_theme + 
  theme(aspect.ratio=1, 
        plot.tag = element_text(face="bold"),
        plot.title = element_text(size=16),
        strip.text = element_text(size = 14),
        plot.caption = element_text(size=13, face="italic",hjust=0)) + 
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Red_Hawk"="Red Hawk")))
  #annotate("text", x = 1.2, y = 1, label = "G2 treatment *", 
           #size = 4, fontface = "italic")
MSU_G2_Root_pcoa


# 2.) MSU G2 Rhizosphere Ordination 

# with constrained ordination (CAP = constrained analysis of principal coordinates)
# CAP with constraint and conditioning
MSU_G2_Rhizo_cap <- capscale(MSU_G2_Rhizo_average_bc ~ 1 + Condition(planting_group), data = meta_MSU_G2_Rhizo)
# make the ordinates
MSU_G2_Rhizo_cap.ord.df = data.frame(vegan::scores(MSU_G2_Rhizo_cap, display="sites")) %>%
  tibble::rownames_to_column(var="SampleID") %>%
  dplyr::select(SampleID, MDS1, MDS2) %>%
  left_join(meta_MSU_G2_Rhizo, by="SampleID") 

MSU_G2_Rhizo_eigvec = vegan::eigenvals(MSU_G2_Rhizo_cap)
MSU_G2_Rhizo_fracvar = round(MSU_G2_Rhizo_eigvec/sum(MSU_G2_Rhizo_eigvec)*100, 2)
# make the plot
MSU_G2_Rhizo_cap.ord.df$Genotype <- "Red Hawk"

MSU_G2_Rhizo_pcoa <- ggplot(MSU_G2_Rhizo_cap.ord.df, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=G1_G2, shape=G1_treatment), size=3, alpha=1) +
  labs(x=paste("PCoA1 [", MSU_G2_Rhizo_fracvar[1], "%]", sep=""),
       y=paste("PCoA2 [", MSU_G2_Rhizo_fracvar[2], "%]", sep="")) +
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment",
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought"), guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  # add a title and delete the automatic caption
  labs(title = "(USA) G2 Rhizosphere") + my_theme + 
  theme(aspect.ratio=1, legend.position="none", 
        plot.tag = element_text(face="bold"), 
        plot.title = element_text(size=16),
        strip.text = element_text(size = 14)) +
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Red_Hawk"="Red Hawk")))
MSU_G2_Rhizo_pcoa



# 2. INRAE Data

# Unrarefied dataset
INRAE_root_rhizo_CD

# Perform Multiple Rarefactions 
set.seed(13)
INRAE_multrare <- phyloseq_mult_raref(INRAE_root_rhizo_CD,
                                    SampSize = 2500,
                                    replace = F,
                                    iter = 100)
INRAE_multrare$`1`

# Perform Multiple Rarefactions and then do average
set.seed(13)
INRAE_multrare_average <- phyloseq_mult_raref_avg(INRAE_root_rhizo_CD,
                                      SampSize = 2500,
                                      replace = F,
                                      iter = 100)
sample_data(INRAE_multrare_average)



# Save the result in the computer
#saveRDS(INRAE_multrare, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/INRAE_multi_rare.rds")

# Read it back the multiraref data
INRAE_multrare <- readRDS("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Drought_multigeneration_study_common_bean_UPDATEDbyAbby/R_analysis_files/INRAE_multi_rare.rds")

############## INRAE: Generation 1 ##############


# 1. Subset the phyloseq list to only include Generation 1 (G1)

INR_G1_pseq <- lapply(INRAE_multrare, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Generation == "G1"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G1_pseq[[1]]

### Subset by Compartment ###

#1.)  Subset the phyloseq list to only include Generation 1 (G1) Root
INR_G1_Root_pseq <- lapply(INR_G1_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Compartment == "Root"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G1_Root_pseq[[1]]
# Make the meta data
meta_INR_G1_Root = data.frame(sample_data(INR_G1_Root_pseq[[1]]))

# 2.) Subset the phyloseq list to only include Generation 1 (G1) Rhizosphere
INR_G1_Rhizo_pseq <- lapply(INR_G1_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Compartment == "Rhizosphere"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G1_Rhizo_pseq[[1]]
# Make the meta data
meta_INR_G1_Rhizo = data.frame(sample_data(INR_G1_Rhizo_pseq[[1]]))

### Subset by Genotype ###

# 3.) Subset the phyloseq list to only include Generation 1 (G1) Root Flavert
INR_G1_Root_Flav_pseq <- lapply(INR_G1_Root_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Genotype == "Flavert"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G1_Root_Flav_pseq[[1]]
# Make the meta data
meta_INR_G1_Root_Flav = data.frame(sample_data(INR_G1_Root_Flav_pseq[[1]]))


#4.)  Subset the phyloseq list to only include Generation 1 (G1) Root Red Hawk
INR_G1_Root_RH_pseq <- lapply(INR_G1_Root_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Genotype == "Red_Hawk"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G1_Root_RH_pseq[[1]]
# Make the meta data
meta_INR_G1_Root_RH = data.frame(sample_data(INR_G1_Root_RH_pseq[[1]]))


# 5.) Subset the phyloseq list to only include Generation 1 (G1) Rhizosphere Flavert
INR_G1_Rhizo_Flav_pseq <- lapply(INR_G1_Rhizo_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Genotype == "Flavert"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G1_Rhizo_Flav_pseq[[1]]
# Make the meta data
meta_INR_G1_Rhizo_Flav = data.frame(sample_data(INR_G1_Rhizo_Flav_pseq[[1]]))


# 6.)  Subset the phyloseq list to only include Generation 1 (G1) Rhizosphere Red Hawk
INR_G1_Rhizo_RH_pseq <- lapply(INR_G1_Rhizo_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Genotype == "Red_Hawk"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G1_Rhizo_RH_pseq[[1]]
# Make the meta data
meta_INR_G1_Rhizo_RH = data.frame(sample_data(INR_G1_Rhizo_RH_pseq[[1]]))


# 4. Calculate Bray-Curtis distances for each iteration

# INRAE G1 Root
set.seed(13)
INR_G1_Root_bray <- mult_dissim(INR_G1_Root_pseq,
                                method = "bray",
                                average = F)
INR_G1_Root_bray[[1]]
# INRAE G1 Rhizosphere
set.seed(13)
INR_G1_Rhizo_bray <- mult_dissim(INR_G1_Rhizo_pseq,
                                 method = "bray",
                                 average = F)
INR_G1_Rhizo_bray[[1]]

# INRAE G1 Root - Flavert
set.seed(13)
INR_G1_Root_Flav_bray <- mult_dissim(INR_G1_Root_Flav_pseq,
                                method = "bray",
                                average = F)
# INRAE G1 Root - Red Hawk
set.seed(13)
INR_G1_Root_RH_bray <- mult_dissim(INR_G1_Root_RH_pseq,
                                     method = "bray",
                                     average = F)

# INRAE G1 Rhizosphere - Flavert
set.seed(13)
INR_G1_Rhizo_Flav_bray <- mult_dissim(INR_G1_Rhizo_Flav_pseq,
                                     method = "bray",
                                     average = F)

# INRAE G1 Rhizosphere - Red Hawk
set.seed(13)
INR_G1_Rhizo_RH_bray <- mult_dissim(INR_G1_Rhizo_RH_pseq,
                                      method = "bray",
                                      average = F)


# Average Beta Diversity over multiple rarefaction iterations

# 1.) INRAE G1 Root
set.seed(13)
INR_G1_Root_average_bc <- mult_dist_average(INR_G1_Root_bray)
INR_G1_Root_average_bc
# 2.) INRAE G1 Rhizosphere
set.seed(13)
INR_G1_Rhizo_average_bc <- mult_dist_average(INR_G1_Rhizo_bray)
INR_G1_Rhizo_average_bc
# 3.) INRAE G1 Root - Flavert
set.seed(13)
INR_G1_Root_Flav_average_bc <- mult_dist_average(INR_G1_Root_Flav_bray)
INR_G1_Root_Flav_average_bc
# 4.) INRAE G1 Root - Red Hawk
set.seed(13)
INR_G1_Root_RH_average_bc <- mult_dist_average(INR_G1_Root_RH_bray)
INR_G1_Root_RH_average_bc
# 5.) INRAE G1 Rhizosphere - Flavert
set.seed(13)
INR_G1_Rhizo_Flav_average_bc <- mult_dist_average(INR_G1_Rhizo_Flav_bray)
INR_G1_Rhizo_Flav_average_bc
# 6.) INRAE G1 Rhizosphere - Red Hawk
set.seed(13)
INR_G1_Rhizo_RH_average_bc <- mult_dist_average(INR_G1_Rhizo_RH_bray)
INR_G1_Rhizo_RH_average_bc


# PERMANOVA test on the average BC distances


# 1.)  INRAE G1 Root
perm_iG1root <- how(nperm = 9999) 
set.seed(13)
INR_G1_Root_perm_average <- adonis2(INR_G1_Root_average_bc ~ Genotype*G1_treatment, data = meta_INR_G1_Root, 
                                    permutations = perm_iG1root, by="term")
INR_G1_Root_perm_average # NS
#                     Df SumOfSqs      R2       F Pr(>F)    
#Genotype               1  0.77307 0.48604 17.0862 0.0001 ***
#G1_treatment           1  0.04857 0.03054  1.0735 0.3092    
#Genotype:G1_treatment  1  0.04499 0.02829  0.9944 0.3329    
#Residual              16  0.72393 0.45514                   

# 2.) INRAE G1 Rhizosphere
perm_iG1rhizo <- how(nperm = 9999) 
set.seed(13)
INR_G1_Rhizo_perm_average <- adonis2(INR_G1_Rhizo_average_bc ~ Genotype*G1_treatment, data = meta_INR_G1_Rhizo, 
                                     permutations = perm_iG1root, by="term")
INR_G1_Rhizo_perm_average # NS
#                       Df SumOfSqs      R2      F Pr(>F)
#Genotype               1  0.14393 0.06119 1.2232 0.2185  
#G1_treatment           1  0.20209 0.08592 1.7175 0.0241 *
#Genotype:G1_treatment  1  0.12340 0.05246 1.0487 0.2943  
#Residual              16  1.88267 0.80043                


# 3.)  INRAE G1 Root - Flavert
perm_iG1rootF <- how(nperm = 9999) 
set.seed(13)
INR_G1_Root_Flav_perm_average <- adonis2(INR_G1_Root_Flav_average_bc ~ G1_treatment, permutations = perm_iG1rootF, 
                                         data = meta_INR_G1_Root_Flav, by="term")
INR_G1_Root_Flav_perm_average #NS
# 4.)  INRAE G1 Root - Red Hawk
perm_iG1rootRH <- how(nperm = 9999) 
set.seed(13)
INR_G1_Root_RH_perm_average <- adonis2(INR_G1_Root_RH_average_bc ~ G1_treatment, permutations = perm_iG1rootRH, 
                                         data = meta_INR_G1_Root_RH, by="term")
INR_G1_Root_RH_perm_average #NS
# 5.)  INRAE G1 Rhizosphere - Flavert
perm_iG1rhizoF <- how(nperm = 9999) 
set.seed(13)
INR_G1_Rhizo_Flav_perm_average <- adonis2(INR_G1_Rhizo_Flav_average_bc ~ G1_treatment, permutations = perm_iG1rhizoF, 
                                         data = meta_INR_G1_Rhizo_Flav, by="term")
INR_G1_Rhizo_Flav_perm_average
#               Df SumOfSqs      R2      F Pr(>F)  
# G1_treatment  1  0.17049 0.16482 1.5787  0.021 *

# 6.)  INRAE G1 Rhizosphere - Red Hawk
perm_iG1rhizoRH <- how(nperm = 9999) 
set.seed(13)
INR_G1_Rhizo_RH_perm_average <- adonis2(INR_G1_Rhizo_RH_average_bc ~ G1_treatment, permutations = perm_iG1rhizoRH, 
                                       data = meta_INR_G1_Rhizo_RH, by="term")
INR_G1_Rhizo_RH_perm_average #NS


###############################################################################


# Figure PCoA Plot

# 1.) INRAE G1 Root - Flavert Ordination
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
INR_G1_Root_Flav_cmd <- cmdscale(INR_G1_Root_Flav_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.INR_G1_Root_Flav <- INR_G1_Root_Flav_cmd$points[,1]
ax2.scores.INR_G1_Root_Flav <- INR_G1_Root_Flav_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.INR_G1_Root_Flav <- INR_G1_Root_Flav_cmd$eig[1]/sum(INR_G1_Root_Flav_cmd$eig)
ax2.INR_G1_Root_Flav <- INR_G1_Root_Flav_cmd$eig[2]/sum(INR_G1_Root_Flav_cmd$eig)
meta_INR_G1_Root_Flav_map <- cbind(meta_INR_G1_Root_Flav,ax1.scores.INR_G1_Root_Flav,ax2.scores.INR_G1_Root_Flav)
# make the plot
INR_G1_Root_Flav_pcoa <- ggplot(data = meta_INR_G1_Root_Flav_map, 
                           aes(x = ax1.scores.INR_G1_Root_Flav, y = ax2.scores.INR_G1_Root_Flav))+
  geom_point(aes(x=ax1.scores.INR_G1_Root_Flav, y=ax2.scores.INR_G1_Root_Flav,color=Treatment, shape=Treatment),
             size=3, alpha=1)+
  labs(title = "(FR) G1 Root", x="PCoA1 [61.4%]", y="PCoA2 [17.2%]") +
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.INR_G1_Root_Flav,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.INR_G1_Root_Flav,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1) ) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + 
  my_theme + theme(legend.position="none", plot.tag = element_text(face="bold"))+
  facet_wrap(~ Genotype)
INR_G1_Root_Flav_pcoa

# 2.) INRAE G1 Root - Red Hawk Ordination
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
INR_G1_Root_RH_cmd <- cmdscale(INR_G1_Root_RH_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.INR_G1_Root_RH <- INR_G1_Root_RH_cmd$points[,1]
ax2.scores.INR_G1_Root_RH <- INR_G1_Root_RH_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.INR_G1_Root_RH <- INR_G1_Root_RH_cmd$eig[1]/sum(INR_G1_Root_RH_cmd$eig)
ax2.INR_G1_Root_RH <- INR_G1_Root_RH_cmd$eig[2]/sum(INR_G1_Root_RH_cmd$eig)
meta_INR_G1_Root_RH_map <- cbind(meta_INR_G1_Root_RH,ax1.scores.INR_G1_Root_RH,ax2.scores.INR_G1_Root_RH)
# make the plot
INR_G1_Root_RH_pcoa <- ggplot(data = meta_INR_G1_Root_RH_map, 
                                aes(x = ax1.scores.INR_G1_Root_RH, y = ax2.scores.INR_G1_Root_RH))+
  geom_point(aes(x=ax1.scores.INR_G1_Root_RH, y=ax2.scores.INR_G1_Root_RH,color=Treatment, shape=Treatment),
             size=3, alpha=1)+
  labs(title = "(FR) G1 Root", x="PCoA1 [40.5%]", y="PCoA2 [19.4%]") +
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.INR_G1_Root_RH,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.INR_G1_Root_RH,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1) ) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + 
  my_theme + theme(legend.position="none", plot.tag = element_text(face="bold"))+
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Red_Hawk"="Red Hawk")))
INR_G1_Root_RH_pcoa

# 3.) INRAE G1 Rhizosphere - Flavert Ordination
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
INR_G1_Rhizo_Flav_cmd <- cmdscale(INR_G1_Rhizo_Flav_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.INR_G1_Rhizo_Flav <- INR_G1_Rhizo_Flav_cmd$points[,1]
ax2.scores.INR_G1_Rhizo_Flav <- INR_G1_Rhizo_Flav_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.INR_G1_Rhizo_Flav <- INR_G1_Rhizo_Flav_cmd$eig[1]/sum(INR_G1_Rhizo_Flav_cmd$eig)
ax2.INR_G1_Rhizo_Flav <- INR_G1_Rhizo_Flav_cmd$eig[2]/sum(INR_G1_Rhizo_Flav_cmd$eig)
meta_INR_G1_Rhizo_Flav_map <- cbind(meta_INR_G1_Rhizo_Flav,ax1.scores.INR_G1_Rhizo_Flav,ax2.scores.INR_G1_Rhizo_Flav)
# make the plot
INR_G1_Rhizo_Flav_pcoa <- ggplot(data = meta_INR_G1_Rhizo_Flav_map, 
                                aes(x = ax1.scores.INR_G1_Rhizo_Flav, y = ax2.scores.INR_G1_Rhizo_Flav))+
  geom_point(aes(x=ax1.scores.INR_G1_Rhizo_Flav, y=ax2.scores.INR_G1_Rhizo_Flav,color=Treatment, shape=Treatment),
             size=3, alpha=1)+
  labs(title = "(FR) G1 Rhizosphere", x="PCoA1 [31.4%]", y="PCoA2 [11.9%]", caption="G1 treatment *") +
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.INR_G1_Rhizo_Flav,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.INR_G1_Rhizo_Flav,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1) ) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + 
  my_theme + theme(legend.position="none", 
                   plot.tag = element_text(face="bold"),
                   plot.caption = element_text(size=13, face="italic",hjust=0))+
  facet_wrap(~ Genotype)
INR_G1_Rhizo_Flav_pcoa

# 4.) INRAE G1 Rhizosphere - Red Hawk Ordination
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
INR_G1_Rhizo_RH_cmd <- cmdscale(INR_G1_Rhizo_RH_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.INR_G1_Rhizo_RH <- INR_G1_Rhizo_RH_cmd$points[,1]
ax2.scores.INR_G1_Rhizo_RH <- INR_G1_Rhizo_RH_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.INR_G1_Rhizo_RH <- INR_G1_Rhizo_RH_cmd$eig[1]/sum(INR_G1_Rhizo_RH_cmd$eig)
ax2.INR_G1_Rhizo_RH <- INR_G1_Rhizo_RH_cmd$eig[2]/sum(INR_G1_Rhizo_RH_cmd$eig)
meta_INR_G1_Rhizo_RH_map <- cbind(meta_INR_G1_Rhizo_RH,ax1.scores.INR_G1_Rhizo_RH,ax2.scores.INR_G1_Rhizo_RH)
# make the plot
INR_G1_Rhizo_RH_pcoa <- ggplot(data = meta_INR_G1_Rhizo_RH_map, 
                              aes(x = ax1.scores.INR_G1_Rhizo_RH, y = ax2.scores.INR_G1_Rhizo_RH))+
  geom_point(aes(x=ax1.scores.INR_G1_Rhizo_RH, y=ax2.scores.INR_G1_Rhizo_RH,color=Treatment, shape=Treatment),
             size=3, alpha=1)+
  labs(title = "(FR) G1 Rhizosphere", x="PCoA1 [32.5%]", y="PCoA2 [14.1%]") +
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.INR_G1_Rhizo_RH,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.INR_G1_Rhizo_RH,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1 Treatment",
                     values = c("#3399FF","#FFCC00"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1) ) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + 
  my_theme + theme(legend.position="none", plot.tag = element_text(face="bold"))+
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Red_Hawk"="Red Hawk")))
INR_G1_Rhizo_RH_pcoa



############## INRAE: Generation 2 ##############


# 1. Subset the phyloseq list to only include Generation 2 (G2)

INR_G2_pseq <- lapply(INRAE_multrare, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Generation == "G2"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G2_pseq[[1]]

### Subset by Compartment ###

#1.)  Subset the phyloseq list to only include Generation 2 (G2) Root
INR_G2_Root_pseq <- lapply(INR_G2_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Compartment == "Root"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G2_Root_pseq[[1]]
# Make the meta data
meta_INR_G2_Root = data.frame(sample_data(INR_G2_Root_pseq[[1]]))

# 2.) Subset the phyloseq list to only include Generation 2 (G2) Rhizosphere
INR_G2_Rhizo_pseq <- lapply(INR_G2_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Compartment == "Rhizosphere"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G2_Rhizo_pseq[[1]]
# Make the meta data
meta_INR_G2_Rhizo = data.frame(sample_data(INR_G2_Rhizo_pseq[[1]]))

### Subset by Genotype ###

# 3.) Subset the phyloseq list to only include Generation 2 (G2) Root Flavert
INR_G2_Root_Flav_pseq <- lapply(INR_G2_Root_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Genotype == "Flavert"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G2_Root_Flav_pseq[[1]]
# Make the meta data
meta_INR_G2_Root_Flav = data.frame(sample_data(INR_G2_Root_Flav_pseq[[1]]))


#4.)  Subset the phyloseq list to only include Generation 2 (G2) Root Red Hawk
INR_G2_Root_RH_pseq <- lapply(INR_G2_Root_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Genotype == "Red_Hawk"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G2_Root_RH_pseq[[1]]
# Make the meta data
meta_INR_G2_Root_RH = data.frame(sample_data(INR_G2_Root_RH_pseq[[1]]))


# 5.) Subset the phyloseq list to only include Generation 2 (G2) Rhizosphere Flavert
INR_G2_Rhizo_Flav_pseq <- lapply(INR_G2_Rhizo_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Genotype == "Flavert"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G2_Rhizo_Flav_pseq[[1]]
# Make the meta data
meta_INR_G2_Rhizo_Flav = data.frame(sample_data(INR_G2_Rhizo_Flav_pseq[[1]]))


# 6.)  Subset the phyloseq list to only include Generation 2 (G2) Rhizosphere Red Hawk
INR_G2_Rhizo_RH_pseq <- lapply(INR_G2_Rhizo_pseq, function(ps) {
  keep_samples <- sample_names(ps)[sample_data(ps)$Genotype == "Red_Hawk"]
  ps_sub <- prune_samples(keep_samples,ps)
  # Remove taxa that are now absent in all remaining samples
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  return(ps_sub)
})
# Check the first subsetted object
INR_G2_Rhizo_RH_pseq[[1]]
# Make the meta data
meta_INR_G2_Rhizo_RH = data.frame(sample_data(INR_G2_Rhizo_RH_pseq[[1]]))


# 4. Calculate Bray-Curtis distances for each iteration

# INRAE G2 Root
set.seed(13)
INR_G2_Root_bray <- mult_dissim(INR_G2_Root_pseq,
                                method = "bray",
                                average = F)
INR_G2_Root_bray[[1]]
# INRAE G2 Rhizosphere
set.seed(13)
INR_G2_Rhizo_bray <- mult_dissim(INR_G2_Rhizo_pseq,
                                 method = "bray",
                                 average = F)
INR_G2_Rhizo_bray[[1]]

# INRAE G2 Root - Flavert
set.seed(13)
INR_G2_Root_Flav_bray <- mult_dissim(INR_G2_Root_Flav_pseq,
                                     method = "bray",
                                     average = F)
# INRAE G2 Root - Red Hawk
set.seed(13)
INR_G2_Root_RH_bray <- mult_dissim(INR_G2_Root_RH_pseq,
                                   method = "bray",
                                   average = F)

# INRAE G2 Rhizosphere - Flavert
set.seed(13)
INR_G2_Rhizo_Flav_bray <- mult_dissim(INR_G2_Rhizo_Flav_pseq,
                                      method = "bray",
                                      average = F)

# INRAE G2 Rhizosphere - Red Hawk
set.seed(13)
INR_G2_Rhizo_RH_bray <- mult_dissim(INR_G2_Rhizo_RH_pseq,
                                    method = "bray",
                                    average = F)


# Average Beta Diversity over multiple rarefaction iterations

# 1.) INRAE G2 Root
set.seed(13)
INR_G2_Root_average_bc <- mult_dist_average(INR_G2_Root_bray)
INR_G2_Root_average_bc
# 2.) INRAE G2 Rhizosphere
set.seed(13)
INR_G2_Rhizo_average_bc <- mult_dist_average(INR_G2_Rhizo_bray)
INR_G2_Rhizo_average_bc
# 3.) INRAE G2 Root - Flavert
set.seed(13)
INR_G2_Root_Flav_average_bc <- mult_dist_average(INR_G2_Root_Flav_bray)
INR_G2_Root_Flav_average_bc
# 4.) INRAE G2 Root - Red Hawk
set.seed(13)
INR_G2_Root_RH_average_bc <- mult_dist_average(INR_G2_Root_RH_bray)
INR_G2_Root_RH_average_bc
# 5.) INRAE G2 Rhizosphere - Flavert
set.seed(13)
INR_G2_Rhizo_Flav_average_bc <- mult_dist_average(INR_G2_Rhizo_Flav_bray)
INR_G2_Rhizo_Flav_average_bc
# 6.) INRAE G2 Rhizosphere - Red Hawk
set.seed(13)
INR_G2_Rhizo_RH_average_bc <- mult_dist_average(INR_G2_Rhizo_RH_bray)
INR_G2_Rhizo_RH_average_bc


# PERMANOVA test on the average BC distances


# 1.)  INRAE G2 Root
perm_iG2root <- how(nperm = 9999) 
set.seed(13)
INR_G2_Root_perm_average <- adonis2(INR_G2_Root_average_bc ~ Genotype*G1_treatment*G2_treatment, 
                                    data = meta_INR_G2_Root, permutations = perm_iG2root, by="term")
INR_G2_Root_perm_average 
# Genotype ***
# G1:G2 *

# 2.) INRAE G2 Rhizosphere
perm_iG2rhizo <- how(nperm = 9999) 
set.seed(13)
INR_G2_Rhizo_perm_average <- adonis2(INR_G2_Rhizo_average_bc ~ Genotype*G1_treatment*G2_treatment, 
                                     data = meta_INR_G2_Rhizo, permutations = perm_iG2root, by="term")
INR_G2_Rhizo_perm_average # 
#                                   Df SumOfSqs      R2      F Pr(>F)    
#Genotype                            1   0.3732 0.06125 2.5661 0.0007 ***
#G1_treatment                        1   0.4197 0.06888 2.8856 0.0002 ***
#G2_treatment                        1   0.3341 0.05483 2.2972 0.0031 ** 
#Genotype:G1_treatment               1   0.3997 0.06560 2.7483 0.0006 ***              


# 3.)  INRAE G2 Root - Flavert
perm_iG2rootF <- how(nperm = 9999) 
set.seed(13)
INR_G2_Root_Flav_perm_average <- adonis2(INR_G2_Root_Flav_average_bc ~ G1_treatment*G2_treatment, permutations = perm_iG2rootF, 
                                         data = meta_INR_G2_Root_Flav, by="term")
INR_G2_Root_Flav_perm_average 
#                           Df SumOfSqs      R2      F Pr(>F)
#G1_treatment:G2_treatment  1  0.46867 0.24490 4.9412 0.0272 *

# Post-hoc
library(pairwiseAdonis)
set.seed(13)
INR_G2_Root_Flav_pairwise <- pairwise.adonis2(INR_G2_Root_Flav_average_bc ~ G1_G2, 
                                              data = meta_INR_G2_Root_Flav, permutations = perm_iG2rootF,by="terms")
INR_G2_Root_Flav_pairwise
# Control_Drought_vs_Drought_Drought  G1_G2 1  0.40364 0.42155 5.1013  0.018 *
# Drought_Control_vs_Drought_Drought  G1_G2 1  0.37498 0.33122 3.962 0.0724 .

# 4.)  INRAE G2 Root - Red Hawk
perm_iG2rootRH <- how(nperm = 9999) 
set.seed(13)
INR_G2_Root_RH_perm_average <- adonis2(INR_G2_Root_RH_average_bc ~ G1_treatment*G2_treatment, permutations = perm_iG2rootRH, 
                                       data = meta_INR_G2_Root_RH, by="term")
INR_G2_Root_RH_perm_average # NS

# Post-hoc
set.seed(13)
INR_G2_Root_RH_pairwise <- pairwise.adonis2(INR_G2_Root_RH_average_bc ~ G1_G2, 
                                             data = meta_INR_G2_Root_RH, permutations = 9999, by="terms")
INR_G2_Root_RH_pairwise #NS


# 5.)  INRAE G2 Rhizosphere - Flavert
perm_iG2rhizoF <- how(nperm = 9999) 
set.seed(13)
INR_G2_Rhizo_Flav_perm_average <- adonis2(INR_G2_Rhizo_Flav_average_bc ~ G1_treatment*G2_treatment, permutations = perm_iG2rhizoF, 
                                          data = meta_INR_G2_Rhizo_Flav, by="term")
INR_G2_Rhizo_Flav_perm_average
#                           Df SumOfSqs      R2      F Pr(>F)    
#G1_treatment               1   0.6704 0.20094 4.3202 0.0001 ***
#G2_treatment               1   0.2530 0.07582 1.6302 0.0497 *  
#G1_treatment:G2_treatment  1   0.2404 0.07205 1.5491 0.0536 .

# Post-hoc
set.seed(13)
INR_G2_Rhizo_Flav_pairwise <- pairwise.adonis2(INR_G2_Rhizo_Flav_average_bc ~ G1_G2, 
                                              data = meta_INR_G2_Rhizo_Flav, permutations = perm_iG2rhizoF,by="terms")
INR_G2_Rhizo_Flav_pairwise
# Control_Control_vs_Control_Drought 1  0.39623 0.2332 2.1289 0.0076 **
# Control_Control_vs_Drought_Drought 1  0.39672 0.33714 4.0688  0.007 **
# Control_Control_vs_Drought_Control 1  0.39467 0.38857 4.4485 0.0079 **
# Control_Drought_vs_Drought_Drought 1   0.5306 0.25483 2.3938 0.0084 **
# Control_Drought_vs_Drought_Control 1  0.52669 0.27442 2.2693 0.0291 *

# 6.)  INRAE G2 Rhizosphere - Red Hawk
perm_iG2rhizoRH <- how(nperm = 9999) 
set.seed(13)
INR_G2_Rhizo_RH_perm_average <- adonis2(INR_G2_Rhizo_RH_average_bc ~ G1_treatment*G2_treatment, permutations = perm_iG2rhizoRH, 
                                        data = meta_INR_G2_Rhizo_RH, by="term")
INR_G2_Rhizo_RH_perm_average 
#                           Df SumOfSqs      R2      F Pr(>F)  
#G1_treatment               1  0.13509 0.05668 0.9956 0.4454  
#G2_treatment               1  0.21431 0.08992 1.5795 0.0525 .
#G1_treatment:G2_treatment  1  0.13431 0.05635 0.9898 0.4634  

# Post-hoc
set.seed(13)
INR_G2_Rhizo_RH_pairwise <- pairwise.adonis2(INR_G2_Rhizo_RH_average_bc ~ G1_G2, 
                                               data = meta_INR_G2_Rhizo_RH, permutations = perm_iG2rhizoRH,by="terms")
INR_G2_Rhizo_RH_pairwise #NS

# Figure PCoA Plot

# 1.) INRAE G2 Root - Flavert Ordination

# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
INR_G2_Root_Flav_cmd <- cmdscale(INR_G2_Root_Flav_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.INR_G2_Root_Flav <- INR_G2_Root_Flav_cmd$points[,1]
ax2.scores.INR_G2_Root_Flav <- INR_G2_Root_Flav_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.INR_G2_Root_Flav <- INR_G2_Root_Flav_cmd$eig[1]/sum(INR_G2_Root_Flav_cmd$eig)
ax2.INR_G2_Root_Flav <- INR_G2_Root_Flav_cmd$eig[2]/sum(INR_G2_Root_Flav_cmd$eig)
meta_INR_G2_Root_Flav_map <- cbind(meta_INR_G2_Root_Flav,ax1.scores.INR_G2_Root_Flav,ax2.scores.INR_G2_Root_Flav)
# make the plot
INR_G2_Root_Flav_pcoa <- ggplot(data = meta_INR_G2_Root_Flav_map, 
                        aes(x = ax1.scores.INR_G2_Root_Flav, y = ax2.scores.INR_G2_Root_Flav))+   
  geom_point(aes(x=ax1.scores.INR_G2_Root_Flav, y=ax2.scores.INR_G2_Root_Flav,color=G1_G2, shape=G1_treatment),
             size=3, alpha=1)+
  labs(title = "(FR) G2 Root", x="PCoA1 [82%]", y="PCoA2 [10.8%]", caption="G1 x G2 treatment *")+
  #caption = "C_C vs D_C **    D_C vs D_D *    C_D vs D_D p=0.06") +
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.INR_G2_Root_Flav,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.INR_G2_Root_Flav,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment",
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"),
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + my_theme + 
  theme(legend.position="none", 
        plot.tag = element_text(face="bold"),
        plot.title = element_text(size=16),
        strip.text = element_text(size = 14),
        plot.caption = element_text(size=13, face="italic",hjust=0)) + 
  facet_wrap(~Genotype)
  #annotate("text", x = -0.3, y = 0.42, label = "G1 x G2 treatment *", 
           #size = 4, fontface = "italic")

INR_G2_Root_Flav_pcoa

# 2.) INRAE G2 Root - Red Hawk Ordination

# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
INR_G2_Root_RH_cmd <- cmdscale(INR_G2_Root_RH_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.INR_G2_Root_RH <- INR_G2_Root_RH_cmd$points[,1]
ax2.scores.INR_G2_Root_RH <- INR_G2_Root_RH_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.INR_G2_Root_RH <- INR_G2_Root_RH_cmd$eig[1]/sum(INR_G2_Root_RH_cmd$eig)
ax2.INR_G2_Root_RH <- INR_G2_Root_RH_cmd$eig[2]/sum(INR_G2_Root_RH_cmd$eig)
meta_INR_G2_Root_RH_map <- cbind(meta_INR_G2_Root_RH,ax1.scores.INR_G2_Root_RH,ax2.scores.INR_G2_Root_RH)
# make the plot
INR_G2_Root_RH_pcoa <- ggplot(data = meta_INR_G2_Root_RH_map, 
                                aes(x = ax1.scores.INR_G2_Root_RH, y = ax2.scores.INR_G2_Root_RH))+   
  geom_point(aes(x=ax1.scores.INR_G2_Root_RH, y=ax2.scores.INR_G2_Root_RH,color=G1_G2, shape=G1_treatment),
             size=3, alpha=1)+
  labs(title = "(FR) G2 Root", x="PCoA1 [86.5%]", y="PCoA2 [6.4%]") +
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.INR_G2_Root_RH,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.INR_G2_Root_RH,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment",
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"),
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), 
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + 
  my_theme + theme(legend.position="none", 
                   plot.tag = element_text(face="bold"),
                   plot.title = element_text(size=16),
                   strip.text = element_text(size = 14)) + 
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Red_Hawk"="Red Hawk")))
INR_G2_Root_RH_pcoa

# 3.) INRAE G2 Rhizosphere - Flavert Ordination

# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
INR_G2_Rhizo_Flav_cmd <- cmdscale(INR_G2_Rhizo_Flav_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.INR_G2_Rhizo_Flav <- INR_G2_Rhizo_Flav_cmd$points[,1]
ax2.scores.INR_G2_Rhizo_Flav <- INR_G2_Rhizo_Flav_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.INR_G2_Rhizo_Flav <- INR_G2_Rhizo_Flav_cmd$eig[1]/sum(INR_G2_Rhizo_Flav_cmd$eig)
ax2.INR_G2_Rhizo_Flav <- INR_G2_Rhizo_Flav_cmd$eig[2]/sum(INR_G2_Rhizo_Flav_cmd$eig)
meta_INR_G2_Rhizo_Flav_map <- cbind(meta_INR_G2_Rhizo_Flav,ax1.scores.INR_G2_Rhizo_Flav,ax2.scores.INR_G2_Rhizo_Flav)
# make the plot
INR_G2_Rhizo_Flav_pcoa <- ggplot(data = meta_INR_G2_Rhizo_Flav_map, 
                                aes(x = ax1.scores.INR_G2_Rhizo_Flav, y = ax2.scores.INR_G2_Rhizo_Flav))+   
  geom_point(aes(x=ax1.scores.INR_G2_Rhizo_Flav, y=ax2.scores.INR_G2_Rhizo_Flav,color=G1_G2, shape=G1_treatment),
             size=3, alpha=1)+
  labs(title = "(FR) G2 Rhizosphere", x="PCoA1 [25.3%]", y="PCoA2 [19.6%]", caption = "G1 treatment ***\nG2 treatment *")+
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.INR_G2_Rhizo_Flav,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.INR_G2_Rhizo_Flav,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment",
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought"), guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"),
                     guide = guide_legend(override.aes = list(size = 3), alpha = 1)) + 
  my_theme + 
  theme(legend.position="none", 
        plot.tag = element_text(face="bold"),
        plot.title = element_text(size=16),
        strip.text = element_text(size = 14),
        plot.caption = element_text(size=13, face="italic",hjust=0)) + 
  facet_wrap(~Genotype)
  #annotate("text", x = -0.7, y = 0.25, label = "G1 treatment ***\nG2 treatment *", 
           #size = 4, fontface = "italic", hjust=0, vjust=1)

INR_G2_Rhizo_Flav_pcoa


# 4.) INRAE G2 Rhizosphere - Red Hawk Ordination

# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
INR_G2_Rhizo_RH_cmd <- cmdscale(INR_G2_Rhizo_RH_average_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.INR_G2_Rhizo_RH <- INR_G2_Rhizo_RH_cmd$points[,1]
ax2.scores.INR_G2_Rhizo_RH <- INR_G2_Rhizo_RH_cmd$points[,2] 
# calculate percent variance explained, then add to plot
ax1.INR_G2_Rhizo_RH <- INR_G2_Rhizo_RH_cmd$eig[1]/sum(INR_G2_Rhizo_RH_cmd$eig)
ax2.INR_G2_Rhizo_RH <- INR_G2_Rhizo_RH_cmd$eig[2]/sum(INR_G2_Rhizo_RH_cmd$eig)
meta_INR_G2_Rhizo_RH_map <- cbind(meta_INR_G2_Rhizo_RH,ax1.scores.INR_G2_Rhizo_RH,ax2.scores.INR_G2_Rhizo_RH)
# make the plot
INR_G2_Rhizo_RH_pcoa <- ggplot(data = meta_INR_G2_Rhizo_RH_map, 
                                 aes(x = ax1.scores.INR_G2_Rhizo_RH, y = ax2.scores.INR_G2_Rhizo_RH))+   
  geom_point(aes(x=ax1.scores.INR_G2_Rhizo_RH, y=ax2.scores.INR_G2_Rhizo_RH,color=G1_G2, shape=G1_treatment),
             size=3, alpha=1)+
  labs(title = "(FR) G2 Rhizosphere", x="PCoA1 [22%]", y="PCoA2 [13.2%]")+
  #scale_x_continuous(name = paste("PCoA1: ",round(ax1.INR_G2_Rhizo_RH,3)*100,"% var. explained", sep=""))+
  #scale_y_continuous(name = paste("PCoA2: ",round(ax2.INR_G2_Rhizo_RH,3)*100,"% var. explained", sep=""))+
  scale_color_manual(name = "G1_G2 Consecutive\nTreatment",
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought"), guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  scale_shape_manual(name = "G1 Treatment",
                     values = c("circle", "triangle"), labels = c("Control", "Drought"), guide = guide_legend(override.aes = list(size = 3), alpha = 1)) +
  my_theme + theme(plot.tag = element_text(face="bold"),
                   legend.position="none", 
                   plot.title = element_text(size=16),
                   strip.text = element_text(size = 14)) + 
  facet_wrap(~ Genotype, nrow=1, scales="free", labeller = labeller(Genotype = c("Red_Hawk"="Red Hawk")))

INR_G2_Rhizo_RH_pcoa


### Full Figures (combined plots - Alpha diversity in Root and Rhizosphere at INRAE and MSU) ###

# G1 Beta Diversity

G1_full_beta <- (INR_G1_Root_Flav_pcoa + INR_G1_Root_RH_pcoa + MSU_G1_Root_pcoa) / 
                (INR_G1_Rhizo_Flav_pcoa + INR_G1_Rhizo_RH_pcoa + MSU_G1_Rhizo_pcoa) + 
                plot_annotation(tag_levels = 'A')

G1_full_beta

ggsave(plot=G1_full_beta, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Figures/NewFigures31_07_25/Fig.6.tif",
       width=9, height=6,
       device = "tiff",
       units= "in", dpi = 400,
       compression="lzw", bg= "white")

# G2 Beta Diversity

G2_full_beta <- (INR_G2_Root_Flav_pcoa + INR_G2_Root_RH_pcoa + MSU_G2_Root_pcoa) / 
                (INR_G2_Rhizo_Flav_pcoa + INR_G2_Rhizo_RH_pcoa + MSU_G2_Rhizo_pcoa) + 
                plot_annotation(tag_levels = 'A') &
                theme(plot.tag=element_text(size = 20))

G2_full_beta


ggsave(plot=G2_full_beta, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Figures/NewFigures31_07_25/Fig.7.tif",
       width=12, height=8,
       device = "tiff",
       units= "in", dpi = 400,
       compression="lzw", bg= "white")


##### Beta Dispersion Analysis #####


### MSU: Generation 1 ###

# 1. MSU G1 Root
set.seed(13)
MSU_G1_Root_betadisp <- betadisper(MSU_G1_Root_average_bc, 
                                   group=as.factor(meta_MSU_G1_Root$G1_treatment), type = "median")
MSU_G1_Root_betadisp
set.seed(13)
permutest(MSU_G1_Root_betadisp, pairwise=TRUE, permutations=999) #NS
TukeyHSD(MSU_G1_Root_betadisp)#NS
# boxplot
boxplot(MSU_G1_Root_betadisp)
# get betadisper data
MSU_G1_Root_BDdata <- get_betadisper_data(MSU_G1_Root_betadisp)
# do some transformations on the data
MSU_G1_Root_BDdata$eigenvalue <- mutate(MSU_G1_Root_BDdata$eigenvalue, percent = eig/sum(eig))
# add convex hull points 
# this could be put in a function
MSU_G1_Root_BDdata$chull <- group_by(MSU_G1_Root_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# combine centroid and eigenvector dataframes for plotting (only need this if we want to plot the centroids)
MSU_G1_Root_BDdata_lines <- merge(dplyr::select(MSU_G1_Root_BDdata$centroids, group, PCoA1, PCoA2), 
                                  dplyr::select(MSU_G1_Root_BDdata$eigenvector, group, PCoA1, PCoA2), 
                                  by = c('group'))
# Now the dataframes are all ready to be completely customisable in ggplot
MSU_G1_Root_betadisp_plot <- ggplot(MSU_G1_Root_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme +
  scale_y_continuous(limits = c(0.2, 0.4))+
  scale_color_manual(name="G1 Treatment", values = c("#3399FF","#FFCC00"), labels = c("Control", 'Drought')) +
  ylab('Distance to Median') + xlab("G1 Treatment") +
  scale_x_discrete(labels = c("Control", 'Drought')) +
  theme(legend.position = "none")
MSU_G1_Root_betadisp_plot 

# 2. MSU G1 Rhizosphere
set.seed(13)
MSU_G1_Rhizo_betadisp <- betadisper(MSU_G1_Rhizo_average_bc, 
                                   group=as.factor(meta_MSU_G1_Rhizo$G1_treatment), type = "median")
MSU_G1_Rhizo_betadisp
set.seed(13)
permutest(MSU_G1_Rhizo_betadisp, pairwise=TRUE, permutations=999) #NS
TukeyHSD(MSU_G1_Rhizo_betadisp)#NS
# boxplot
boxplot(MSU_G1_Rhizo_betadisp)
# get betadisper data
MSU_G1_Rhizo_BDdata <- get_betadisper_data(MSU_G1_Rhizo_betadisp)
# do some transformations on the data
MSU_G1_Rhizo_BDdata$eigenvalue <- mutate(MSU_G1_Rhizo_BDdata$eigenvalue, percent = eig/sum(eig))
# add convex hull points 
# this could be put in a function
MSU_G1_Rhizo_BDdata$chull <- group_by(MSU_G1_Rhizo_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# combine centroid and eigenvector dataframes for plotting (only need this if we want to plot the centroids)
MSU_G1_Rhizo_BDdata_lines <- merge(dplyr::select(MSU_G1_Rhizo_BDdata$centroids, group, PCoA1, PCoA2), 
                                  dplyr::select(MSU_G1_Rhizo_BDdata$eigenvector, group, PCoA1, PCoA2), 
                                  by = c('group'))
# Now the dataframes are all ready to be completely customisable in ggplot
MSU_G1_Rhizo_betadisp_plot <- ggplot(MSU_G1_Rhizo_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme + scale_y_continuous(limits = c(0.2, 0.4))+
  scale_color_manual(name="G1 Treatment", values = c("#3399FF","#FFCC00"), labels = c("Control", 'Drought')) +
  ylab('Distance to Median') + xlab("G1 Treatment") +
  scale_x_discrete(labels = c("Control", 'Drought')) +
  theme(legend.position = "none")
MSU_G1_Rhizo_betadisp_plot 


### INRAE: Generation 1 ###

# 1. INRAE G1 Root - Flavert
INR_G1_Root_Flav_betadisp <- betadisper(INR_G1_Root_Flav_average_bc, 
                                   group=as.factor(meta_INR_G1_Root_Flav$G1_treatment), type = "median")
INR_G1_Root_Flav_betadisp
set.seed(13)
permutest(INR_G1_Root_Flav_betadisp, pairwise=TRUE, permutations=999) #NS
TukeyHSD(INR_G1_Root_Flav_betadisp)#NS
# boxplot
boxplot(INR_G1_Root_Flav_betadisp)
# get betadisper data
INR_G1_Root_Flav_BDdata <- get_betadisper_data(INR_G1_Root_Flav_betadisp)
# do some transformations on the data
INR_G1_Root_Flav_BDdata$eigenvalue <- mutate(INR_G1_Root_Flav_BDdata$eigenvalue, percent = eig/sum(eig))
# add convex hull points 
# this could be put in a function
INR_G1_Root_Flav_BDdata$chull <- group_by(INR_G1_Root_Flav_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# combine centroid and eigenvector dataframes for plotting (only need this if we want to plot the centroids)
INR_G1_Root_Flav_BDdata_lines <- merge(dplyr::select(INR_G1_Root_Flav_BDdata$centroids, group, PCoA1, PCoA2), 
                                  dplyr::select(INR_G1_Root_Flav_BDdata$eigenvector, group, PCoA1, PCoA2), 
                                  by = c('group'))
# Now the dataframes are all ready to be completely customisable in ggplot
INR_G1_Root_Flav_betadisp_plot <- ggplot(INR_G1_Root_Flav_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme +
  scale_color_manual(name="G1 Treatment", values = c("#3399FF","#FFCC00"), labels = c("Control", 'Drought')) +
  ylab('Distance to Median') + xlab("G1 Treatment") +
  scale_x_discrete(labels = c("Control", 'Drought')) +
  theme(legend.position = "none")
INR_G1_Root_Flav_betadisp_plot 

# 2. INRAE G1 Root - Red Hawk
INR_G1_Root_RH_betadisp <- betadisper(INR_G1_Root_RH_average_bc, 
                                        group=as.factor(meta_INR_G1_Root_RH$G1_treatment), type = "median")
INR_G1_Root_RH_betadisp
set.seed(13)
permutest(INR_G1_Root_RH_betadisp, pairwise=TRUE, permutations=999) #NS
TukeyHSD(INR_G1_Root_RH_betadisp)#NS
# boxplot
boxplot(INR_G1_Root_RH_betadisp)
# get betadisper data
INR_G1_Root_RH_BDdata <- get_betadisper_data(INR_G1_Root_RH_betadisp)
# do some transformations on the data
INR_G1_Root_RH_BDdata$eigenvalue <- mutate(INR_G1_Root_RH_BDdata$eigenvalue, percent = eig/sum(eig))
# add convex hull points 
# this could be put in a function
INR_G1_Root_RH_BDdata$chull <- group_by(INR_G1_Root_RH_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# combine centroid and eigenvector dataframes for plotting (only need this if we want to plot the centroids)
INR_G1_Root_RH_BDdata_lines <- merge(dplyr::select(INR_G1_Root_RH_BDdata$centroids, group, PCoA1, PCoA2), 
                                       dplyr::select(INR_G1_Root_RH_BDdata$eigenvector, group, PCoA1, PCoA2), 
                                       by = c('group'))
# Now the dataframes are all ready to be completely customisable in ggplot
INR_G1_Root_RH_betadisp_plot <- ggplot(INR_G1_Root_RH_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme +
  scale_color_manual(name="G1 Treatment", values = c("#3399FF","#FFCC00"), labels = c("Control", 'Drought')) +
  ylab('Distance to Median') + xlab("G1 Treatment") +
  scale_x_discrete(labels = c("Control", 'Drought')) +
  theme(legend.position = "none")
INR_G1_Root_RH_betadisp_plot 

# 3. INRAE G1 Rhizosphere - Flavert
INR_G1_Rhizo_Flav_betadisp <- betadisper(INR_G1_Rhizo_Flav_average_bc, 
                                        group=as.factor(meta_INR_G1_Rhizo_Flav$G1_treatment), type = "median")
INR_G1_Rhizo_Flav_betadisp
set.seed(13)
permutest(INR_G1_Rhizo_Flav_betadisp, pairwise=TRUE, permutations=999) 
INR_G1_Rhizo_Flav_permut <- permutest(INR_G1_Rhizo_Flav_betadisp, pairwise=TRUE, permutations=999) 
set.seed(13)
p.adjust(INR_G1_Rhizo_Flav_permut$pairwise$permuted, method = "fdr") #      0.038 
TukeyHSD(INR_G1_Rhizo_Flav_betadisp)# p-val=0.08
#          Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups  1 0.0060854 0.0060854 4.0062    999  0.046 *
#Pairwise comparisons:
# Control vs. Drought           0.034 (permuted p-value)

# boxplot
boxplot(INR_G1_Rhizo_Flav_betadisp)
# get betadisper data
INR_G1_Rhizo_Flav_BDdata <- get_betadisper_data(INR_G1_Rhizo_Flav_betadisp)
# do some transformations on the data
INR_G1_Rhizo_Flav_BDdata$eigenvalue <- mutate(INR_G1_Rhizo_Flav_BDdata$eigenvalue, percent = eig/sum(eig))
# add convex hull points 
# this could be put in a function
INR_G1_Rhizo_Flav_BDdata$chull <- group_by(INR_G1_Rhizo_Flav_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# combine centroid and eigenvector dataframes for plotting (only need this if we want to plot the centroids)
INR_G1_Rhizo_Flav_BDdata_lines <- merge(dplyr::select(INR_G1_Rhizo_Flav_BDdata$centroids, group, PCoA1, PCoA2), 
                                       dplyr::select(INR_G1_Rhizo_Flav_BDdata$eigenvector, group, PCoA1, PCoA2), 
                                       by = c('group'))
# Now the dataframes are all ready to be completely customisable in ggplot
iG1rhizoF_BD_anno <- data.frame(x1 = c(1), x2 = c(2), 
                                y1 = c(0.42), y2 = c(0.43),
                                xstar = c(1.5), ystar = c(0.433),
                                lab = c("*"))

INR_G1_Rhizo_Flav_betadisp_plot <- ggplot(INR_G1_Rhizo_Flav_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme + labs(title="(FR) G1 Rhizosphere - Flavert") +
  labs(caption="G1 treatment *") + 
  theme(plot.tag = element_text(face="bold"),
        plot.caption = element_text(size=13, face="italic", hjust=0),
        legend.position = "none", 
        axis.title.x = element_blank())+
  scale_color_manual(name="G1 Treatment", values = c("#3399FF","#FFCC00"), labels = c("Control", 'Drought')) +
  ylab('Distance to Median') + xlab("G1 Treatment") +
  scale_x_discrete(labels = c("Control", 'Drought')) +
  geom_text(inherit.aes=FALSE, data = iG1rhizoF_BD_anno, aes(x = xstar,  y = ystar, label = lab), size=6) +
  geom_segment(inherit.aes=FALSE, data =iG1rhizoF_BD_anno, aes(x = x1, xend = x1, y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes=FALSE, data = iG1rhizoF_BD_anno, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes=FALSE, data = iG1rhizoF_BD_anno, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black")
INR_G1_Rhizo_Flav_betadisp_plot 

# 4. INRAE G1 Rhizosphere - Red Hawk
INR_G1_Rhizo_RH_betadisp <- betadisper(INR_G1_Rhizo_RH_average_bc, 
                                      group=as.factor(meta_INR_G1_Rhizo_RH$G1_treatment), type = "median")
INR_G1_Rhizo_RH_betadisp
set.seed(13)
permutest(INR_G1_Rhizo_RH_betadisp, pairwise=TRUE, permutations=999) #NS
TukeyHSD(INR_G1_Rhizo_RH_betadisp)#NS
# boxplot
boxplot(INR_G1_Rhizo_RH_betadisp)
# get betadisper data
INR_G1_Rhizo_RH_BDdata <- get_betadisper_data(INR_G1_Rhizo_RH_betadisp)
# do some transformations on the data
INR_G1_Rhizo_RH_BDdata$eigenvalue <- mutate(INR_G1_Rhizo_RH_BDdata$eigenvalue, percent = eig/sum(eig))
# add convex hull points 
# this could be put in a function
INR_G1_Rhizo_RH_BDdata$chull <- group_by(INR_G1_Rhizo_RH_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# combine centroid and eigenvector dataframes for plotting (only need this if we want to plot the centroids)
INR_G1_Rhizo_RH_BDdata_lines <- merge(dplyr::select(INR_G1_Rhizo_RH_BDdata$centroids, group, PCoA1, PCoA2), 
                                     dplyr::select(INR_G1_Rhizo_RH_BDdata$eigenvector, group, PCoA1, PCoA2), 
                                     by = c('group'))
# Now the dataframes are all ready to be completely customisable in ggplot
INR_G1_Rhizo_RH_betadisp_plot <- ggplot(INR_G1_Rhizo_RH_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme +
  scale_color_manual(name="G1 Treatment", values = c("#3399FF","#FFCC00"), labels = c("Control", 'Drought')) +
  ylab('Distance to Median') + xlab("G1 Treatment") +
  scale_x_discrete(labels = c("Control", 'Drought')) +
  theme(legend.position = "none")
INR_G1_Rhizo_RH_betadisp_plot


### MSU: Generation 2 ###

# Generation 2, test G1 treatment, G2 treatment, and G1_G2 groups

# 1. MSU G2 Root
MSU_G2_Root_betadisper_G1 <- betadisper(MSU_G2_Root_average_bc, group=as.factor(meta_MSU_G2_Root$G1_treatment), type = "median")
MSU_G2_Root_betadisper_G2 <- betadisper(MSU_G2_Root_average_bc, group=as.factor(meta_MSU_G2_Root$G2_treatment), type = "median")
MSU_G2_Root_betadisper_G1G2 <- betadisper(MSU_G2_Root_average_bc, group=as.factor(meta_MSU_G2_Root$G1_G2), type = "median")
set.seed(13)
permutest(MSU_G2_Root_betadisper_G1, pairwise=TRUE, permutations=999) #NS
set.seed(13)
permutest(MSU_G2_Root_betadisper_G2, pairwise=TRUE, permutations=999) #NS
set.seed(13)
permutest(MSU_G2_Root_betadisper_G1G2, pairwise=TRUE, permutations=999) #NS
TukeyHSD(MSU_G2_Root_betadisper_G1G2)
# Nothing is significant
boxplot(MSU_G2_Root_betadisper_G1G2)

# get betadisper data ####
MSU_G2_Root_BDdata <- get_betadisper_data(MSU_G2_Root_betadisper_G1G2)
MSU_G2_Root_BDdata$eigenvalue <- mutate(MSU_G2_Root_BDdata$eigenvalue, percent = eig/sum(eig))
MSU_G2_Root_BDdata$chull <- group_by(MSU_G2_Root_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# Now the dataframes are all ready to be completely customisable in ggplot
MSU_G2_Root_betadisp_plot <- ggplot(MSU_G2_Root_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme +
  scale_color_manual(name="G1_G2 Consecutive Treatment", 
                     values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), 
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) +
  ylab('Distance to Median') + xlab("G1_G2 Consecutive Treatment") +
  scale_x_discrete(labels = c("C_C", 'C_D', "D_C", "D_D")) +
  theme(legend.position = "none")
MSU_G2_Root_betadisp_plot

# 2. MSU G2 Rhizosphere
MSU_G2_Rhizo_betadisper_G1 <- betadisper(MSU_G2_Rhizo_average_bc, group=as.factor(meta_MSU_G2_Rhizo$G1_treatment), type = "median")
MSU_G2_Rhizo_betadisper_G2 <- betadisper(MSU_G2_Rhizo_average_bc, group=as.factor(meta_MSU_G2_Rhizo$G2_treatment), type = "median")
MSU_G2_Rhizo_betadisper_G1G2 <- betadisper(MSU_G2_Rhizo_average_bc, group=as.factor(meta_MSU_G2_Rhizo$G1_G2), type = "median")
set.seed(13)
permutest(MSU_G2_Rhizo_betadisper_G1, pairwise=TRUE, permutations=999) #NS
set.seed(13)
permutest(MSU_G2_Rhizo_betadisper_G2, pairwise=TRUE, permutations=999)
#            Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.0030534 0.00305342 5.1359    999  0.028 *
#Pairwise comparisons:
# Control vs. Drought           0.032 (permuted p-value)
set.seed(13)
permutest(MSU_G2_Rhizo_betadisper_G1G2, pairwise=TRUE, permutations=999) #NS
#Pairwise comparisons:
# C_D vs. D_C           0.023 (permuted p-value)
TukeyHSD(MSU_G2_Rhizo_betadisper_G1G2) #NS
boxplot(MSU_G2_Rhizo_betadisper_G1G2)

# get betadisper data ####
MSU_G2_Rhizo_BDdata <- get_betadisper_data(MSU_G2_Rhizo_betadisper_G1G2)
MSU_G2_Rhizo_BDdata$eigenvalue <- mutate(MSU_G2_Rhizo_BDdata$eigenvalue, percent = eig/sum(eig))
MSU_G2_Rhizo_BDdata$chull <- group_by(MSU_G2_Rhizo_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# Now the dataframes are all ready to be completely customisable in ggplot

mG2rhizo_BD_anno <- data.frame(x1 = c(2), x2 = c(3), 
                               y1 = c(0.48), y2 = c(0.485), 
                               xstar = c(2.5), ystar = c(0.487),
                               lab = c("*"))
mG2rhizo_BD_anno

MSU_G2_Rhizo_betadisp_plot <- ggplot(MSU_G2_Rhizo_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme + labs(title="(USA) G2 Rhizosphere - Red Hawk") +
  scale_color_manual(name="G1_G2 Consecutive\nTreatment", values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"),
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) +
  ylab('') + xlab("") +
  #scale_y_continuous(limits = c(0.3, 0.5))+
  scale_x_discrete(labels = c("C_C", 'C_D', "D_C", "D_D")) + 
  labs(caption="G2 treatment *") + 
  theme(plot.tag = element_text(face="bold"),
        plot.caption = element_text(size=13, face="italic", hjust=0)) 
MSU_G2_Rhizo_betadisp_plot
MSU_G2_Rhizo_betadisp_plot2 <- MSU_G2_Rhizo_betadisp_plot + 
  geom_text(inherit.aes=FALSE, data = mG2rhizo_BD_anno, aes(x = xstar,  y = ystar, label = lab), size=6) +
  geom_segment(inherit.aes=FALSE, data = mG2rhizo_BD_anno, aes(x = x1, xend = x1, 
                                                               y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes=FALSE, data = mG2rhizo_BD_anno, aes(x = x2, xend = x2, 
                                                               y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes=FALSE, data = mG2rhizo_BD_anno, aes(x = x1, xend = x2, 
                                                               y = y2, yend = y2),
               colour = "black")
MSU_G2_Rhizo_betadisp_plot2


### INRAE: Generation 2 ###

# 1. INRAE G2 Root Flavert
INR_G2_Root_Flav_betadisper_G1 <- betadisper(INR_G2_Root_Flav_average_bc, group=as.factor(meta_INR_G2_Root_Flav$G1_treatment), type = "median")
INR_G2_Root_Flav_betadisper_G2 <- betadisper(INR_G2_Root_Flav_average_bc, group=as.factor(meta_INR_G2_Root_Flav$G2_treatment), type = "median")
INR_G2_Root_Flav_betadisper_G1G2 <- betadisper(INR_G2_Root_Flav_average_bc, group=as.factor(meta_INR_G2_Root_Flav$G1_G2), type = "median")
set.seed(13)
permutest(INR_G2_Root_Flav_betadisper_G1, pairwise=TRUE, permutations=999) # NS
set.seed(13)
permutest(INR_G2_Root_Flav_betadisper_G2, pairwise=TRUE, permutations=999) # NS
set.seed(13)
permutest(INR_G2_Root_Flav_betadisper_G1G2, pairwise=TRUE, permutations=999) # NS
TukeyHSD(INR_G2_Root_Flav_betadisper_G1G2) #NS
boxplot(INR_G2_Root_Flav_betadisper_G1G2)
# get betadisper data ####
INR_G2_Root_Flav_BDdata <- get_betadisper_data(INR_G2_Root_Flav_betadisper_G1G2)
INR_G2_Root_Flav_BDdata$eigenvalue <- mutate(INR_G2_Root_Flav_BDdata$eigenvalue, percent = eig/sum(eig))
INR_G2_Root_Flav_BDdata$chull <- group_by(INR_G2_Root_Flav_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# the plot
INR_G2_Root_Flav_betadisper_plot <- ggplot(INR_G2_Root_Flav_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme + labs(title="(FR) G2 Roots - Flavert ") +
  scale_color_manual(name="G1_G2 Consecutive Treatment", values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"),
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) +
  ylab('Distance to Median') + xlab("") +
  scale_x_discrete(labels = c("C_C", 'C_D', "D_C", "D_D")) +
  theme(legend.position = "none", plot.tag = element_text(face="bold")) 
INR_G2_Root_Flav_betadisper_plot

# 2. INRAE G2 Root Red Hawk
INR_G2_Root_RH_betadisper_G1 <- betadisper(INR_G2_Root_RH_average_bc, group=as.factor(meta_INR_G2_Root_RH$G1_treatment), type = "median")
INR_G2_Root_RH_betadisper_G2 <- betadisper(INR_G2_Root_RH_average_bc, group=as.factor(meta_INR_G2_Root_RH$G2_treatment), type = "median")
INR_G2_Root_RH_betadisper_G1G2 <- betadisper(INR_G2_Root_RH_average_bc, group=as.factor(meta_INR_G2_Root_RH$G1_G2), type = "median")
set.seed(13)
permutest(INR_G2_Root_RH_betadisper_G1, pairwise=TRUE, permutations=999) # NS
set.seed(13)
permutest(INR_G2_Root_RH_betadisper_G2, pairwise=TRUE, permutations=999) # NS
set.seed(13)
permutest(INR_G2_Root_RH_betadisper_G1G2, pairwise=TRUE, permutations=999) # NS
#Pairwise comparisons:
# C_D vs. D_C 0.01* 
TukeyHSD(INR_G2_Root_RH_betadisper_G1G2) #NS
boxplot(INR_G2_Root_RH_betadisper_G1G2)
# get betadisper data ####
INR_G2_Root_RH_BDdata <- get_betadisper_data(INR_G2_Root_RH_betadisper_G1G2)
INR_G2_Root_RH_BDdata$eigenvalue <- mutate(INR_G2_Root_RH_BDdata$eigenvalue, percent = eig/sum(eig))
INR_G2_Root_RH_BDdata$chull <- group_by(INR_G2_Root_RH_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# Now the dataframes are all ready to be completely customisable in ggplot

iG2rootRH_BD_anno <- data.frame(x1 = c(2), x2 = c(3), 
                                y1 = c(0.26), y2 = c(0.28),
                                xstar = c(2.5), ystar = c(0.3),
                                lab = c("*"))

INR_G2_Root_RH_betadisper_plot <- ggplot(INR_G2_Root_RH_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme +
  scale_color_manual(name="G1_G2 Consecutive Treatment", values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"), 
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) +
  #ylab('Distance to Median') +
  labs(title="(FR) G2 Root - Red Hawk") +
  scale_x_discrete(labels = c("C_C", 'C_D', "D_C", "D_D")) +
  theme(legend.position = "none", axis.title = element_blank(), plot.tag = element_text(face="bold")) +
  geom_text(inherit.aes=FALSE, data = iG2rootRH_BD_anno, aes(x = xstar,  y = ystar, label = lab), size=6) +
  geom_segment(inherit.aes=FALSE, data = iG2rootRH_BD_anno, aes(x = x1, xend = x1, y = y1, yend = y2), 
               colour = "black") +
  geom_segment(inherit.aes=FALSE, data = iG2rootRH_BD_anno, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes=FALSE, data = iG2rootRH_BD_anno, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black")

INR_G2_Root_RH_betadisper_plot

# 3. INRAE G2 Rhizosphere Flavert
INR_G2_Rhizo_Flav_betadisper_G1 <- betadisper(INR_G2_Rhizo_Flav_average_bc, group=as.factor(meta_INR_G2_Rhizo_Flav$G1_treatment), type = "median")

INR_G2_Rhizo_Flav_betadisper_G2 <- betadisper(INR_G2_Rhizo_Flav_average_bc, group=as.factor(meta_INR_G2_Rhizo_Flav$G2_treatment), type = "median")

INR_G2_Rhizo_Flav_betadisper_G1G2 <- betadisper(INR_G2_Rhizo_Flav_average_bc, group=as.factor(meta_INR_G2_Rhizo_Flav$G1_G2), type = "median")

set.seed(13)
permutest(INR_G2_Rhizo_Flav_betadisper_G1, pairwise=TRUE, permutations=999) # NS
set.seed(13)
permutest(INR_G2_Rhizo_Flav_betadisper_G2, pairwise=TRUE, permutations=999) 
#           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.072719 0.072719 4.0549    999  0.048 *
#Pairwise comparisons:
# Control vs. Drought 0.044*
set.seed(13)
permutest(INR_G2_Rhizo_Flav_betadisper_G1G2, pairwise=TRUE, permutations=999)
#           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     3 0.15615 0.052051 4.1156    999   0.02 *
#Pairwise comparisons:
# C_C vs C_D *
# C_C vs D_C *
# C_C vs D_D **
TukeyHSD(INR_G2_Rhizo_Flav_betadisper_G1G2) 
# C_C vs C_D *
boxplot(INR_G2_Rhizo_Flav_betadisper_G1G2)

# get betadisper data ####
INR_G2_Rhizo_Flav_BDdata <- get_betadisper_data(INR_G2_Rhizo_Flav_betadisper_G1G2)
INR_G2_Rhizo_Flav_BDdata$eigenvalue <- mutate(INR_G2_Rhizo_Flav_BDdata$eigenvalue, percent = eig/sum(eig))
INR_G2_Rhizo_Flav_BDdata$chull <- group_by(INR_G2_Rhizo_Flav_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# Now the dataframes are all ready to be completely customisable in ggplot
iG2rhizoF_BD_anno <- data.frame(x1 = c(1, 1, 1), x2 = c(2, 3, 4), 
                                y1 = c(0.2, 0.15, 0.1), y2 = c(0.18, 0.13, 0.08), 
                                xstar = c(1.5, 2, 2.5), ystar = c(0.150, 0.100, 0.050),
                                lab = c("*", "*", "**"))
iG2rhizoF_BD_anno

INR_G2_Rhizo_Flav_betadisper_plot <- ggplot(INR_G2_Rhizo_Flav_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme + labs(title="(FR) G2 Rhizosphere - Flavert") +
  scale_color_manual(name="G1_G2 Consecutive Treatment", values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"),
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) +
  ylab('Distance to Median') +
  xlab("") +
  labs(caption="G1 x G2 treatment *") + 
  scale_y_continuous(limits = c(0, 0.8))+
  scale_x_discrete(labels = c("C_C", 'C_D', "D_C", "D_D")) +
  theme(legend.position = "none", plot.tag = element_text(face="bold"),
        plot.caption = element_text(size=13, face="italic", hjust=0))+ 
  geom_text(inherit.aes=FALSE, data = iG2rhizoF_BD_anno, aes(x = xstar,  y = ystar, label = lab), size=6) +
  geom_segment(inherit.aes=FALSE, data =iG2rhizoF_BD_anno, aes(x = x1, xend = x1, y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes=FALSE, data = iG2rhizoF_BD_anno, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes=FALSE, data = iG2rhizoF_BD_anno, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black")
INR_G2_Rhizo_Flav_betadisper_plot

# 4. INRAE G2 Rhizosphere Red Hawk
INR_G2_Rhizo_RH_betadisper_G1 <- betadisper(INR_G2_Rhizo_RH_average_bc, group=as.factor(meta_INR_G2_Rhizo_RH$G1_treatment), type = "median")
INR_G2_Rhizo_RH_betadisper_G2 <- betadisper(INR_G2_Rhizo_RH_average_bc, group=as.factor(meta_INR_G2_Rhizo_RH$G2_treatment), type = "median")
INR_G2_Rhizo_RH_betadisper_G1G2 <- betadisper(INR_G2_Rhizo_RH_average_bc, group=as.factor(meta_INR_G2_Rhizo_RH$G1_G2), type = "median")
set.seed(13)
permutest(INR_G2_Rhizo_RH_betadisper_G1, pairwise=TRUE, permutations=999) # NS
set.seed(13)
permutest(INR_G2_Rhizo_RH_betadisper_G2, pairwise=TRUE, permutations=999) # NS
set.seed(13)
permutest(INR_G2_Rhizo_RH_betadisper_G1G2, pairwise=TRUE, permutations=999) # NS
TukeyHSD(INR_G2_Rhizo_RH_betadisper_G1G2) #NS
boxplot(INR_G2_Rhizo_RH_betadisper_G1G2)
# get betadisper data ####
INR_G2_Rhizo_RH_BDdata <- get_betadisper_data(INR_G2_Rhizo_RH_betadisper_G1G2)
INR_G2_Rhizo_RH_BDdata$eigenvalue <- mutate(INR_G2_Rhizo_RH_BDdata$eigenvalue, percent = eig/sum(eig))
INR_G2_Rhizo_RH_BDdata$chull <- group_by(INR_G2_Rhizo_RH_BDdata$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()

# Now the dataframes are all ready to be completely customisable in ggplot
INR_G2_Rhizo_RH_betadisper_plot <- ggplot(INR_G2_Rhizo_RH_BDdata$distances, aes(group, distances, col = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.55)) +
  geom_point(size=3, shape=21, fill="white") +
  my_theme +
  scale_color_manual(name="G1_G2 Consecutive Treatment", values = c("#3399FF","#99CCFF", "#FFCC00", "#FF6600"),
                     labels = c("Control_Control", "Control_Drought", "Drought_Control", "Drought_Drought")) +
  ylab('Distance to Median') + xlab("G1_G2 Consecutive Treatment") +
  scale_x_discrete(labels = c("C_C", 'C_D', "D_C", "D_D")) +
  theme(legend.position = "none")
INR_G2_Rhizo_RH_betadisper_plot 


### Save the full beta dispersion plots ###

beta_dispersion_full <- INR_G1_Rhizo_Flav_betadisp_plot +
                        INR_G2_Root_RH_betadisper_plot + INR_G2_Rhizo_Flav_betadisper_plot + 
                        MSU_G2_Rhizo_betadisp_plot2 + plot_annotation(tag_levels = 'A')
beta_dispersion_full

ggsave(plot=beta_dispersion_full, "/Users/emiliedehon/Nextcloud/Documents/PROJECT/Bean_Drought_Multigenerations/Figures/NewFigures31_07_25/Fig.8.tif",
       width=10, height=8,
       device = "tiff",
       units= "in", dpi = 400,
       compression="lzw", bg= "white")






