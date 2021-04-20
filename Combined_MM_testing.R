# This script will combine the analysis and create final figures for my final combined data. 

getwd()
setwd("/home/coreyschultz/1.Projects/2.Maize.Endophyte.Project/Final_Data_and_Figs")
getwd()
library(tidyr)
library(tidyverse)
library(broom)
library(ggplot2)
library(car)
library(data.table)
library(ggpubr)

Herb_data = read.csv("Final_Herb_Data.csv", sep =",")
Burk_data = read.csv("Final_Burk_Data.csv", sep =",")
Serendip_data = read.csv("Final_Sbecsii_Data.csv", sep =",")

############## Do a variance decomposition loop with the mixed modeol no rep or residuals
HerbResultsT1 <- data.frame(Chlorophyll1=numeric(4),Chlorophyll2=numeric(4),Chlorophyll3=numeric(4),PlantHeight=numeric(4),LeafArea=numeric(4),RootLength=numeric(4),RootVolume=numeric(4))
rownames(HerbResultsT1) <- c("Group","Genotype","Inoculation","Genotype:Inoculation")

# Sum of Squares analysis loop

myvars <- names(Herb_data[5:11]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Herb_data[[m]] ~ (1 | Genotype) + (1 | Condition) + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Herb_data, REML = FALSE)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
}




SereResultsT1 <- data.frame(PlantHeight=numeric(4),RootLength=numeric(4),RootMass=numeric(4),ShootMass=numeric(4))
rownames(SereResultsT1) <- c("Group","Genotype","Inoculation","Genotype:Inoculation")

myvars <- names(Serendip_data[5:8]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  #linmod <- lmer(Serendip_data[[m]] ~ Group_or_Date + Genotype + Condition + Genotype:Condition + (1 | Group_or_Date/Genotype), data = Serendip_data, REML = FALSE)
  #linmod <- lmer(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition + (1 | Group_or_Date), data = Serendip_data, REML = TRUE)
  linmod <- lmer(Serendip_data[[m]] ~  (1 | Genotype) + (1 | Condition) + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  HPH1 <- anova(linmod, type = c("I"))
  print(HPH1) #anova
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
  #print(summary(linmod)) # then print mm summary for residuals
  SereResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2])
  itlist = HPH1[,6]
  Signif_list <- append(Signif_list, itlist)
}

# All random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~  (1 | Genotype) + (1 | Condition) + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
  print("/n")
  print(summary(linmod))
  print(deviance(linmod))
  print("################################################################")
}

# Only group random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  print(varcomp, comp = "Variance")
  print("/n")
  print(summary(linmod))
  print(deviance(linmod))
  print("################################################################")
  
}

# Anova just group random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition + Group_or_Date + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  HPH1 <- anova(linmod, type=c("I"))
  print(HPH1)
  print(summary(linmod))
}

# Anova plain lm
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Serendip_data[[m]] ~ (1 |Group_or_Date) + Genotype + Condition + Genotype:Condition, data = Serendip_data)
  HPH1 <- anova(linmod)
  print(HPH1)
}

# Interaction +  group random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ Genotype + Condition + (1 |Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
}

# Interaction + genotype + group random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ (1 |Genotype) + Condition + (1 |Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
}

# Just Genotype and group
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ (1 |Genotype) + Condition + Genotype:Condition + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
}


library(pander)

linmod <- lmer(PlantHeight ~  (1 | Genotype) + (1 | Condition) + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data)
HPH1 <- anova(linmod, type = c("I"))
pander(HPH1)
varcomp <- VarCorr(linmod)
print(varcomp, comp = "Variance")
get_variance(linmod)


linmod <- lmer(PlantHeight ~ Genotype + Condition + Genotype:Condition + Group_or_Date + (1 | Group_or_Date), data = Serendip_data)
HPH1 <- anova(linmod, type = c("I"))
pander(HPH1)

linmod <- lmer(PlantHeight ~ Genotype + Condition + Genotype:Condition + Group_or_Date + (1 | Group_or_Date), data = Serendip_data)
HPH1 <- anova(linmod, type = c("III"))
pander(HPH1)

print(aov(PlantHeight ~ Genotype + Condition + Genotype:Condition + Error(Group_or_Date), data = Serendip_data))

library(nlme)
m1 <- lme(PlantHeight ~ Genotype + Condition + Genotype:Condition, random=~1|Group_or_Date, data = Serendip_data)
anova(m1)
summary(m1)


# tommorow: show which models do or dont converge. Which residuals and variance you can get, and what you seem to not understande. Jason will help. 