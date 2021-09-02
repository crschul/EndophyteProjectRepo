# This script will test for variance assumptions

# shapiro - wilk normality test - less than .05 than normality is rejected - the least informative test of the bunch
# qq plots and residual histograms are used to see how far the data deviates from norm quantiles
# Levens test - less sensitive to normality - less than .05 means data variance is not equal

getwd()       # PLEASE CHANGE YOUR WD 
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

names(Herb_data) <- c("Genotype", "Condition", "Rep", "Endophyte","Chlorophyll1","Chlorophyll2","Chlorophyll3",
                      "PlantHeight","LeafArea","RootLength","RootVolume","Group_or_Date")
names(Burk_data) <- c("Genotype", "Condition", "Rep", "Endophyte",
                      "PlantHeight","LeafArea","RootLength","RootVolume")
names(Serendip_data) <- c("Genotype", "Condition", "Rep", "Endophyte",
                          "PlantHeight","RootLength","RootMass","ShootMass")

Serendip_data$Condition <- gsub('C', 'Control', Serendip_data$Condition)
Serendip_data$Condition <- gsub('I', 'Inoculated', Serendip_data$Condition)
Herb_data$Condition <- gsub('C', 'Control', Herb_data$Condition)
Herb_data$Condition <- gsub('I', 'Inoculated', Herb_data$Condition)
Burk_data$Condition <- gsub('C', 'Control', Burk_data$Condition)
Burk_data$Condition <- gsub('I', 'Inoculated', Burk_data$Condition)


# number of variables in our model + residuals
n = 4

# Serendipita #Drop rep
SereResultsT1 <- data.frame()
SereResultsT1 <- data.frame(PlantHeight=numeric(n),RootLength=numeric(n),RootMass=numeric(n),ShootMass=numeric(n))
rownames(SereResultsT1) <- c("Genotype","Inoculation","Genotype:Inoculation", "Residuals")

myvars <- names(Serendip_data[5:8]) # create a list of traits

for( m in myvars){
  print(m)
  print(as.name(m))
  qqPlot(Serendip_data[[m]])
  print(shapiro.test(Serendip_data[[m]]))
  print(leveneTest(Serendip_data[[m]] ~ Genotype*Condition, data = Serendip_data))
  linmod <- lm(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Serendip_data)
  linmod2 <- lm(Serendip_data[[m]] ~ Condition + Genotype + Condition:Genotype, data = Serendip_data)
  hist(residuals(linmod))
  print(shapiro.test(residuals(linmod)))
  anv <- anova(linmod)
  print(anv) # type 1
  print(anova(linmod2)) # type 1 switched
  print(Anova(lm(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Serendip_data),type = 2)) # Type 2
  print(Anova(lm(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Serendip_data),type = 3)) # type 3
  print(" ")
  print(" ******************************************************************************************")
}

# Herbaspirillum

# Average Chlorophyll together
Herb_data$Chlorophyll <- rowMeans(Herb_data[ , c("Chlorophyll1","Chlorophyll2","Chlorophyll3")], na.rm=TRUE)

HerbResultsT1 <- data.frame(Chlorophyll=numeric(n),PlantHeight=numeric(n),LeafArea=numeric(n),RootLength=numeric(n),RootVolume=numeric(n))
rownames(HerbResultsT1) <- c("Genotype","Inoculation","Genotype:Inoculation", "Residuals")

herb_cols <- c(13, 8:11)
myvars <- names(Herb_data[,herb_cols]) # create a list of traits

for( m in myvars){
  print(m)
  print(as.name(m))
  qqPlot(Herb_data[[m]])
  print(shapiro.test(Herb_data[[m]]))
  print(leveneTest(Herb_data[[m]] ~ Genotype*Condition, data = Herb_data))
  linmod <- lm(Herb_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Herb_data)
  linmod2 <- lm(Herb_data[[m]] ~ Condition + Genotype + Condition:Genotype, data = Herb_data)
  hist(residuals(linmod))
  print(shapiro.test(residuals(linmod)))
  anv <- anova(linmod)
  print(anv) # type 1
  print(anova(linmod2)) # type 1 switched
  print(Anova(lm(Herb_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Herb_data),type = 2)) # Type 2
  print(Anova(lm(Herb_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Herb_data),type = 3)) # type 3
  print(" ")
  print(" ******************************************************************************************")
  
}

# Burkholderia
Burk_data$RootLength <- as.numeric(Burk_data$RootLength)
Burk_data$RootVolume <- as.numeric(Burk_data$RootVolume)

burkResultsT1 <- data.frame(PlantHeight=numeric(n),LeafArea=numeric(n),RootLength=numeric(n),RootVolume=numeric(n))
rownames(burkResultsT1) <- c("Genotype","Inoculation","Genotype:Inoculation", "Residuals")

myvars <- names(Burk_data[5:8]) # create a list of traits

for( m in myvars){
  print(m)
  print(as.name(m))
  qqPlot(Burk_data[[m]] )
  print(shapiro.test(Burk_data[[m]] ))
  print(leveneTest(Burk_data[[m]]  ~ Genotype*Condition, data = Burk_data))
  linmod <- lm(Burk_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Burk_data)
  linmod2 <- lm(Burk_data[[m]] ~ Condition + Genotype + Condition:Genotype, data = Burk_data)
  hist(residuals(linmod))
  print(shapiro.test(residuals(linmod)))
  anv <- anova(linmod)
  print(anv) # type 1
  print(anova(linmod2)) # type 1 switched
  print(Anova(lm(Burk_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Burk_data),type = 2)) # Type 2
  print(Anova(lm(Burk_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Burk_data),type = 3)) # type 3
  print(" ")
  print(" **************************************************************")
}

