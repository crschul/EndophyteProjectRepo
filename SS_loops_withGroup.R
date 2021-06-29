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


#get sd for error bars
data_summary <- function(data, varname, groupnames){
  require(plyr)
  
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



# Fig 4 - Sum Square analysis
# With loops comment out below

# number of variables in our model + residuals
n = 4

# Serendipita
SereResultsT1 <- data.frame(PlantHeight=numeric(n),RootLength=numeric(n),RootMass=numeric(n),ShootMass=numeric(n))
rownames(SereResultsT1) <- c("Group","Inoculation","Geno:Inoculation", "Residuals")

myvars <- names(Serendip_data[5:8]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Serendip_data[[m]] ~ Group_or_Date + Condition + Genotype:Condition, data = Serendip_data)
  HPH1 <- anova(linmod)
  print(HPH1)
  SereResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2])
  itlist = HPH1[,5]
  Signif_list <- append(Signif_list, itlist)
}

#Get rid of the NAs from residuals
Signif_list[is.na(Signif_list)] = 1
Signif_list

# Convert these p values to *
for( s in 1:length(Signif_list)){
  #print(s)
  if (Signif_list[s] < .05){
    Signif_list[s] <- '*'
    #print("True")
  } else {
    Signif_list[s] <- " "
    #print("False")
  }
}
Signif_list

# Now Normalize all the columns
ST1 <- SereResultsT1

ST1[] <- lapply(ST1[], function(x) x/sum(x))
setDT(ST1, keep.rownames = TRUE)[]

SereT1_long <- ST1 %>%
  gather(ST1, value,PlantHeight:ShootMass)

SereT1_long$Signif <- Signif_list


# Herbaspirillum
HerbResultsT1 <- data.frame(Chlorophyll1=numeric(n),Chlorophyll2=numeric(n),Chlorophyll3=numeric(n),PlantHeight=numeric(n),LeafArea=numeric(n),RootLength=numeric(n),RootVolume=numeric(n))
rownames(HerbResultsT1) <- c("Group","Inoculation","Geno:Inoculation", "Residuals")

myvars <- names(Herb_data[5:11]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Group_or_Date + Condition + Genotype:Condition, data = Herb_data)
  HPH1 <- anova(linmod)
  print(HPH1)
  HerbResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2])
  itlist = HPH1[,5]
  Signif_list <- append(Signif_list, itlist)
}

#Get rid of the NAs from residuals
Signif_list[is.na(Signif_list)] = 1
Signif_list

# Convert these p values to *
for( s in 1:length(Signif_list)){
  #print(s)
  if (Signif_list[s] < .05){
    Signif_list[s] <- '*'
    #print("True")
  } else {
    Signif_list[s] <- " "
    #print("False")
  }
}
Signif_list

# Now Normalize all the columns
HerbT1 <- HerbResultsT1

HerbT1[] <- lapply(HerbT1[], function(x) x/sum(x))
setDT(HerbT1, keep.rownames = TRUE)[]

HerbT1_long <- HerbT1 %>%
  gather(HerbT1, value,Chlorophyll1:RootVolume)

HerbT1_long$Signif <- Signif_list


# Burkholderia
Burk_data$RootLength <- as.numeric(Burk_data$RootLength)
Burk_data$RootVolume <- as.numeric(Burk_data$RootVolume)

burkResultsT1 <- data.frame(PlantHeight=numeric(n),LeafArea=numeric(n),RootLength=numeric(n),RootVolume=numeric(n))
rownames(burkResultsT1) <- c("Genotype","Inoculation","Geno:Inoculation", "Residuals")

myvars <- names(Burk_data[5:8]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Burk_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Burk_data)
  HPH1 <- anova(linmod)
  print(HPH1)
  burkResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2])
  itlist = HPH1[,5]
  Signif_list <- append(Signif_list, itlist)
}

#Get rid of the NAs from residuals
Signif_list[is.na(Signif_list)] = 1
Signif_list

# Convert these p values to *
for( s in 1:length(Signif_list)){
  #print(s)
  if (Signif_list[s] < .05){
    Signif_list[s] <- '*'
    #print("True")
  } else {
    Signif_list[s] <- " "
    #print("False")
  }
}
Signif_list

# Now Normalize all the columns
BurkT1 <- burkResultsT1

BurkT1[] <- lapply(BurkT1[], function(x) x/sum(x))
setDT(BurkT1, keep.rownames = TRUE)[]

BurkT1_long <- BurkT1 %>%
  gather(BurkT1, value,PlantHeight:RootVolume)

BurkT1_long$Signif <- Signif_list


# Endophytes included in Name: don't like that but without empty faceting
SereT1_long["Endophyte"] <- "Serendipita"
HerbT1_long["Endophyte"] <- "Herbaspirillum"
BurkT1_long["Endophyte"] <- "Burkholderia"


names(SereT1_long)[names(SereT1_long) == "ST1"] <- "Phenotype"
names(HerbT1_long)[names(HerbT1_long) == "HerbT1"] <- "Phenotype"
names(BurkT1_long)[names(BurkT1_long) == "BurkT1"] <- "Phenotype"

SandHandB <- rbind(HerbT1_long,BurkT1_long,SereT1_long)
SandH <-rbind(HerbT1_long,SereT1_long)

cbPalette <- c("#999999", "forestgreen", "tan3", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fig3_cap = "Figure 4: Sum square analysis of phenotypes for all genotypes, separated by endophyte. Asterisks denote significance based on ANOVA analysis. This Sum Square analysis allows us to visualize where the variance in phenotype is orginating from."

ggplot(SandH, aes(x = Phenotype, y = value, fill = forcats::fct_rev(rn), label = Signif)) + 
  geom_col(position=position_stack()) + theme(axis.text.x = element_text( size = 12)) + 
  theme(axis.text.y = element_text(size = 14)) + 
  labs(fill = "Variables") + geom_text(aes(label = Signif), size = 10, position = position_stack(vjust = .5), hjust = .25)+ 
  ggtitle("Breakdown of Sum Squares Analysis") + xlab("Measured Phenotypes") + coord_flip() +
  ylab("Proportion of SS") + scale_fill_manual(values=cbPalette) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white"), 
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) + 
  facet_grid(rows = vars(Endophyte), drop = TRUE, space = "free", scales = "free") + 
  theme(legend.position = "right") + theme(axis.line = element_line(colour = "white"), 
                                           panel.border = element_blank()) +
  labs(caption=str_wrap(fig3_cap, 100)) + 
  theme(plot.caption = element_text(hjust = 0.5,size = 16))

# Calculate total variance due to genotype 
SandHandB
# Need to change Genotype: Inoculation name - just changed column names above
# SandHandB <- as.data.frame(lapply(SandHandB, function(x) {
#   gsub("Genotype:Inoculation", "GenxInoc", x)
# }))

count = 0.0
sum = 0.0

for (row in 1:nrow(SandH)){
  if(grepl("Group", SandH[row, 1])){
    count = count+1
    sum = sum + (SandH[row,3])
    #print(SandHandB[row,3])
  }}
print(sum/count) 
# Without Group Genotype Variance = .55
# With Group Genotype Variance = .39
# Group Variance = .188
# Group without Genotype = 