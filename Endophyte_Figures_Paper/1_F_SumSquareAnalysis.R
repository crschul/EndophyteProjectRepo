# This script will combine the analysis and create final figures for my final combined data. 

getwd()                      # PLEASE CHANGE WD
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
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Serendip_data)
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

# Average Chlorophyll together
Herb_data$Chlorophyll <- rowMeans(Herb_data[ , c("Chlorophyll1","Chlorophyll2","Chlorophyll3")], na.rm=TRUE)

HerbResultsT1 <- data.frame(Chlorophyll=numeric(n),PlantHeight=numeric(n),LeafArea=numeric(n),RootLength=numeric(n),RootVolume=numeric(n))
rownames(HerbResultsT1) <- c("Genotype","Inoculation","Genotype:Inoculation", "Residuals")

herb_cols <- c(13, 8:11)
myvars <- names(Herb_data[,herb_cols]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Genotype + Condition + Genotype:Condition, data = Herb_data)
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
  gather(HerbT1, value,Chlorophyll:RootVolume)

HerbT1_long$Signif <- Signif_list


# Burkholderia
Burk_data$RootLength <- as.numeric(Burk_data$RootLength)
Burk_data$RootVolume <- as.numeric(Burk_data$RootVolume)

burkResultsT1 <- data.frame(PlantHeight=numeric(n),LeafArea=numeric(n),RootLength=numeric(n),RootVolume=numeric(n))
rownames(burkResultsT1) <- c("Genotype","Inoculation","Genotype:Inoculation", "Residuals")

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
SereT1_long["Endophyte"] <- "Experiment 3 \n (Serendipita)"
HerbT1_long["Endophyte"] <- "Experiment 1 \n (Herbaspirillum)"
BurkT1_long["Endophyte"] <- "Experiment 2 \n (Burkholderia)"


names(SereT1_long)[names(SereT1_long) == "ST1"] <- "Phenotype"
names(HerbT1_long)[names(HerbT1_long) == "HerbT1"] <- "Phenotype"
names(BurkT1_long)[names(BurkT1_long) == "BurkT1"] <- "Phenotype"

SandHandB <- rbind(HerbT1_long,BurkT1_long,SereT1_long)

cbPalette <- c("#999999", "forestgreen", "tan3", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Add spaces to phenotypes

SandHandB$Phenotype[SandHandB$Phenotype == "PlantHeight"] <- "Plant Height"
SandHandB$Phenotype[SandHandB$Phenotype == "LeafArea"] <- "Leaf Area"
SandHandB$Phenotype[SandHandB$Phenotype == "RootLength"] <- "Root Length"
SandHandB$Phenotype[SandHandB$Phenotype == "RootVolume"] <- "Root Volume"
SandHandB$Phenotype[SandHandB$Phenotype == "ShootMass"] <- "Shoot Mass"
SandHandB$Phenotype[SandHandB$Phenotype == "RootMass"] <- "Root Mass"

# Change names in legends to fit the manuscript text
SandHandB$rn[SandHandB$rn == "Genotype"] <- "Maize Genotype"
SandHandB$rn[SandHandB$rn == "Inoculation"] <- "Endophyte Inoculation"
SandHandB$rn[SandHandB$rn == "Genotype:Inoculation"] <- "Genome-by-Genome Interaction"

# reorder legend 
SandHandB$rn <- factor(SandHandB$rn, levels = c("Maize Genotype","Genome-by-Genome Interaction","Endophyte Inoculation","Residuals"))

# why do I need to use another vjust? cus its turned sideways? that is stupid
fig_1 <- ggplot(SandHandB, aes(x = Phenotype, y = value, fill = forcats::fct_rev(rn), label = Signif)) + 
  geom_col(position=position_stack()) + theme(axis.text.x = element_text( size = 14)) + 
  theme(axis.text.y = element_text(size = 14)) + 
  labs(fill = "Variables") + geom_text(aes(label = Signif), size = 10, position = position_stack(vjust = 0.5), 
                                       vjust = .75)+ 
  xlab("Measured Phenotypes") + coord_flip() +
  ylab("Proportion of SS") + scale_fill_manual(values=cbPalette) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white"), 
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) + 
  facet_grid(rows = vars(Endophyte), drop = TRUE, space = "free", scales = "free", switch = NULL) + 
  theme(legend.position = "right") + theme(axis.line = element_line(colour = "white"), 
                                           panel.border = element_blank()) + 
  theme(legend.title = element_text(size=18)) + theme(legend.text = element_text(size=14)) + 
  theme(strip.text.y = element_text(size = 14)) + theme(plot.title = element_text(size=20)) +
  theme(axis.title.x = element_text(30)) + theme(axis.title.y = element_text(20)) + 
  guides(fill = guide_legend(reverse = TRUE)) + theme(legend.position = "bottom")


ggsave("Fig1_Endo_Submission", fig_1, device = "png", width = 10, height = 5, dpi = 800)

# Need to change Genotype: Inoculation name - just changed column names above
SandHandB <- as.data.frame(lapply(SandHandB, function(x) {
  gsub("Genotype:Inoculation", "GenxInoc", x)
}))

# Calculate total variance due to genotype 
SandHandB

# Need to change Genotype: Inoculation name - just changed column names above
SandHandB <- as.data.frame(lapply(SandHandB, function(x) {
   gsub("Genotype:Inoculation", "GenxInoc", x)
 }))

count = 0.0
sum = 0.0

for (row in 1:nrow(SandHandB)){
  if(grepl("Genome-by-Genome Interaction", SandHandB[row, 1])){
    count = count + 1
    sum = sum + as.numeric(as.character(SandHandB[row,3]))
    #print(SandHandB[row,3])
  }}

interaction_avg = (sum/count) 
print(interaction_avg)

count = 0.0
sum = 0.0

for (row in 1:nrow(SandHandB)){
  if(grepl("Genotype", SandHandB[row, 1])){
    count = count + 1
    sum = sum + as.numeric(as.character(SandHandB[row,3]))
    #print(SandHandB[row,3])
  }}

Genotype_avg = (sum/count) 
print(Genotype_avg)

print("Interaction over Genotype")
print(interaction_avg / Genotype_avg)

