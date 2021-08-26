# This script will combine the analysis and create final figures for my final combined data. 

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

##################################################################################################
# Real Time - qPCR    
###################################################################################################
# Note: each sample had 2 technical reps! Experiment 3 had a different dna concentration due to differences in 
# extracted DNA concnetrations. (they were too low compared to exps 1 and 2)
library(Rmisc)
qpcr_data <- read.csv("CombinedRaw_qPCR_clean.csv", header = TRUE, sep = ",")
# Toss any rows with missing values
qpcr_data = qpcr_data[complete.cases(qpcr_data),]

#add empty column
qpcr_data['calcs'] <- NA

# delta ct: experimental gene - housekeeping gene - this is for a single sample
qpcr_delta1 <- transform.data.frame(qpcr_data, calcs = Cp_ITS - Cp_CDK)

# average all the delta ct values by genotype and condition so you can do the second delta from the combined samples
qpcr_sum <- summarySE(data = qpcr_delta1, measurevar = "calcs", groupvars = c("Genotype","Condition"), na.rm = TRUE)

qgenotypes <- unique(qpcr_sum$Genotype)
length(qgenotypes)

# find the double delta for each genotype
doubledelta_df <- data.frame("Genotype" = qgenotypes, "Double_Delta" = NA)

# Inoculated - Control 
for( i in 1:length(qgenotypes)){
  df <-qpcr_sum[which (qpcr_sum$Genotype == (unique(qpcr_sum$Genotype)[i])),c(1:7)]
  doubled <- df[2,4]-df[1,4]
  doubledelta_df[i,2] = doubled
}

#Add experiments so we can fill by them
doubledelta_df['Grow_Number'] <- c('3','3','3','1','3','2','1','2','1','1','1','3')

# Drop CML69 and MO18W - we have no phenotype data due to greenhouse fertilization by a worker
doubledelta_df <- doubledelta_df[!(doubledelta_df$Genotype=="CML69"),]
doubledelta_df <- doubledelta_df[!(doubledelta_df$Genotype=="MO18W"),]


# Transform to fold change
foldchange_df = transform.data.frame(doubledelta_df, Double_Delta = log(2^-Double_Delta))
#foldchange_df = transform.data.frame(doubledelta_df, Double_Delta = 2^(-Double_Delta))

# Add phenotype significance from t - tests .1 - .05 : This significance values are calculated for the actual phenotype in 2_T_GrowthSignif_TTest and applied here
foldchange_df['Significance'] <- c('~',' ',' ',' ','*','~',' ','*','*','*')


#Heritability
qpcrFinModel <- lm(qpcr_delta1$calcs~ qpcr_delta1$Genotype*qpcr_delta1$Condition, data=qpcr_delta1)
qPCRFintab <- anova(qpcrFinModel)
colonizationFinal.het <- qPCRFintab[3,2]/(qPCRFintab[3,2] + qPCRFintab[4,2])


foldchange_df$Grow_Number <- gsub('1', 'Grow 1', foldchange_df$Grow_Number)
foldchange_df$Grow_Number <- gsub('3', 'Grow 2', foldchange_df$Grow_Number)


# With Faceting:
fig4_qpcr <- ggplot(foldchange_df, aes(x = reorder(Genotype, -Double_Delta), Double_Delta, fill = Grow_Number, label = Significance)) + 
  geom_bar(stat = "summary", fun = "mean", position="dodge") + 
  theme(plot.title = element_text(size=18),axis.text.x = element_text(size=18, angle = 90), axis.text.y = element_text(size=18), axis.title = element_text(size = 18)) + 
  scale_fill_manual(values = c("#56B4E9", "tan3")) + 
  ggtitle("Differences in Serendipita Abundance by Genotype") + xlab("Genotype") + 
  ylab("Log Fold Change Fungal Abundance") + geom_text(aes(label = Significance), size = 10, position = position_stack(vjust = .5)) + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white"), 
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "light grey"))  +
  facet_grid(cols = vars(Grow_Number), drop = TRUE, space = "free", scales = "free") + #ylim(0,6) +
  theme(legend.position = "right") + theme(axis.line = element_line(colour = "white"), 
                                           panel.border = element_blank())  +
  theme(strip.text.x = element_text(size = 18)) + theme(legend.title = element_blank()) + theme(legend.position = "none")
# labs(caption=str_wrap(fig4_cap, 100)) + 
# theme(plot.caption = element_text(hjust = 0.5,size = 16))
# annotate("text", x = 3, y = 5, size=5, label = " P-Values for a phenotype less than:\n ~ = 0.1,    * = 0.05 ") +


dat_txt <- data.frame(
  label = c("  * p < 0.05  \n  ~ p < 0.1 "),
  Grow_Number = c("Grow 2")
)

fig4_qpcr + geom_text(data = dat_txt,
                      mapping = aes(x = 4.5, y = 9.0, label = label),
                      size = 8) + coord_cartesian(xlim = c(1,5), clip = 'off')
