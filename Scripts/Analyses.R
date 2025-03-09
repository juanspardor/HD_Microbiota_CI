
#Package installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("graph")
BiocManager::install("RBGL")
BiocManager::install("Rgraphviz")

#Packages
library(pcalg)
library(micd) #For mixed data tests
library(Rgraphviz) #Graph visualization
library(readr)
library(factoextra)
library(tibble)
library(dplyr)

#Directory
setwd("/Users/juanse/Documents/Tesis/IIND/HD_Microbiota_CI")

#Data Reading
{
  #Anthropometric
  anthro_data = read_delim("./Data/Consolidated/anthro_data.meta",
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
  
  #OTU observations
  otu_data = read_delim("./Data/Consolidated/otu_data.otus",
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
  
  #Nutrient data
  nutrients_data = read_delim("./Data/Consolidated/nutrients_data.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE)
  unwanted_columns_nut = c("Codalt")
  nutrients_data = nutrients_data[, !names(nutrients_data) %in% unwanted_columns_nut]
  
  #Food group data 
  foodGroups_data = read.csv("./Data/Consolidated/food_groups_u24h.csv")
  #Remove certain columns
  unwanted_columns_fg = c("Code", "Codalt", "City", "r24h")
  foodGroups_data = foodGroups_data[, !names(foodGroups_data) %in% unwanted_columns_fg]
}

#Numeric confounder histograms -> for checking normality
{
  hist(anthro_data$age) #Not normal-looking
  
  hist(anthro_data$bmi) #Gaussian-like
  
  hist(anthro_data$systolic_bp) #Gaussian-like
  
  hist(anthro_data$diastolic_bp) #Gaussian-like
  
  hist(anthro_data$triglycerides) #Not normal-looking
  
  hist(anthro_data$HDL) #Gaussian-like
  
  hist(anthro_data$glucose) #Not normal looking
  
  hist(anthro_data$HOMA_IR) #Highly affected by outliers
  anthro_data_no_out = anthro_data[anthro_data$HOMA_IR < 10,]
  hist(anthro_data_no_out$HOMA_IR) #Gaussian-like (not really)
  
  hist(anthro_data$hsCRP) #Highly affected by outliers
  antrho_data_no_hscrp_out = anthro_data[anthro_data$hsCRP < 4,]
  hist(antrho_data_no_hscrp_out$hsCRP) #Not normal-looking
}

#Normalization of non-gaussian-looking confounders
{
  
}

#Creation of dummy variables
{
  sex_yes = ifelse(anthro_data$sex=="Male", 1, 0)
  city_BGA = ifelse(anthro_data$city=="Bucaramanga", 1, 0)
  city_MED = ifelse(anthro_data$city=="Medellin", 1, 0)
  city_BAR = ifelse(anthro_data$city=="Barranquilla", 1, 0)
  city_CLI = ifelse(anthro_data$city=="Cali", 1, 0)
}

#Identification of 0 observation OTUs
{
  #OTU observation count
  total_otu_observation_counts = colSums(otu_data[,-1])
  
  #OTUs with no observations
  zero_observation_otus = which(total_otu_observation_counts == 0)
  length(zero_observation_otus) #There are 18 OTUs with no observations -> we dont need to use them
  zero_observation_otus
  
  #Add 1 for column indices (to account for ID)
  zero_observation_otus = zero_observation_otus + 1
  
  #Select only OTUs with observations
  observed_otus= otu_data[,-zero_observation_otus]
  
  #Sanity check
  observed_OTU_test = colSums(observed_otus[,-1])
  length(which(observed_OTU_test == 0) == 0) == 0

}

#OTU transformations -> Relative abundances + log transform
{
  
}