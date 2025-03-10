
#Package installation (some packages that cant be installed with install.pacakges())
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("graph")
#BiocManager::install("RBGL")
#BiocManager::install("Rgraphviz")

#Packages
library(pcalg)
library(micd) #For mixed data tests
library(Rgraphviz) #Graph visualization
library(readr)
library(factoextra)
library(tibble)
library(dplyr)

#Directory
setwd("XXX")

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

#PCA preparation
{
  #PCA Nutrients
  nutrients_components = prcomp(nutrients_data[,-1], scale. = TRUE)
  fviz_eig(nutrients_components)
  #nutrient_values = nutrients_components$rotation
  pca_nutrients <- as.data.frame(nutrients_components$x)
  nutrients_pca_data = data.frame(ID = nutrients_data$ID, 
                                  PC_nutri_1 = pca_nutrients$PC1,
                                  PC_nutri_2 =pca_nutrients$PC2,
                                  PC_nutri_3 =pca_nutrients$PC3)
  
  #PCA Food groups (sin HEI y SCORE GABAS)
  fg_components = prcomp(foodGroups_data[,2:12], scale. = TRUE)
  fviz_eig(fg_components)
  #fg_rotations = fg_components$rotation
  pca_fg = as.data.frame(fg_components$x)
  fg_pca_data = data.frame(ID = foodGroups_data$ID,
                           PC_fg_1 = pca_fg$PC1,
                           PC_fg_2 = pca_fg$PC2,
                           PC_fg_3 = pca_fg$PC3)
}

#Numeric confounder histograms -> for checking normality
{
  hist(anthro_data$age) #Not normal-looking
  
  hist(anthro_data$bmi) #Gaussian-like
  
  hist(anthro_data$body_fat) #Gaussian-like
  
  hist(anthro_data$systolic_bp) #Gaussian-like
  
  hist(anthro_data$diastolic_bp) #Gaussian-like
  
  hist(anthro_data$triglycerides) #Not normal-looking
  
  hist(anthro_data$HDL) #Gaussian-like
  
  hist(anthro_data$glucose) #Not normal looking
  
  hist(anthro_data$HOMA_IR) #Highly affected by outliers
  anthro_data_no_out = anthro_data[anthro_data$HOMA_IR < 10,]
  hist(anthro_data_no_out$HOMA_IR) #Gaussian-like (not really)
  
  hist(anthro_data$hsCRP) #Highly affected by outliers
  antrho_data_no_hscrp_out = anthro_data[anthro_data$hsCRP < 3,]
  hist(antrho_data_no_hscrp_out$hsCRP) #Not normal-looking
  
  hist(nutrients_pca_data$PC_nutri_1) #Gaussian-like
  hist(nutrients_pca_data$PC_nutri_2) #Gaussian-like
  hist(nutrients_pca_data$PC_nutri_3) #Seems affected by outliers
  nutrients_pca_data_no_out = nutrients_pca_data[nutrients_pca_data$PC_nutri_3<5, ]
  hist(nutrients_pca_data_no_out$PC_nutri_3) #Gaussian-like
  
  hist(fg_pca_data$PC_fg_1) #Gaussian-like (barely, not so much)
  hist(fg_pca_data$PC_fg_2) #Seems affected by otuliers
  fg_pca_no_out = fg_pca_data[fg_pca_data$PC_fg_2 > (0-5),]
  hist(fg_pca_no_out$PC_fg_2) #Gaussian-like
  hist(fg_pca_data$PC_fg_3) #Gaussiasn-like
}

#Normalization of non-gaussian-looking confounders
{
  anthro_data$age = scale(anthro_data$age)
  anthro_data$triglycerides = scale(anthro_data$triglycerides)
  anthro_data$glucose = scale(anthro_data$glucose)
  anthro_data$hsCRP = scale(anthro_data$hsCRP)
}

#Creation of dummy variables
{
  sex_yes = ifelse(anthro_data$sex=="Male", 1, 0)
  city_BGA = ifelse(anthro_data$city=="Bucaramanga", 1, 0)
  city_MED = ifelse(anthro_data$city=="Medellin", 1, 0)
  city_BAR = ifelse(anthro_data$city=="Barranquilla", 1, 0)
  city_CLI = ifelse(anthro_data$city=="Cali", 1, 0)
}

#Creation of confounding variables matrices: Anthro confounders, food-related confounders
{
  anthro_confounders <- data.frame(
    ID = anthro_data$ID,
    sex_yes = sex_yes,
    city_BAR = city_BAR,
    city_MED = city_MED,
    city_BGA = city_BGA,
    city_CLI = city_CLI,
    age = anthro_data$age,
    bmi = anthro_data$bmi,
    body_fat = anthro_data$body_fat,
    systolic_bp = anthro_data$systolic_bp,
    diastolic_bp = anthro_data$diastolic_bp,
    triglycerides = anthro_data$triglycerides,
    HDL = anthro_data$HDL,
    glucose = anthro_data$glucose,
    HOMA_IR = anthro_data$HOMA_IR,
    hsCRP = anthro_data$hsCRP
  )
  
  food_related_confounders = inner_join(nutrients_pca_data, fg_pca_data, by = "ID")
}

#Identification of low observation OTUs
{
  #OTU observation count
  total_otu_observation_counts = colSums(otu_data[,-1])
  
  #OTUs with at most 1 observation
  zero_observation_otus = which(total_otu_observation_counts < 2)
  length(zero_observation_otus) #There are 18 OTUs with no observations -> we dont need to use them
  zero_observation_otus
  
  #Add 1 for column indices (to account for ID)
  zero_observation_otus = zero_observation_otus + 1
  
  #Select only OTUs with observations
  observed_otus= otu_data[,-zero_observation_otus]
  
  #Sanity checks
  observed_OTU_test = colSums(observed_otus[,-1])
  length(which(observed_OTU_test == 0)) == 0
  
  observed_OTU_test = colSums(observed_otus[,-1])
  length(which(observed_OTU_test == 1)) == 0
  
}

#OTU transformations -> Relative abundances 
{
  #Relative abundance
  otu_no_ids = observed_otus[,-1]
  
  abundances = otu_no_ids / rowSums(otu_no_ids) * 100
  otu_relative_abundances = cbind(observed_otus[,1], abundances)
}

#Creation of complete data matrix
{
  #Join all data
  complete_data_confounders = inner_join(anthro_confounders, food_related_confounders, by ="ID")
  complete_data = inner_join(complete_data_confounders, otu_relative_abundances, by = "ID")
  
  #Since we're working with Systolic BP, we remove MI_354_H because it has NA
  complete_data = complete_data[!(complete_data$ID %in% c("MI_354_H")),]
  
  #Remove ID Column
  complete_data = as.data.frame(complete_data[,-1])
  
}

#Causal Model Testing: Gradually increase p to see increase in computation time

#Causal model 1: 30 variables -> 9 OTUs (0.136s)
{
  #Matrix
  mc1 = complete_data[,c(1:30)]
  
  #Causal discovery and time to compute
  start_time = proc.time()
  
  mc1.fit = pc(suffStat = mc1, indepTest = mixCItest, alpha = 0.01,
               labels = colnames(mc1), skel.method = "stable.fast",
               numCores = 5)
  
  end_time = proc.time()
  print("mc1 fit time:")
  (end_time - start_time)
  
  #Estimation of Causal effect of arbitrary OTU over BMI -> local reduce computation
  mc1.causalEffect = ida(29, 7, cov(mc1), mc1.fit@graph, method = "local")
  mc1.causalEffect
  
  #Verification if graph is a valid cpdag
  mc1.adjMatrix = as(mc1.fit@graph, "matrix")
  isValidGraph(mc1.adjMatrix, type ="cpdag", verbose = TRUE)
  
  #Plot
  plot(mc1.fit@graph, main ="MC1 Fit")
}

#Causal model 2: 150 variables (2.6s)
{
  #Matrix
  mc2 = complete_data[,c(1:150)]
  
  #Causal discovery and time to compute
  start_time = proc.time()
  
  mc2.fit = pc(suffStat = mc2, indepTest = mixCItest, alpha = 0.01,
               labels = colnames(mc2), skel.method = "stable.fast",
               numCores = 5)
  
  end_time = proc.time()
  print("mc2 fit time:")
  (end_time - start_time)
  
  #Estimation of Causal effect of arbitrary OTU over BMI -> local reduce computation
  mc2.causalEffect = ida(115, 7, cov(mc2), mc2.fit@graph, method = "local")
  mc2.causalEffect
  
  #Verification if graph is a valid cpdag
  mc2.adjMatrix = as(mc2.fit@graph, "matrix")
  isValidGraph(mc2.adjMatrix, type ="cpdag", verbose = TRUE)
  
}

#Causal model 3: 500 variables (46.276s)
{
  #Matrix
  mc3 = complete_data[,c(1:500)]
  
  #Causal discovery and time to compute
  start_time = proc.time()
  
  mc3.fit = pc(suffStat = mc3, indepTest = mixCItest, alpha = 0.01,
               labels = colnames(mc3), skel.method = "stable.fast",
               numCores = 5)
  
  end_time = proc.time()
  print("mc3 fit time:")
  (end_time - start_time)
  
  #Estimation of Causal effect of arbitrary OTU over BMI -> local reduce computation
  mc3.causalEffect = ida(359, 7, cov(mc3), mc3.fit@graph, method = "local")
  mc3.causalEffect
  
  #Verification if graph is a valid cpdag
  mc3.adjMatrix = as(mc3.fit@graph, "matrix")
  isValidGraph(mc3.adjMatrix, type ="cpdag", verbose = TRUE)
  
}

#Causal model 4: 1000 variables (205.193s)
{
  #Matrix
  mc4 = complete_data[,c(1:1000)]
  
  #Causal discovery and time to compute
  start_time = proc.time()
  
  mc4.fit = pc(suffStat = mc4, indepTest = mixCItest, alpha = 0.01,
               labels = colnames(mc4), skel.method = "stable.fast",
               numCores = 5)
  
  end_time = proc.time()
  print("mc4 fit time:")
  (end_time - start_time)
  
  #Estimation of Causal effect of arbitrary OTU over BMI -> local reduce computation
  mc4.causalEffect = ida(678, 7, cov(mc4), mc4.fit@graph, method = "local")
  mc4.causalEffect
  
  #Verification if graph is a valid cpdag
  mc4.adjMatrix = as(mc4.fit@graph, "matrix")
  isValidGraph(mc4.adjMatrix, type ="cpdag", verbose = TRUE)
  
}

#Causal model 5: All variables
{
  #Matrix
  mc5 = complete_data
  
  #Causal discovery and time to compute
  start_time = proc.time()
  
  mc4.fit = pc(suffStat = mc5, indepTest = mixCItest, alpha = 0.01,
               labels = colnames(mc5), skel.method = "stable.fast",
               numCores = 5)
  
  end_time = proc.time()
  print("mc5 fit time:")
  (end_time - start_time)
  
  #Estimation of Causal effect of arbitrary OTU over BMI -> local reduce computation
  mc5.causalEffect = ida(678, 7, cov(mc5), mc5.fit@graph, method = "local")
  mc5.causalEffect
  
  #Verification if graph is a valid cpdag
  mc5.adjMatrix = as(mc5.fit@graph, "matrix")
  isValidGraph(mc5.adjMatrix, type ="cpdag", verbose = TRUE)
  
}


