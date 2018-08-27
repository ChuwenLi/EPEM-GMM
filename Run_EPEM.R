library(tidyr)
require(gtools)
library(matrixcalc)
library(mvtnorm)
library(MASS)
library(ggplot2)
library(gridExtra) 
library(Matrix)
library(magic)
library(nlme)
library(lme4)
library(data.table)
library(mvnfast) 
library(emulator)
library(sas7bdat)
library(tibble)


##################################################################################
## Source program that holds all functions used to run the program
source("C:/Users/dlu7/Documents/thesiswork/GMM/dlu7-thesis/SimulationStudyFinal/PenEM/programs/GitHubFiles/EM_GMM_PenLik.R")

#### These need to be set by the user!!
epsilon= 1E-4  ## Convergence criteria for FE parameters
B      = 10000 ## Number of max iterations for EM algorithm until convergence
p      = 2     ## FE polynomial order 
q      = 2     ## RE polynomial order

dat.use<-read.csv(file="SimData_R4_N400_K4_p2_q2.csv",header=T)
head(dat.use)

timepoints = seq(1,10,1)
Ruse       = 4

## Data matrix needs to have the following column names
## RepNum
## gene
## V1.1 to V1.10 (where the unique time-points are after the period)

## Running EPEM using all available data (does not average over replicates)
test<-EPEM_function(dat.use=dat.use, ## Input dataset
              timepoints=timepoints, ## vector of unique time-points used in the dataset
              Ruse=4,                ## Nubmer of replicates
              iterview=TRUE)         ## If you don't want to see all the results at each iteration set to FALSE

## Running EPEM by averaging over replicates 
test2<-EPEM_function(dat.use=dat.use,      
                    timepoints=timepoints, 
                    Ruse=1,          ## If you want to average over replicates, set R=1          
                    iterview=FALSE)

