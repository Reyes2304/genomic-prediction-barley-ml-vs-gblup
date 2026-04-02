
###-----------------------------------------------------------------------
#      phenotypic data analysis for heritability and BLUEs
###-----------------------------------------------------------------------
## set working directory
setwd(".:/....")

# libraries
library(asreml)
library(data.table)

# load the input raw data
spring <- fread("input_raw_data_spring.txt")

spring$Env <- factor(spring$Env)
spring$replication_org<- factor(spring$replication_org)
spring$block<- factor(spring$block)
spring$dummy_rep<-factor(spring$dummy_rep)
spring$dummy_unrep<-factor(spring$dummy_unrep)
spring$genotypes<-factor(spring$genotypes)

# function for BLUEs and heritability calculation
## Slight adjustment has to be done for the model to analysis other traits if the model cannot converge
## a license for asreml-R is needed

BLUE <- function(data, trait) {
  
  # Convert trait column to numeric
  data[[trait]] <- as.numeric(data[[trait]])
  
  # Calculate heritability components
  rec <- length(data$genotypes[!is.na(data[[trait]])])       			         	    
  Nr.genos <- length(unique(data$genotypes[!is.na(data[[trait]])]))            
  Nr.rep.ac <- rec / Nr.genos                                                                 
  Nr.env <- length(unique(data$Env[!is.na(data[[trait]])])) 
  
  # Fit the model for heritability calculation
  asr <- asreml(fixed = as.formula(paste(trait, "~ 1")), 
                random = ~ Env + genotypes + dummy_rep:Env:replication_org + dummy_rep:Env:replication_org:block + genotypes:Env,
                na.action = na.method(x = "omit"),
                data = data, maxit = 100)
  
  varcomp <- summary(asr)$varcomp
  
  # Extract variance components
  sigma.Gen <- varcomp["genotypes", "component"]
  sigma.Env <- varcomp["Env", "component"]
  sigma.Gen.x.Env <- varcomp["genotypes:Env", "component"] 
  sigma.e <- varcomp["units!R", "component"]
  
  # Calculate heritability
  h2.Gen <- sigma.Gen / (sigma.Gen + sigma.Gen.x.Env / Nr.env + sigma.e / Nr.rep.ac)
  
  Heritability <- data.frame("Trait" = trait,
                             "Sigma.Gen" = sigma.Gen,
                             "Sigma.Gen.x.Env" = sigma.Gen.x.Env,
                             "Sigma.Env" = sigma.Env,
                             "Sigma.e" = sigma.e,
                             "h2.Gen" = h2.Gen,
                             "Nr.Rep" = Nr.rep.ac,
                             "Nr.Env" = Nr.env)
  
  # Fit model for BLUEs calculation
  asrl <- asreml(fixed = as.formula(paste(trait, "~ genotypes")), 
                 random = ~ Env + dummy_rep:Env:replication_org + dummy_rep:Env:replication_org:block + genotypes:Env,
                 na.action = na.method(x = "omit"),
                 data = data, maxit = 100)
  
  predict.fix <- predict(asrl,
                         classify = "genotypes",
                         asreml.options(maxit = 50, workspace = 64e6, pworkspace = 64e6),
                         data = data)
  
  Blues <- predict.fix$pvals[, c("genotypes", "predicted.value")]
  
  # output the result
  write.csv(Blues, file = paste0("varcomp_", trait, "_spring.csv"), row.names = FALSE)
  write.csv(Blues, file = paste0("BLUEs_", trait, "_spring.csv"), row.names = FALSE)
  
}

#run the model
BLUE(spring, "puccinia_hordei")
