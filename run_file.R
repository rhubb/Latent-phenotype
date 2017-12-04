##--------------------------------------------------------------------------##
## Apply Bayesian latent phenotyping model to pediatric T2DM
##
## Last updated: 12/4/17
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
## Load Packages, Functions and Data
##--------------------------------------------------------------------------##
library(runjags)
library(rjags)

source("source_functions.R") # source functions and JAGS models

##--------------------------------------------------------------------------##
## Run Bayesian Model
##--------------------------------------------------------------------------##

# Set parameters to monitor
monitor <- c("rho","hba1c_b_int","hba1c_b_dm","hba1c_sigma","rpg_b_int",
             "rpg_b_dm","rpg_sigma","codeDM2_b_int","codeDM2_b_dm",
             "codeEN_b_int","codeEN_b_dm","codeMT_b_int","codeMT_b_dm",
             "codeIN_b_int","codeIN_b_dm","a0","a1","a2","a3","a_DM_miss","b0",
             "b1","b2","b3","b_DM_miss","r1","r2","r3")

# Run latent phenotyping model defined by jagsmod_peds using data in pedsdata
peds_bayes<-run_bayes_peds(pedsdata, jagsmod =jagsmod_peds,Monitor=monitor, Burnin=1000, 
                      sample=5000, method='simple', adapt=1000)

# Save MCMC object and data set sorted by missingness pattern
outdata <- peds_bayes[[2]]
peds_bayes <- peds_bayes[[1]]

# Extract posterior probabilities of DM
rho <- peds_bayes[[1]][,which(regexpr("rho",colnames(peds_bayes[[1]]))>0)]

# Extract other model parameters
parm <- peds_bayes[[1]][,which(regexpr("rho",colnames(peds_bayes[[1]]))<=0)]


##--------------------------------------------------------------------------##
# Summarize Posterior Means and CIs
##--------------------------------------------------------------------------##

# Transform parameters for missing HbA1c or RPG to OR
parm[,c("a_DM_miss","b_DM_miss")] <- exp(parm[,c("a_DM_miss","b_DM_miss")])

postmeans <- apply(parm,2,mean)
postci    <-  apply(parm,2,quantile, probs = c(0.025, 0.975))

codesens <- apply(parm[,c("codeDM2_b_int", "codeEN_b_int", "codeMT_b_int", 
                  "codeIN_b_int")] + parm[,c("codeDM2_b_dm", "codeEN_b_dm",
                  "codeMT_b_dm", "codeIN_b_dm")],2,function(x){mean(expit(x))})
codesens.ci <- apply(parm[,c("codeDM2_b_int", "codeEN_b_int", "codeMT_b_int", 
                  "codeIN_b_int")] + parm[,c("codeDM2_b_dm", "codeEN_b_dm", 
                  "codeMT_b_dm", "codeIN_b_dm")],2,function(x){quantile(expit(x), 
                  probs = c(0.025, 0.975))})
codespec <- apply(parm[,c("codeDM2_b_int", "codeEN_b_int", "codeMT_b_int", 
                  "codeIN_b_int")],2,function(x){mean(1-expit(x))})
codespec.ci <- apply(parm[,c("codeDM2_b_int", "codeEN_b_int", "codeMT_b_int", 
                  "codeIN_b_int")],2,function(x){quantile(1-expit(x), 
                  probs = c(0.025, 0.975))})

# Compute posterior mean probability of T2DM
postDM <- apply(rho,2,mean)


