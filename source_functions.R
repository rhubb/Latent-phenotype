##--------------------------------------------------------------------------##
## Functions and JAGS model specifications for Bayesian latent phenotyping 
## model
##
## Last updated: 12/4/17
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
## Helper Functions
##--------------------------------------------------------------------------##

expit <- function(x){
  exp(x)/(1 + exp(x))
}

# function for identifying patterns of missingness and converting 
# these data accordingly into lists for use in JAGS estimation
parse_miss<-function(dat){
  # Get names of columns that contain missingness
  miss_var_names<-colnames(dat)[ apply(dat,2,function(x) sum(is.na(x)))>0 ]
  miss_dat<-dat[,miss_var_names]
  dat$comb<-apply(miss_dat,1,function(x) paste0(is.na(x),collapse = '-') )
  
  map<-unique(dat$comb)
  map<-1:length(map)
  # Map is all the combinations of hba1c-rpg missingness
  names(map)<-c("FALSE-FALSE","TRUE-FALSE","FALSE-TRUE","TRUE-TRUE")
  dat$comb<-map[dat$comb]
  split_dat<-split(dat,dat$comb)
  dl<-lapply(split_dat, as.list)
  ll<-do.call('c',dl)
  names(ll)<-paste0(substring(text = names(ll),
                              first = as.numeric(gregexpr("\\.",names(ll)))+1, 
                              last=nchar(names(ll)) ),
                    substring(text = names(ll),
                              first = 1, 
                              last=as.numeric(gregexpr("\\.",names(ll)))-1 ) )
  ind<-unlist(lapply(ll,function(l) length(l)!=sum(is.na(l)) ))
  ll<-ll[ ind ]
  return(ll)
}

# counts of missingness in hba1c and rpg in a data.frame
count_miss<-function(d){
  miss_hba1c<-mean(is.na(d$hba1c))
  miss_rpg<-mean(is.na(d$rpg))
  return(c(mhba1c=miss_hba1c, mrpg=miss_rpg))
}


##--------------------------------------------------------------------------##
## JAGS model for estimation for PEDSnet data 
##--------------------------------------------------------------------------##

jagsmod_peds <- "
model{

  # Specify likelihood
  for(i in 1:sum(Npat)){
    
    # Clinical Codes and Medications #
    logit(pDM2[i]) <- codeDM2_b_int + codeDM2_b_dm*DM[i] # T2D
    logit(pEN[i]) <- codeEN_b_int + codeEN_b_dm*DM[i] # Endocrinologist visit
    logit(pMT[i]) <- codeMT_b_int + codeMT_b_dm*DM[i] # Metformin prescription
    logit(pIN[i]) <- codeIN_b_int + codeIN_b_dm*DM[i] # Insulin prescription
    
    T2D[i] ~ dbern(pDM2[i])
    endocrinologist[i] ~ dbern(pEN[i])
    metformin[i] ~ dbern(pMT[i])
    insulin[i] ~ dbern(pIN[i])
    
    # Missingness in biomarkers #
    logit(pMiss_hba1c[i])<- b0 + b2*age[i] + b3*bmi[i] + b_DM_miss*DM[i]
    logit(pMiss_rpg[i])<- a0 + a1*highrisk[i] + a2*age[i] + a3*bmi[i] +
                          a_DM_miss*DM[i]
    
    hba1c_miss_ind[i] ~ dbern(pMiss_hba1c[i])
    rpg_miss_ind[i] ~ dbern(pMiss_rpg[i])
    
    # Probability of latent T2DM #
    logit(rho[i])<-r0[i] + r1*highrisk[i] + r2*age[i] + r3*bmi[i]
  
  }
  
  # Loop over subjects with non-missing HbA1c # 
  for (i in 1:(Npat[1]+Npat[2])){
    hba1c[i] ~ dnorm(hba1c_b_int + hba1c_b_dm*DM[i], hba1c_tau)
  }
  
  # Loop over subjects with non-missing RPG
  for (i in 1:Npat[1]){
    rpg[i] ~ dnorm(rpg_b_int + rpg_b_dm*DM[i], rpg_tau)
  }
  for (i in (Npat[1]+Npat[2]+1):(Npat[1]+Npat[2]+Npat[3])){
    rpg[i] ~ dnorm(rpg_b_int + rpg_b_dm*DM[i], rpg_tau)
  }
  
  ## Specify priors
  hba1c_b_int ~ dnorm(6,0.1)
  hba1c_b_dm ~ dnorm(2.9,1) # informative prior based on AUC = 0.95
  hba1c_tau ~ dgamma(0.001,1000)
  hba1c_sigma <- 1/hba1c_tau    
  
  rpg_b_int ~ dnorm(85,0.1)
  rpg_b_dm ~ dnorm(75.6,0.1) # informative prior based on AUC = 0.95
  rpg_tau ~ dgamma(0.001,1000)
  rpg_sigma <- 1/rpg_tau
  
  codeDM2_b_int  ~ dunif(-7,-5)
  codeDM2_b_dm ~ dnorm(0,0.1)
  
  codeEN_b_int  ~ dunif(-7,-1)
  codeEN_b_dm ~ dnorm(0,0.1)
  
  codeMT_b_int  ~ dunif(-7,-1)
  codeMT_b_dm ~ dnorm(0,0.1)
  
  codeIN_b_int  ~ dunif(-7,-1)
  codeIN_b_dm ~ dnorm(0,0.1)
  
  a0 ~ dnorm(0,0.1)
  a1 ~ dnorm(0,0.1)
  a2 ~ dnorm(0,0.1)
  a3 ~ dnorm(0,0.1)
  a_DM_miss ~ dnorm(0,0.1)
  
  b0 ~ dnorm(0,0.1)
  b1 ~ dnorm(0,0.1)
  b2 ~ dnorm(0,0.1)
  b3 ~ dnorm(0,0.1)
  b_DM_miss ~ dnorm(0,0.1)
  
  r1 ~ dnorm(0,10)
  r2 ~ dnorm(0,10)
  r3 ~ dnorm(0,10)
  
  for (i in 1:sum(Npat)){
    DM[i] ~ dbern(rho[i])
    r0[i] ~ dunif(-7,-1)
  }
}
"


##--------------------------------------------------------------------------##
## Estimation Functions
##--------------------------------------------------------------------------##

# run_bayes_peds() is a wrapper function for run.jags()
# it takes as input the data set as well as arguments to be passed to run.jags 
# it returns posterior samples for model probability of T2DM and other model
# parameters specificed by Monitor
# function also returns input data set sorted by missingness pattern
run_bayes_peds<-function(pedsdata, jagsmod, Burnin, sample, method, adapt, Monitor){

    # Store original data set
    temp <- pedsdata
    
    # Select only columns of data set needed by JAGS
    pedsdata <- pedsdata[, c("person_id","gender","endocrinologist",
                             "metformin","insulin","T2D")]

    # Center continuous variables
    pedsdata$age <- temp$age_at_baseline-mean(temp$age_at_baseline) # want the continuous variables to be centered
    pedsdata$bmi <- temp$first.bmi_zscore-mean(temp$first.bmi_zscore)
    
    # Create numeric binary indicator variables
    pedsdata$highrisk <- as.numeric(temp$high.risk.race.ethnicity)
    pedsdata$highrisk <- ifelse(is.na(pedsdata$highrisk)==TRUE, 0, pedsdata$highrisk)
    pedsdata$endocrinologist <- as.numeric(pedsdata$endocrinologist)
    pedsdata$metformin <- as.numeric(pedsdata$metformin)
    pedsdata$insulin <- as.numeric(pedsdata$insulin)
    pedsdata$T2D <- as.numeric(pedsdata$T2D)
    
    # Rename predictors
    pedsdata$hba1c <- temp$avg.hemoglobin.A1c_percentage
    pedsdata$rpg <- temp$avg.random.glucose_mg.dl
    
    # Create missingness indicators
    pedsdata$hba1c_miss_ind <- ifelse(is.na(temp$avg.hemoglobin.A1c_percentage)==TRUE, 1, 0)
    pedsdata$rpg_miss_ind <- ifelse(is.na(temp$avg.random.glucose_mg.dl)==TRUE, 1, 0)
    
    # Add vector with number of occurences of each missingness pattern to data
    misspat <- paste(pedsdata$hba1c_miss_ind,pedsdata$rpg_miss_ind)
    pedsdata <- pedsdata[order(misspat),]
    outdata <- data.frame(pedsdata,T1D = temp[order(misspat),"T1D"],eMERGE = temp[order(misspat),"eMERGE"])
  
    dat_list <- list()
    dat_list$`Npat` <- as.vector(table(misspat))
    dat_list <- c(dat_list, as.list(pedsdata))
    
    # Run JAGS
    results<-run.jags(jagsmod_peds, 
                      monitor=Monitor, 
                      data = dat_list,
                      burnin = Burnin, sample = sample, n.chains = 1,
                      method = method, adapt = adapt, silent.jags = TRUE) 

    
    # Format MCMC samples
    results_mcmc<-as.mcmc.list(results)

    output <- list(mcmc_list=results_mcmc, outdata = outdata) 
    
    return(output)
}