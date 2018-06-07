##Required functions and packages
expit = function(x){
  exp(x)/(1+exp(x))
}

require(MCMCpack)

###Design function 
PS.design = function(treatment, #treatment vector
                     covardataset, #n x p covariate matrix where p is the number of covariates
                     bayes = TRUE, #draw from PS posterior distribution or estimate MLE
                     K=1000){ #number of draws from the PS posterior distribution
  
  dataset = data.frame(treatment,covardataset)
  n = length(treatment) #sample size
  
  ###Bayesian PS estimation
  if(bayes){
    #draw from posterior distribution of alpha
    posteriors = as.matrix(MCMClogit(treatment~.,data=dataset,burnin=5*K,mcmc=K*50,thin=50))
    
    covars = as.matrix(cbind(rep(1,n),covardataset))
    ps = expit(covars%*%t(posteriors)) #PS estimated from draws of alpha
  } else {
    
    #########Frequentist PS estimation
    #estimate PS using logistic regression
    ps = predict(glm(treatment~.,data=dataset,family=binomial),type="response")
  
  }
  
    return(ps) #returns a n x K matrix of covariates, K = 1 for frequentist PSA
}