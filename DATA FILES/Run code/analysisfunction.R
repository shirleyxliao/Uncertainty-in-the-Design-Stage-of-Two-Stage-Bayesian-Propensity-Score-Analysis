###Analysis code

#required packages 
require(MCMCpack)
require(MatchIt)
require(survey) 

#necessary side functions
expit = function(x){
  exp(x)/(1+exp(x))
}

logit = function(x){
  log(x/(1-x))
}

mu1fxn = function(treatment,outcome,psvector,y1){
  ((treatment*outcome/psvector) -
     ((treatment-psvector)*y1/psvector))
  
}

mu0fxn = function(treatment,outcome,psvector,y0){
  
  ((1-treatment)*outcome/(1-psvector) +
     (treatment-psvector)*y0/(1-psvector))
}

drsdfxn = function(psvector,dr_mu0,dr_mu1,ey1,ey0,y1,y0,outcome,treatment){
  n = length(psvector)
  ipw_y1 = (treatment*outcome/psvector)
  ipw_y0 = ((1-treatment)*outcome/(1-psvector))
  ipw_mu1 = sum(ipw_y1)/sum(treatment/psvector)
  ipw_mu0 = sum(ipw_y0)/sum((1-treatment)/(1-psvector))
  
  sigma_ipw2 = sum((treatment*(outcome-ipw_mu1)/psvector - (1-treatment)*(outcome-ipw_mu0)/(1-psvector))^2)/n
  
  #sigma_ipw2 = mean((y1 - dr_mu1)^2/psvector +(y0 - dr_mu0)^2/(1-psvector))
  
  sqrt( sigma_ipw2 -
          (mean(sqrt((1-psvector)/psvector)*(y1 - dr_mu1) +
                  sqrt(psvector/(1-psvector))*(y0-dr_mu0)))^2)/sqrt(length(psvector))
}


#fxn for DR 
thetafxn = function(treatment,outcome,psvector,y1,y0){
  ((treatment*outcome/psvector) -
     ((treatment-psvector)*y1/psvector)) -
    ((1-treatment)*outcome/(1-psvector) +
       (treatment-psvector)*y0/(1-psvector))}


covar.bal = function(treated,control,covar.names,weights1,weights0){
  tryCatch({mean((apply((treated[,covar.names]*weights1),2,mean,na.rm=TRUE) -
          apply((control[,covar.names]*weights0),2,mean,na.rm=TRUE))/apply((treated[,covar.names]*weights1),2,sd,na.rm=TRUE))} , error=function(e){NA})
}

#analysis function
PS.analysis = function(psvector, ##psvector = input vector of ps 
                       treatment,  ##treatment = input indicator vector of treated/control
                       outcome,   ##outcome = input vector of outcomes
                       covardataset,  ##covardataset = matrix of covariates, ncol = # of covariates, nrows = # of observations
                       method,    ##method = takes "nn.matching", "caliper.matching", "ipw", "dr", and "stratification
                       matching.caliper = 0.5, 
                       strat.cut = 5,   ##strat.cut = number of strata for stratification analysis
                       bpsa = 1, ##bpsa-1 uses asymptotic methods (+ sandwich estimator), bpsa-2 samples from the conditional posterior distribution of delta
                       match.iter = 1, #how many times to perform matching conditional on a single PS
                       matching.ratio = 1, #1-how many matching?
                       sandwich.est = FALSE #sandwich estimators for weighted analyses?
                       ){ 
  #data structure creation
  dataset = data.frame(outcome,treatment,psvector,covardataset)
  
  #variable definition
  S = 200 #draws from the posterior distribution of delta
  K = 1000 #draws from the posterior distribution of PS
  ncovar = dim(covardataset)[2] #number of covariates
  n = length(psvector) #number of observations
  
  #fill in dummy answers so function doesn't break
  percent.data.used = 1
  estimated.effect = NA
  estimated.se = NA
  nu = NA
  n.treat = NA
  n.control=NA
  strata.prop = rep(NA,strat.cut)
  strata.TE = rep(NA,strat.cut)
  cov.balance = NA
  sandwich.se = NA
  
###############################################################
   #Caliper MATCHING with REPLACEMENT

if(method=="caliper.matching"){ #caliper matching with replacement
  #caliper matching, 1-1 with replacement, caliper = 0.5*ps sd, distance used is logit(ps)  
        
    #bpsa 1: uses sandwich estimator for weights + one nu per PS + single estimate of condtional mean ATE and SE
    if(bpsa==1){
    
    #perform random matching within caliper
    matchobj= suppressWarnings(matchit(treatment~psvector, data = dataset, method="nearest",replace=TRUE, 
                                         caliper=matching.caliper, 
                                          ratio = matching.ratio,
                                         distance=logit(dataset$psvector)))
      
      #matrix of matched set
    matchdata = match.data(matchobj)
    
    #number of observations matched
    nmatched = n - sum(matchobj$weights==0)
    
    #obtain ATE + sandwitch SE estimator
    des = svydesign(~0,data=matchdata,weights=matchdata$weights)
    matchfit = svyglm(outcome~treatment,design=des)
      
    #load in estimated ATE and SE
    estimated.effect = coef(matchfit)[2]
    estimated.se = summary(matchfit)$coef[2,2]
    percent.data.used = nmatched/n
    nu = matrix(0,nrow=n,ncol=1)
    rownames(nu) = 1:n
    nu[rownames(matchdata),1] = matchdata$weights
    n.control = matchobj$nn[2,1]
    n.treat = matchobj$nn[2,2]
    
    cov.balance = covar.bal(matchdata[matchdata$treatment==1,],
                            matchdata[matchdata$treatment==0,],
                            names(matchdata)[!names(matchdata)%in%c("outcome","treatment","psvector","distance","weights")],
                            matchdata$weights[matchdata$treatment==1],
                            matchdata$weights[matchdata$treatment==0])

     #bpsa 2: uses unweighted SE estimator + allows multiple nu calculated on 
    ## a single set of PS + outputs ATE as draws form the marginal distrib.
 } else{
   
   #create data objects which account for multiple implementations performed
   ##on a single PS
   nu = matrix(0,nrow=n,ncol=match.iter)
   rownames(nu) = 1:n
   perc.used = rep(NA,match.iter)
   effect = rep(NA,match.iter)
   sd.effect = rep(NA,match.iter)
   cond.n.match = matrix(NA,nrow=match.iter,ncol=2)
   all.cov.balance = rep(NA,match.iter)
   sandwich.ses = rep(NA,match.iter)
   
   #perform multiple implementations conditional on a single PS
   for(r in 1:match.iter){
     
     #use matching algorithm
     matchobj= suppressWarnings(matchit(treatment~psvector, data = dataset, method="nearest",replace=TRUE, 
                                        caliper=matching.caliper, 
                                        ratio = matching.ratio,
                                        distance=logit(dataset$psvector)))
     
     #matrix of matched set
     matchdata = match.data(matchobj)
     
     #save nu
     nu[rownames(matchdata),r] = matchdata$weights
     
     #number of matches
     nmatched = n - sum(matchobj$weights==0)
     perc.used[r] = nmatched/n
     
     #use sandwich estimator to calculate conditional var?

     des = svydesign(~0,data=matchdata,weights=matchdata$weights)
     matchfit = svyglm(outcome~treatment,design=des)
     sandwich.ses[r] = summary(matchfit)$coefficients[2,2]
     
    matchfit = lm(outcome~treatment,data=matchdata,weights=matchdata$weights)
     
     
     #load in estimated ATE and SE
     effect[r] = coef(matchfit)[2]
     sd.effect[r] = summary(matchfit)$coefficients[2,2]
     all.cov.balance[r] = covar.bal(matchdata[matchdata$treatment==1,],
                                    matchdata[matchdata$treatment==0,],
                                    names(matchdata)[!names(matchdata)%in%c("outcome","treatment","psvector","distance","weights")],
                                    matchdata$weights[matchdata$treatment==1],
                                    matchdata$weights[matchdata$treatment==0])
     cond.n.match[r,] =  c(matchobj$nn[2,1],matchobj$nn[2,2])
   }
     
   percent.data.used = mean(perc.used)
   
   #ONLY for caliper matching: ATE is outputted as draws from the marginal
   ##posterior distirbution of Delta
   estimated.effect = matrix(rnorm(S*match.iter,effect,sd.effect),ncol=match.iter)
   
   if(sandwich.est){
     estimated.se = mean(sandwich.ses) } else{
   estimated.se = mean(sd.effect)}
   
   sandwich.se = mean(sandwich.ses)
   
   #cond.ev.nu = apply(nu,1,mean)
   n.treat = mean(cond.n.match[,2])
   n.control = mean(cond.n.match[,1])
   cov.balance = mean(all.cov.balance)
 }
  
###########################################################################
  #NN MATCHING WITH REPLACEMENT
        
} else if(method=="nn.matching"){

  #get matches
  rownames(dataset)<-1:dim(dataset)[1]
  
  matchobj= suppressWarnings(matchit(treatment~psvector, data = dataset, method="nearest",replace=TRUE, 
                    caliper=matching.caliper, 
                    ratio = matching.ratio,
                    mahvars = c("psvector"),
                    distance=as.numeric(logit(dataset$psvector))))
  
  #matrix of match data
  matchdata = match.data(matchobj)
  
  ##number of observations in the matched set
  nmatched = n - sum(matchobj$weights==0)

  #save nu (frequency weights)
  nu = matrix(0,nrow=n,ncol=1)
  rownames(nu) = 1:n
 nu[rownames(matchdata),1] = matchdata$weights
  
 #fit weighted linear regression to obtain ATE and (sandwich optional) SE

  des = svydesign(~0,data=matchdata,weights=matchdata$weights)
  matchfit = svyglm(outcome~treatment,design=des)
  sandwich.se = summary(matchfit)$coefficients[2,2]
  
  if(!sandwich.est){
    matchfit = lm(outcome~treatment,data=matchdata,weights=matchdata$weights)
  }

  #save estimated ATE and SE
 estimated.effect = coef(matchfit)[2]
 estimated.se = summary(matchfit)$coefficients[2,2]
    #estimated.effect = rnorm(S,coef(matchfit)[2],estimated.se)
    
    percent.data.used= nmatched/n

    n.control = matchobj$nn[2,1]
    n.treat = matchobj$nn[2,2]
    cov.balance = covar.bal(matchdata[matchdata$treatment==1,],
              matchdata[matchdata$treatment==0,],
              names(matchdata)[!names(matchdata)%in%c("outcome","treatment","psvector","distance","weights")],
              matchdata$weights[matchdata$treatment==1],
              matchdata$weights[matchdata$treatment==0])
    
  
##################################################################################
  #DOUBLY ROBUST ESTIMATION
} else if(method =="dr"){
  
  regdataset = data.frame(outcome,treatment,covardataset)
  treatment.form = formula(paste("treatment ~ ",paste(names(data.frame(covardataset)),
                                                      collapse=" + "),sep=""))
  outcome.form = formula(paste("outcome ~ ",paste(names(data.frame(covardataset)),
                                                  collapse=" + "),sep=""))
  
  #calculate inverse probability weights
  weights = treatment/psvector + (1-treatment)/(1-psvector) 
  
  #fit potential outcome regression model, impute counterfactual outcomes
  
  #our DR estimate
  fit1 = suppressWarnings(lm(outcome.form,data=regdataset[treatment==1,]))
  fit0 = suppressWarnings(lm(outcome.form,data=regdataset[treatment==0,]))
  
  #imputed outcomes
  y0 = predict(fit0,newdata=regdataset)
  y1 = predict(fit1,newdata=regdataset)
  
  #theta = thetafxn(treatment,outcome,psvector,y0,y1)
  mu1 = mu1fxn(treatment,outcome,psvector,y1)
  mu0 = mu0fxn(treatment,outcome,psvector,y0)
  #theta = thetafxn(treatment,outcome,psvector,y0,y1)
  estimated.effect = mean(mu1)-mean(mu0)
  estimated.se = drsdfxn(psvector,dr_mu0=mean(mu0),dr_mu1=mean(mu1),ey1=mu1,ey0=mu0,y1,y0,outcome,treatment)
  
  nu = weights
  cov.balance = covar.bal(dataset[treatment==1,],
                          dataset[treatment==0,],
                          names(dataset)[!names(dataset)%in%c("outcome","treatment","psvector")],
                          weights[treatment==1],
                          weights[treatment==0])
  
  
  
  #fit = suppressWarnings(lm(outcome~.,data=regdataset))
    #atecovars1 =  regdataset #counterfactual "all treatment" world
    #atecovars1$treatment = rep(1,n)

    #atecovars0 =  regdataset #counterfactual "all control" world
    #atecovars0$treatment = rep(0,n)
    
    #imputed outcomes
    #y0 = predict(fit,newdata=atecovars0)
    #y1 = predict(fit,newdata=atecovars1)
  
    #theta = thetafxn(treatment,outcome,psvector,y0,y1)
   # estimated.effect = mean(theta)
    #estimated.se = sqrt(sum((theta-estimated.effect)^2)/n^2)
    #nu = weights
    #cov.balance = covar.bal(dataset[treatment==1,],
     #         dataset[treatment==0,],
     #         names(dataset)[!names(dataset)%in%c("outcome","treatment","psvector")],
      #        weights[treatment==1],
      #        weights[treatment==0])
    
#######################################################################
    #IPW
}else if(method=="ipw"){
  
  #calculate inverse probability weights
  weights = treatment/psvector + (1-treatment)/(1-psvector) 
  
  #fit weighted regression with or without sandwich estimators
  
  #Independent Sampling design (with replacement)
  des = svydesign(~0,data=dataset,weights=weights)
  fit = svyglm(outcome~treatment,design=des) #sandwich estimator 
  sandwich.se = summary(fit)$coef[2,2]
  
  if(!sandwich.est){ 
    fit = lm(outcome~treatment,data=dataset,weights=weights)
    }

  
  #save estimated ATE and SE
  estimated.effect = summary(fit)$coef[2,1]
  estimated.se = summary(fit)$coef[2,2]
  nu = weights
  cov.balance = covar.bal(dataset[treatment==1,],
                          dataset[treatment==0,],
                          names(dataset)[!names(dataset)%in%c("outcome","treatment","psvector")],
                          weights[treatment==1],
                          weights[treatment==0])
  
  
 ###############################################################  
  ##Stratification, quintiles (redefine every time)
  
  } else if(method == "stratification"){
    
    #separate observations into quartile strata
    dataset$quartile = NA
    quartilebounds = quantile(dataset$psvector,prob=seq(0,1,by=1/strat.cut))
    goodstrata = 1:strat.cut
    n1 = rep(NA,strat.cut)
    n0 = rep(NA,strat.cut)
    all.cov.balance = rep(NA,strat.cut)
    
   for(k in 1:strat.cut){
    dataset$quartile[quartilebounds[k]<=dataset$psvector & dataset$psvector<=quartilebounds[k+1]] = k
    n1[k] = sum(dataset$quartile==k & dataset$treatment==1,na.rm=TRUE)
    n0[k] = sum(dataset$quartile==k & dataset$treatment==0,na.rm=TRUE)
    
    #identify strata with no datapoints
    if(n1[k]==0 | n0[k]==0){
      goodstrata[k] = NA} else{
      strata.TE[k] = mean(dataset$outcome[dataset$quartile==k & dataset$treatment==1],na.rm=TRUE) -
        mean(dataset$outcome[dataset$quartile==k & dataset$treatment==0],na.rm=TRUE) 
      all.cov.balance[k] =  covar.bal(dataset[dataset$quartile==k & dataset$treatment==1,],
                              dataset[dataset$quartile==k & dataset$treatment==0,],
                              names(dataset)[!names(dataset)%in%c("outcome","treatment","psvector","quartile")],weights1=rep(1,dim(treated)[1]),weights0=rep(1,dim(control)[1]))
    }
    
    
   }
    cov.balance = mean(all.cov.balance,na.rm=TRUE)
    
    datasetp = dataset[dataset$quartile %in% goodstrata,]
    #change quintiles to indicator variables
    datasetp$quartile = as.factor(datasetp$quartile)
    proportions = as.matrix(table(datasetp$quartile)[-1]/n)
    
  
######################## BPSA-1 : MLE linear regression with interactions
    if(bpsa==1){
      
      fit = lm(outcome ~treatment*quartile,data=datasetp)
      A = matrix(c(1,proportions),nrow=1)
      indicies = (sum(!is.na(goodstrata))+2):dim(summary(fit)$coef)[1]
      estimated.effect = A%*%matrix(summary(fit)$coef[c(2,indicies),1],ncol=1)
       
      covar.matrix = vcov(fit)[c(2,indicies),c(2,indicies)]
      estimated.se = sqrt(A%*%covar.matrix%*%t(A))
      
    } else{ #####################################BPSA-2
         #perform bayesian linear regression 
      tryCatch({
      posteriors = as.matrix(MCMCregress(outcome~treatment*quartile,data=datasetp,burnin=S*5,mcmc=S*50,thin=50))
      posteriors = data.frame(posteriors)
      relevantcovar = as.matrix(posteriors[, grepl( "treatment." , names( posteriors ) ) ])
      ates = posteriors$treatment + t(relevantcovar%*%proportions)
      
      estimated.effect = mean(ates)
      estimated.se = sd(ates)

      } , error=function(e){})
      
    }
    
    percent.data.used = sum(n1[goodstrata],n0[goodstrata])/n 
    nu = dataset$quartile
    strata.prop = n1/(n/strat.cut) #proportion of treated observations in each strata

#################################################################################
} else{return("Invalid method")}

  return(list("ATE" = estimated.effect,"Average SE" = estimated.se,
              "PercDataUsed" = percent.data.used,"Nu" = nu, "Ntreat" = n.treat, 
              "Ncontrol" = n.control,
              "Prop.Treated.Strata" = strata.prop, "Indi.Strata.TE" = strata.TE,
              "CovariateBalance" = cov.balance,"Sandwich SE" = sandwich.se))
}