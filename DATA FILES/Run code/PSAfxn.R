###continuousfreq fxn 
require(MASS)
#require(survey)

PSA.fxn = function(input.data,savefile,foldername){

  #load required functions
  source(paste('/n/home10/silverfoxflower/',foldername,'/designfunction.R',sep=""))
  source(paste('/n/home10/silverfoxflower/',foldername,'/analysisfunction.R',sep=""))
  
  #####grab required information from inputted dataset
  nmethods = 6  #strat,nnmatch, caliper, caliper.5, ipw, dr
  n=length(input.data[[1]][["outcome"]]) #sample size
  p=dim(input.data[[1]][["covariates"]])[2] #number of covariates
  M=length(input.data) - 1 #number of datasets 
  
  ###make data structures
  strata.PROP.list = list()
  
  nn.matchingNU = matrix(NA,nrow=n,ncol=M)
  caliper.matchingNU = matrix(NA,nrow=n,ncol=M) 
  caliper.matchingNU.5 = matrix(NA,nrow=n,ncol=M) 
  ipwNU = matrix(NA,nrow=n,ncol=M)
  stratNU = matrix(NA,nrow=n,ncol=M)
  
  PSA.est.PS = matrix(NA,nrow=n,ncol=M) #save estimated PS
  
  PSA.sim.results = data.frame(ATE = rep(NA,M*nmethods))
  PSA.sim.results$method = NA
  PSA.sim.results$percused = NA
  PSA.sim.results$se = NA
  
  for(m in 1:M){
    #grab treatment/outcome/covariate information
    treatment = input.data[[m]][["treatment"]]
    outcome = input.data[[m]][["outcome"]]
    covardataset = input.data[[m]][["covariates"]]
    
    ###frequentist propensity score estimation
    PSA.est.PS[,m] = PS.design(treatment,covardataset,bayes=FALSE)
    
    ###estimation via caliper.matching
    caliper.matching = PS.analysis(PSA.est.PS[,m],treatment,outcome,covardataset,"caliper.matching",sandwich.est=TRUE)

    ###estimation via caliper.matching
    caliper.matching.5 = PS.analysis(PSA.est.PS[,m],treatment,outcome,covardataset,"caliper.matching",matching.ratio=5,sandwich.est=TRUE)
    
    ###estimation via nn.matching's method
    nn.matching = PS.analysis(PSA.est.PS[,m],treatment,outcome,covardataset,"nn.matching",sandwich.est=TRUE)
    
    ###estimation via DR
    dr = PS.analysis(PSA.est.PS[,m],treatment,outcome,covardataset,"dr")
    
    ###estimation via IPW
    ipw = PS.analysis(PSA.est.PS[,m],treatment,outcome,covardataset,"ipw",sandwich.est=TRUE)
    
    ###estmation via stratification
    strat = PS.analysis(PSA.est.PS[,m],treatment,outcome,covardataset,method="stratification")
    
    
    
    #helpful to load data into PSA.sim.results
    startit = 1+(nmethods*(m-1))
    endit = startit+nmethods-1
    
    #load data into PSA.sim.results
    PSA.sim.results$ATE[startit:endit] = c(mean(caliper.matching$ATE),
                                           mean(caliper.matching.5$ATE),
                                           nn.matching$ATE,
                                            dr$ATE,
                                            ipw$ATE,
                                           strat$ATE)

    PSA.sim.results$method[startit:endit] = c("caliper.matching",
                                              "caliper.matching.5",
                                            
                                              "nn.matching",   
                                          "dr",
                                          "ipw",
                                          "strat")
    
    PSA.sim.results$se[startit:endit] = c(caliper.matching[["Average SE"]],
                                          caliper.matching.5[["Average SE"]],
                                      nn.matching[["Average SE"]],
                                       dr[["Average SE"]],
                                      ipw[["Average SE"]],
                                      strat[["Average SE"]])
                                      
    
    PSA.sim.results$percused[startit:endit] = c(caliper.matching$PercDataUsed,
                                                caliper.matching.5$PercDataUsed,
                                            nn.matching$PercDataUsed,
                                              1,1,strat$PercDataUsed)
    
    PSA.sim.results$ncontrol[startit:endit] = c(nn.matching$Ncontrol,
                                                 caliper.matching$Ncontrol,
                                                caliper.matching.5$Ncontrol,
                                                 NA,NA,NA)
    
    PSA.sim.results$ntreat[startit:endit] = c(nn.matching$Ntreat,
                                              caliper.matching$Ntreat,
                                              caliper.matching.5$Ntreat,
                                              NA,NA,NA)
    
    PSA.sim.results$covar.balance[startit:endit] = c(nn.matching$CovariateBalance,
                                                     caliper.matching$CovariateBalance,
                                                     caliper.matching.5$CovariateBalance,
                                      ipw$CovariateBalance,dr$CovariateBalance,strat$CovariateBalance)
    
    nn.matchingNU[,m] = nn.matching$Nu
    caliper.matchingNU[,m] = caliper.matching$Nu
    caliper.matchingNU.5[,m] = caliper.matching.5$Nu
    stratNU[,m] = strat$Nu
    ipwNU[,m] = ipw$Nu
    
    
    strata.PROP.list[[m]] = strat$Prop.Treated.Strata

  }
  PSA.output = list(PSA.sim.results,nn.matchingNU,caliper.matchingNU,caliper.matchingNU.5,ipwNU,stratNU,strata.PROP.list)
  save(PSA.output,file = paste("/n/home10/silverfoxflower/",foldername,"/",savefile,".R",sep=""))
}