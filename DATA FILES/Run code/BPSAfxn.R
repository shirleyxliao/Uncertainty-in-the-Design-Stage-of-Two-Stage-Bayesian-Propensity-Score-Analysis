#wrapper function for BPSA
require(MASS)
require(survey)

###rubin's se function
library(parallel)

# Calculate the number of cores
no_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
#no_cores = 4

# Initiate cluster
cl <- makeCluster(no_cores,type="FORK",outfile="")


#function
BPSA.fxn = function(input.data,foldername){
  
  nmethods = 6  #strat,nn, caliper, caliper.5, ipw, dr
  K =1000 #MCMC iterations
  n=length(input.data[[1]][["outcome"]]) #sample size
  p=dim(input.data[[1]][["covariates"]])[2] #number of covariates
  M=length(input.data) - 1 #number of datasets 
  
  ###Create data structures to keep results in
  BPSA.est.PS = list() #holds estimated PS
  
  #BPSA.sim.results = data.frame(ATE = rep(NA,M*nmethods))
  BPSA.sim.results = data.frame(ATE = rep(NA,nmethods))
  #PS.sd = matrix(NA,nrow=n,ncol=M)
  #strata.PROP.list = list()
  #strata.TE.list = list()

  #iterate through the datasets
  parSapply(cl,as.character(1:M),function(m,input.data){
    rubinse = function(ates,sigmas){
      K = length(ates)
      sqrt(mean(sigmas^2,na.rm=TRUE) + var(ates,na.rm=TRUE))
    }
    
    #call on design and analysis fxns
    source(paste('/n/home10/silverfoxflower/',foldername,'/designfunction.R',sep=""))
    source(paste('/n/home10/silverfoxflower/',foldername,'/analysisfunction.R',sep=""))
    
    
    m = as.numeric(m)
  #for(m in 1:M){
    #grab simulated treatment/outcome/covariate information
    treatment = input.data[[m]][["treatment"]]
    outcome = input.data[[m]][["outcome"]]
    covardataset = input.data[[m]][["covariates"]]
    
    ###bayesian PS estimation
    BPSA.est.PS = PS.design(treatment,covardataset,bayes=TRUE,K=K)
    #PS.sd[,m] = apply(BPSA.est.PS[[m]],1,sd)
    
    ##data objects
    caliper.matchingATE = rep(NA,K)
    caliper.matching.5ATE = rep(NA,K)
    ipwATE = rep(NA,K)
    nn.matchingATE = rep(NA,K)
    stratATE = rep(NA,K)
    regATE = rep(NA,K)
    drATE = rep(NA,K)
    
    caliper.matchingSE = rep(NA,K)
    caliper.matching.5SE = rep(NA,K)
    ipwSE = rep(NA,K)
    nn.matchingSE = rep(NA,K)
    stratSE = rep(NA,K)
    regSE = rep(NA,K)
    drSE = rep(NA,K)
    
    caliper.matchingSANDWICHSE = rep(NA,K)
    caliper.matching.5SANDWICHSE = rep(NA,K)
    ipwSANDWICHSE = rep(NA,K)
    nn.matchingSANDWICHSE = rep(NA,K)

    caliper.matchingPU = rep(NA,K)
    caliper.matching.5PU = rep(NA,K)
    ipwPU = rep(NA,K)
    stratPU = rep(NA,K)
    regPU = rep(NA,K)
    nn.matchingPU = rep(NA,K)
    drPU = rep(NA,K)
    
    caliper.matchingNU = matrix(NA,nrow=n,ncol=K)
    caliper.matching.5NU = matrix(NA,nrow=n,ncol=K)
    ipwNU = matrix(NA,nrow=n,ncol=K)
    stratNU = matrix(NA,nrow=n,ncol=K)
    regNU = matrix(NA,nrow=n,ncol=K)
    nn.matchingNU =matrix(NA,nrow=n,ncol=K)
    drNU = matrix(NA,nrow=n,ncol=K)
    
    caliper.matchingNM= matrix(NA,nrow=K,ncol=2)
    caliper.matching.5NM= matrix(NA,nrow=K,ncol=2)
    nn.matchingNM = matrix(NA,nrow=K,ncol=2)

    
    #strataPROP = matrix(NA,nrow=5,ncol=K)
    
    
    #caliper.matchingCB = rep(NA,K)
    #caliper.matching.5CB = rep(NA,K)
    #ipwCB = rep(NA,K)
    #stratCB = rep(NA,K)
    #regCB = rep(NA,K)
    #nn.matchingCB = rep(NA,K)
    #drCB = rep(NA,K)
    
    
    for(k in 1:K){
      
      ###estimation via 1-1 nn matching with replacement
      nn.matching = PS.analysis(BPSA.est.PS[,k],treatment,outcome,covardataset,"nn.matching",sandwich.est=FALSE)
      
      nn.matchingATE[k] = nn.matching$ATE
      nn.matchingSE[k] = nn.matching[["Average SE"]]
      nn.matchingPU[k] = nn.matching$PercDataUsed
      nn.matchingNU[,k] = nn.matching$Nu
      nn.matchingNM[k,] = c(nn.matching$Ncontrol,nn.matching$Ntreat)
      #nn.matchingCB[k] = nn.matching$CovariateBalance
      nn.matchingSANDWICHSE[k] = nn.matching[["Sandwich SE"]]
      
      ###1-1 caliper matching with replacement (BPSA-2)
      caliper.matching = PS.analysis(BPSA.est.PS[,k],treatment,outcome,covardataset,"caliper.matching",bpsa=2,matching.ratio=1,sandwich.est=FALSE)
      
      caliper.matchingATE[k] = mean(caliper.matching$ATE)
      caliper.matchingSE[k] = caliper.matching[["Average SE"]]
      caliper.matchingPU[k] = caliper.matching$PercDataUsed
      caliper.matchingNU[,k] = caliper.matching$Nu
      caliper.matchingNM[k,] = c(caliper.matching$Ncontrol,caliper.matching$Ntreat)
      #caliper.matchingCB[k] = caliper.matching$CovariateBalance
      caliper.matchingSANDWICHSE[k] = caliper.matching[["Sandwich SE"]]
      
      ###1-5 caliper matching with replacement (BPSA-2)
      caliper.matching.5 = PS.analysis(BPSA.est.PS[,k],treatment,outcome,covardataset,"caliper.matching",bpsa=2,matching.ratio=5,sandwich.est=FALSE)
      
      caliper.matching.5ATE[k] = mean(caliper.matching.5$ATE)
      caliper.matching.5SE[k] = caliper.matching.5[["Average SE"]]
      caliper.matching.5PU[k] = caliper.matching.5$PercDataUsed
      caliper.matching.5NU[,k] = caliper.matching.5$Nu
      caliper.matching.5NM[k,] = c(caliper.matching.5$Ncontrol,caliper.matching.5$Ntreat)
      #caliper.matching.5CB[k] = caliper.matching.5$CovariateBalance
      caliper.matching.5SANDWICHSE[k] = caliper.matching.5[["Sandwich SE"]]
            
      ###estimation via DR
      dr = PS.analysis(BPSA.est.PS[,k],treatment,outcome,covardataset,"dr")
      
      drATE[k] = dr$ATE
      drSE[k] = dr[["Average SE"]]
      drPU[k] = dr$PercDataUsed
      drNU[,k] = dr$Nu
      #drCB[k] = dr$CovariateBalance
      
      ###estimation via IPW
      ipw = PS.analysis(BPSA.est.PS[,k],treatment,outcome,covardataset,"ipw",sandwich.est=FALSE)
      
      ipwATE[k] = ipw$ATE
      ipwSE[k] = ipw[["Average SE"]]
      ipwPU[k] = ipw$PercDataUsed
      ipwNU[,k] = ipw$Nu
      ipwSANDWICHSE[k] = ipw[["Sandwich SE"]]
      #ipwCB[k] = ipw$CovariateBalance
      
      ###estmation via stratification
      strat = PS.analysis(BPSA.est.PS[,k],treatment,outcome,covardataset,method="stratification",bpsa=2)
      
      stratATE[k] = strat$ATE
      stratSE[k] = strat[["Average SE"]]
      stratPU[k] = strat$PercDataUsed
      stratNU[,k] = strat$Nu
      #strataPROP[,k] = strat$Prop.Treated.Strata
      #stratCB[k] = strat$CovariateBalance
    }
    #for the purpose of loading information into the BPSA.sim.results dataframe 
    #startit = 1+(nmethods*(m-1))
    #endit = startit+nmethods-1
    
    #load data into BPSA.sim.results
   # BPSA.sim.results$MCMC.SE = c(sd(nn.matchingMCMC,na.rm=TRUE),
   #                                       sd(matchingMCMC,na.rm=TRUE),
   #                                       sd(drMCMC,na.rm=TRUE),
   #                                       sd(ipwMCMC,na.rm=TRUE),
   #                                       sd(stratMCMC,na.rm=TRUE))
                                      
    
    BPSA.sim.results$ATE = c(mean(nn.matchingATE,na.rm=TRUE),
                                             mean(caliper.matchingATE,na.rm=TRUE),
                                            mean(caliper.matching.5ATE,na.rm=TRUE), 
                                            mean(drATE,na.rm=TRUE),
                                             mean(ipwATE,na.rm=TRUE),
                                             mean(stratATE,na.rm=TRUE))
                                        
    BPSA.sim.results$between.design.SE = c(sd(nn.matchingATE,na.rm=TRUE),
                             sd(caliper.matchingATE,na.rm=TRUE),
                             sd(caliper.matching.5ATE,na.rm=TRUE), 
                             sd(drATE,na.rm=TRUE),
                             sd(ipwATE,na.rm=TRUE),
                             sd(stratATE,na.rm=TRUE))
    
    BPSA.sim.results$within.design.SE = c(mean(nn.matchingSE,na.rm=TRUE),
                                                   mean(caliper.matchingSE,na.rm=TRUE),
                                                   mean(caliper.matching.5SE,na.rm=TRUE),
                                                   mean(drSE,na.rm=TRUE),
                                        mean(ipwSE,na.rm=TRUE),
                                        mean(stratSE,na.rm=TRUE))
    
    BPSA.sim.results$within.design.sandwich.SE = c(mean(nn.matchingSANDWICHSE,na.rm=TRUE),
                                    mean(caliper.matchingSANDWICHSE,na.rm=TRUE),
                                    mean(caliper.matching.5SANDWICHSE,na.rm=TRUE),
                                    NA,
                                    mean(ipwSANDWICHSE,na.rm=TRUE),
                                    NA)
                          
    
    BPSA.sim.results$method = c("nn.matching",
                                           "caliper.matching",
                                           "caliper.matching.5",
                                           "dr",
                                           "ipw","strat")
    
    BPSA.sim.results$rubinse = c(rubinse(nn.matchingATE,nn.matchingSE),
                                            rubinse(caliper.matchingATE,caliper.matchingSE),
                                            rubinse(caliper.matching.5ATE,caliper.matchingSE),
                                            rubinse(drATE,drSE),
                                            rubinse(ipwATE,ipwSE),
                                            rubinse(stratATE,stratSE))
                                         
    
    BPSA.sim.results$interval.lower = BPSA.sim.results$ATE - 1.96*BPSA.sim.results$rubinse
                                              
    
    BPSA.sim.results$interval.upper = BPSA.sim.results$ATE + 1.96*BPSA.sim.results$rubinse
                                                
    
    BPSA.sim.results$percused =c(mean(nn.matchingPU,na.rm=TRUE),
                                            mean(caliper.matchingPU,na.rm=TRUE),
                                            mean(caliper.matching.5PU,na.rm=TRUE),
                                            mean(drPU,na.rm=TRUE),
                                            mean(ipwPU,na.rm=TRUE), 
                                            mean(stratPU,na.rm=TRUE))
  
    
 #   BPSA.sim.results$num.na =c(mean(is.na(nn.matchingATE)),
  #                                        mean(is.na(caliper.matchingATE)),
   #                                       mean(is.na(drATE)),
   #                                       mean(is.na(ipwATE)),
   #                                       mean(is.na(stratATE)))
    
    BPSA.sim.results$marginal.SD.nu = c(mean(apply(nn.matchingNU,1,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(caliper.matchingNU,1,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(caliper.matching.5NU,1,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(drNU,1,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(ipwNU,1,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(stratNU,1,sd,na.rm=TRUE),na.rm=TRUE))
    
    BPSA.sim.results$marginal.SD.nu.treated = c(mean(apply(nn.matchingNU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(caliper.matchingNU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(caliper.matching.5NU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(drNU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(ipwNU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(stratNU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE))
    
    BPSA.sim.results$marginal.SD.nu.control =  c(mean(apply(nn.matchingNU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(caliper.matchingNU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(caliper.matching.5NU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(drNU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(ipwNU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(stratNU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE))
    
    BPSA.sim.results$ncontrol = c(mean(nn.matchingNM[,1],na.rm=TRUE),
                                                 mean(caliper.matchingNM[,1],na.rm=TRUE),mean(caliper.matching.5NM[,1],na.rm=TRUE),
                                                 NA,NA,NA)
    
    BPSA.sim.results$ntreat = c(mean(nn.matchingNM[,2],na.rm=TRUE),
                                               mean(caliper.matchingNM[,2],na.rm=TRUE),
                                               mean(caliper.matching.5NM[,2],na.rm=TRUE),
                                               NA,NA,NA)
    
    #BPSA.sim.results$cov.balance =c(mean(nn.matchingCB,na.rm=TRUE),
    #                                mean(caliper.matchingCB,na.rm=TRUE),
    #                                mean(caliper.matching.5CB,na.rm=TRUE),
     #                               mean(drCB,na.rm=TRUE),
     #                               mean(ipwCB,na.rm=TRUE), 
      #                              mean(stratCB,na.rm=TRUE))
    
    #strata.PROP.list[[m]] = strataPROP
    
    return(BPSA.sim.results)
    
  },input.data = input.data) 
  
  #BPSA.output = list(BPSA.sim.results,nn.matchingNU,caliper.matchingNU,caliper.matching.5NU,ipwNU,stratNU,PS.sd,strata.PROP.list)
  #save(BPSA.output,file = paste("/n/home10/silverfoxflower/",foldername,"/",savefile,".R",sep=""))
}