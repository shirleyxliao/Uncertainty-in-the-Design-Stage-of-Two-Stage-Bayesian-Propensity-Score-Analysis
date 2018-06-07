#wrapper function for BPSA
require(MASS)
require(survey)

###rubin's se function
rubinse = function(ates,sigmas){
  K = length(ates)
  sqrt(mean(sigmas^2,na.rm=TRUE) + var(ates,na.rm=TRUE))
}

#wrapper function which calls analysisfxn
BPSA.fxn = function(input.data,savefile,foldername){
  nmethods = 6  #caliper matching with 1, 2, 5 ratios and nn matching with 1, 2, 5 ratios
  K =200 #MCMC iterations
  n=length(input.data[[1]][["outcome"]]) #sample size
  p=dim(input.data[[1]][["covariates"]])[2] #number of covariates
  M=length(input.data) - 1 #number of datasets 
  S = 200 #number of draws taken from the posterior distribution of delta
  R = 15 #number of times caliper matching is repeated on the same PS
  
  #call on design and analysis fxns
  source(paste('/n/home10/silverfoxflower/',foldername,'/designfunction.R',sep=""))
  source(paste('/n/home10/silverfoxflower/',foldername,'/analysisfunction.R',sep=""))
  
  ###Create data structures to keep results in
  BPSA.est.PS = list() #holds estimated PS
  BPSA.sim.results = data.frame(ATE = rep(NA,M*nmethods)) #holds summary stats
  PS.sd = matrix(NA,nrow=n,ncol=M) #holds SD of the estimated PS
  
  #iterate through the datasets
  for(m in 1:M){
    #grab simulated treatment/outcome/covariate information
    treatment = input.data[[m]][["treatment"]]
    outcome = input.data[[m]][["outcome"]]
    covardataset = input.data[[m]][["covariates"]]
    
    ###bayesian PS estimation
    BPSA.est.PS[[m]] = PS.design(treatment,covardataset,bayes=TRUE,K=K)
    PS.sd[,m] = apply(BPSA.est.PS[[m]],1,sd) #calculate SD of PS within each observation
    
    ##data objects
    #ATE
    #NOTE: for caliper matching, ATEs are outputted as SxR draws from the posterior
    ##distribution of delta, for nn matching, ATEs are outputted as the mean
    ###of the posterior distribution
    caliper.matching.5ATE = array(NA,dim=c(K,S,R))
    nn.matching.5ATE = rep(NA,K)
    caliper.matching.2ATE = array(NA,dim=c(K,S,R))
    nn.matching.2ATE = rep(NA,K)
    caliper.matching.1ATE = array(NA,dim=c(K,S,R))
    nn.matching.1ATE = rep(NA,K)
    
    #avg SE
    caliper.matching.5SE = rep(NA,K)
    nn.matching.5SE = rep(NA,K)
    caliper.matching.2SE = rep(NA,K)
    nn.matching.2SE = rep(NA,K)
    caliper.matching.1SE = rep(NA,K)
    nn.matching.1SE = rep(NA,K)
    
    #percent of data used
    caliper.matching.5PU = rep(NA,K)
    nn.matching.5PU = rep(NA,K)
    caliper.matching.2PU = rep(NA,K)
    nn.matching.2PU = rep(NA,K)
    caliper.matching.1PU = rep(NA,K)
    nn.matching.1PU = rep(NA,K)
    
    #saves nu in a matrix or array
    caliper.matching.5NU = array(NA,dim=c(K,n,R))
    nn.matching.5NU = matrix(NA,nrow=n,ncol=K)
    caliper.matching.2NU = array(NA,dim=c(K,n,R))
    nn.matching.2NU = matrix(NA,nrow=n,ncol=K)
    caliper.matching.1NU = array(NA,dim=c(K,n,R))
    nn.matching.1NU = matrix(NA,nrow=n,ncol=K)
    
    #number matched, treated/control
    caliper.matching.5NM= matrix(NA,nrow=K,ncol=2)
    nn.matching.5NM = matrix(NA,nrow=K,ncol=2)
    caliper.matching.2NM= matrix(NA,nrow=K,ncol=2)
    nn.matching.2NM = matrix(NA,nrow=K,ncol=2)
    caliper.matching.1NM= matrix(NA,nrow=K,ncol=2)
    nn.matching.1NM = matrix(NA,nrow=K,ncol=2)
    
    for(k in 1:K){
      
      ###estimation via nn matching with replacement (ratio = 5)
      nn.matching = PS.analysis(BPSA.est.PS[[m]][,k],treatment,outcome,
                                covardataset,"nn.matching",bpsa=2,
                                matching.ratio=5)
      
      nn.matching.5ATE[k] = nn.matching$ATE
      nn.matching.5SE[k] = nn.matching[["Average SE"]]
      nn.matching.5PU[k] = nn.matching$PercDataUsed
      nn.matching.5NU[,k] = nn.matching$Nu
      nn.matching.5NM[k,] = c(nn.matching$Ncontrol,nn.matching$Ntreat)
      
      ###estimation via nn matching with replacement (ratio = 1)
      nn.matching = PS.analysis(BPSA.est.PS[[m]][,k],treatment,outcome,
                                covardataset,"nn.matching",bpsa=2,
                                matching.ratio=1)
      
      nn.matching.1ATE[k] = nn.matching$ATE
      nn.matching.1SE[k] = nn.matching[["Average SE"]]
      nn.matching.1PU[k] = nn.matching$PercDataUsed
      nn.matching.1NU[,k] = nn.matching$Nu
      nn.matching.1NM[k,] = c(nn.matching$Ncontrol,nn.matching$Ntreat)
      
   ## estimation via nn matching with replacement (ratio = 2)
      nn.matching = PS.analysis(BPSA.est.PS[[m]][,k],treatment,outcome,covardataset,
                                "nn.matching",bpsa=2,matching.ratio=2)
      
      nn.matching.2ATE[k] = nn.matching$ATE
      nn.matching.2SE[k] = nn.matching[["Average SE"]]
      nn.matching.2PU[k] = nn.matching$PercDataUsed
      nn.matching.2NU[,k] = nn.matching$Nu
      nn.matching.2NM[k,] = c(nn.matching$Ncontrol,nn.matching$Ntreat)
      
      ###estimation via caliper matching with replacement (ratio = 5)
      caliper.matching = PS.analysis(BPSA.est.PS[[m]][,k],treatment,
                                     outcome,covardataset,"caliper.matching",
                                     match.iter=R,bpsa=2,matching.ratio=5)
      
      caliper.matching.5ATE[k,,] = caliper.matching$ATE
      caliper.matching.5SE[k] = caliper.matching[["Average SE"]]
      caliper.matching.5PU[k] = caliper.matching$PercDataUsed
      caliper.matching.5NU[k,,] = caliper.matching$Nu
      caliper.matching.5NM[k,] = c(caliper.matching$Ncontrol,caliper.matching$Ntreat)
      

      #estimation via caliper matching with replacement (ratio = 1)
      caliper.matching = PS.analysis(BPSA.est.PS[[m]][,k],treatment,outcome,
                                     covardataset,"caliper.matching",
                                     match.iter=R,bpsa=2,matching.ratio=1)
      
      caliper.matching.1ATE[k,,] = caliper.matching$ATE
      caliper.matching.1SE[k] = caliper.matching[["Average SE"]]
      caliper.matching.1PU[k] = caliper.matching$PercDataUsed
      caliper.matching.1NU[k,,] = caliper.matching$Nu
      caliper.matching.1NM[k,] = c(caliper.matching$Ncontrol,caliper.matching$Ntreat)
    

      #estimation via caliper matching with replacement (ratio = 2)
      caliper.matching = PS.analysis(BPSA.est.PS[[m]][,k],treatment,outcome,
                                     covardataset,"caliper.matching",match.iter=R,
                                     matching.ratio=2,bpsa=2)
      
      caliper.matching.2ATE[k,,] = caliper.matching$ATE
      caliper.matching.2SE[k] = caliper.matching[["Average SE"]]
      caliper.matching.2PU[k] = caliper.matching$PercDataUsed
      caliper.matching.2NU[k,,] = caliper.matching$Nu
      caliper.matching.2NM[k,] = c(caliper.matching$Ncontrol,caliper.matching$Ntreat)
      
          }
    #for the purpose of loading information into the BPSA.sim.results dataframe 
    startit = 1+(nmethods*(m-1))
    endit = startit+nmethods-1
  
    #ATE
    BPSA.sim.results$ATE[startit:endit] = c(mean(nn.matching.5ATE,na.rm=TRUE),
                                            mean(caliper.matching.5ATE,na.rm=TRUE),
                                            mean(nn.matching.1ATE,na.rm=TRUE),
                                             mean(caliper.matching.1ATE,na.rm=TRUE),
                                             mean(nn.matching.2ATE,na.rm=TRUE),
                                              mean(caliper.matching.2ATE,na.rm=TRUE))

    #average conditional SE
    BPSA.sim.results$avg.est.SE[startit:endit] = c(mean(nn.matching.5SE,na.rm=TRUE),
                                                   mean(caliper.matching.5SE,na.rm=TRUE),
                                                   mean(nn.matching.1SE,na.rm=TRUE),
                                                   mean(caliper.matching.1SE,na.rm=TRUE),
                                                   mean(nn.matching.2SE,na.rm=TRUE),
                                                   mean(caliper.matching.2SE,na.rm=TRUE))
    
    #method
    BPSA.sim.results$method[startit:endit] = c("nn.matching",
                                               "caliper.matching","nn.matching",
                                               "caliper.matching","nn.matching",
                                               "caliper.matching")
    
    #ratio
    BPSA.sim.results$ratio[startit:endit] = c(5,5,1,1,2,2)
                                               
    
    #marginal SE calculated via rubin's method or empirical SD of KxRxS draws
    ##from the posterior distribution of delta
    BPSA.sim.results$rubinse[startit:endit] = c(rubinse(nn.matching.5ATE,nn.matching.5SE),
                                                sd(caliper.matching.5ATE),
                                                rubinse(nn.matching.1ATE,nn.matching.1SE),
                                                sd(caliper.matching.1ATE),
                                                rubinse(nn.matching.2ATE,nn.matching.2SE),
                                                sd(caliper.matching.2ATE))
    
    #lower bound of credible interval
    BPSA.sim.results$interval.lower[startit:endit] = BPSA.sim.results$ATE[startit:endit] - 1.96*BPSA.sim.results$rubinse[startit:endit]
    
    #upper bound of credible interval
    BPSA.sim.results$interval.upper[startit:endit] = BPSA.sim.results$ATE[startit:endit] + 1.96*BPSA.sim.results$rubinse[startit:endit]
    
    #percent of data used
    BPSA.sim.results$percused[startit:endit] =c(mean(nn.matching.5PU,na.rm=TRUE),
                                                mean(caliper.matching.5PU,na.rm=TRUE),
                                                mean(nn.matching.1PU,na.rm=TRUE),
                                                mean(caliper.matching.1PU,na.rm=TRUE),
                                                mean(nn.matching.2PU,na.rm=TRUE),
                                                mean(caliper.matching.2PU,na.rm=TRUE))
    
    #avg SD of the marginal distribution of nu
    BPSA.sim.results$marginal.SD.nu[startit:endit] = c(mean(apply(nn.matching.5NU,1,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(caliper.matching.5NU,2,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(nn.matching.1NU,1,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(caliper.matching.1NU,2,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(nn.matching.2NU,1,sd,na.rm=TRUE),na.rm=TRUE),
                                                       mean(apply(caliper.matching.2NU,2,sd,na.rm=TRUE),na.rm=TRUE))
    
    #avg SD of the marginal distribution of nu for treated observations                                                  
    BPSA.sim.results$marginal.SD.nu.treated[startit:endit] = c(mean(apply(nn.matching.5NU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(caliper.matching.5NU[,treatment==1,],2,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(nn.matching.1NU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(caliper.matching.1NU[,treatment==1,],2,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(nn.matching.2NU[treatment==1,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                               mean(apply(caliper.matching.2NU[,treatment==1,],2,sd,na.rm=TRUE),na.rm=TRUE))
    
    #avg SD of the marginal distribution of nu for control observations
    BPSA.sim.results$marginal.SD.nu.control[startit:endit] =  c(mean(apply(nn.matching.5NU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(caliper.matching.5NU[,treatment==0,],2,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(nn.matching.1NU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(caliper.matching.1NU[,treatment==0,],2,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(nn.matching.2NU[treatment==0,],1,sd,na.rm=TRUE),na.rm=TRUE),
                                                                mean(apply(caliper.matching.2NU[,treatment==0,],2,sd,na.rm=TRUE),na.rm=TRUE))
    
    #avg number of control observations in the matched set
     BPSA.sim.results$ncontrol[startit:endit] = c(mean(nn.matching.5NM[,1],na.rm=TRUE),
                                               mean(caliper.matching.5NM[,1],na.rm=TRUE),
                                               mean(nn.matching.1NM[,1],na.rm=TRUE),
                                               mean(caliper.matching.1NM[,1],na.rm=TRUE),
                                               mean(nn.matching.2NM[,1],na.rm=TRUE),
                                               mean(caliper.matching.2NM[,1],na.rm=TRUE))
    
     #avg number of treated observations in the matched set
    BPSA.sim.results$ntreat[startit:endit] = c(mean(nn.matching.5NM[,2],na.rm=TRUE),
                                                 mean(caliper.matching.5NM[,2],na.rm=TRUE),
                                                 mean(nn.matching.1NM[,2],na.rm=TRUE),
                                                 mean(caliper.matching.1NM[,2],na.rm=TRUE),
                                                 mean(nn.matching.2NM[,2],na.rm=TRUE),
                                                 mean(caliper.matching.2NM[,2],na.rm=TRUE))
    
    ##writes out new file at every iteration
    #BPSA.output = list(BPSA.sim.results,caliper.matching.5NU,caliper.matching.1NU,caliper.matching.2NU,PS.sd)
   # save(BPSA.output,file = paste("/n/home10/silverfoxflower/",foldername,"/",savefile,".R",sep=""))
  } 
  BPSA.output = list(BPSA.sim.results,caliper.matching.1NU,caliper.matching.2NU,caliper.matching.5NU,PS.sd)
  save(BPSA.output,file = paste("/n/home10/silverfoxflower/",foldername,"/",savefile,".R",sep=""))
}