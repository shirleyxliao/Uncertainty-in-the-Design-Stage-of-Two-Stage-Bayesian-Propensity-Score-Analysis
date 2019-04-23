#load packages + required functions
require(MASS)

expit = function(x){
  exp(x)/(1 + exp(x))
}

###########
n =1000 #sample size
M = 500 #replications

instru = 5 #number of instruments
prog = 5 #number of prognostic variables
noise = 5
confound = 5
p = confound+instru+prog+noise

identify = c(rep(0,confound),rep(1,instru),rep(2,prog),rep(3,noise))
alpha = matrix(c(0,rep(.75,confound),rep(0.75,instru),rep(0,prog+noise)),nrow=1) #covariates for PS model
beta = matrix(c(0.1,seq(0.2,0.4,length=confound),rep(0,instru),rep(0.5,prog),rep(0,noise)),nrow=1) #covariates for outcome model
D = diag(rep(1,p)) #correlation matrix for covaritates: all covariates are independent
sigma = 1 #variance for error term of the outcome model

#empty data structures

truePS = matrix(NA,n,M)
outcomes = matrix(NA,n,M)
treatments = matrix(NA,n,M)

continuousdata = list()

##iterate to create each dataset
for(m in 1:M){

###draw covariates from a multivariate normal distribution
X = matrix( mvrnorm(n,mu=rep(0,p),Sigma=D), ncol = p)
X = cbind(rep(1,n),X) #intercept

#true propensity score calculation
truePS[,m] = expit(alpha%*%t(X)) 

#Treatment drawn based on a Bernoulli with probability = truePS
treatmentmatrix = matrix(rbinom(n,1,truePS[,m]),nrow=n,ncol=1)
treatments[,m] = treatmentmatrix

###continuous outcome simulation
error = rnorm(n,0,sigma) #random error term for outcome
trueATE = 1.5 
outcome = trueATE*treatmentmatrix + X%*%t(beta)+error
outcomes[,m] = outcome

#save simulated data
continuousdata[[m]] = list()

continuousdata[[m]][["outcome"]] = outcome
X = X[,-1]
continuousdata[[m]][["covariates"]] = X[,1:(confound)] #vary this to include instruments in the analysis
continuousdata[[m]][["treatment"]] = treatmentmatrix
continuousdata[[m]][["truePS"]] = truePS[,m]

}


continuousdata[["trueATE"]] = trueATE

save(continuousdata,file="/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Data generation/instru_prog_noise_data.R")

####
