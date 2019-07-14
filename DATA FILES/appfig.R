######Application graphs
require(ggplot2)

###################################################
## Figure 4

load("~/Dropbox/Shirley's Dissertation/Non-Github folder/1st paper/Final Code/Application/Code for application/app.data.R")
source('~/Dropbox/Shirley\'s Dissertation/Non-Github folder/1st paper/Final Code/Simulation/Run code/designfunction.R')

PS.matrix = PS.design(app.data[[1]]$treatment,data.frame(app.data[[1]]$covariates),FALSE)

PSandtreat = data.frame(treatment=app.data[[1]]$treatment,psvector=PS.matrix)
PSandtreat$treatment = as.factor(PSandtreat$treatment)

tiff(file = "~/Dropbox/Shirley's Dissertation/Non-Github folder/1st paper/Final figures/app.hist.tiff")
ggplot(PSandtreat, aes(x=psvector, fill=treatment)) +
  geom_histogram(alpha=0.2, position="identity")  +
  labs(title="MLE PS estimates for treated and control", x = "PS",
       y = "Count")
dev.off()


fullappdata = data.frame(PSandtreat,data.frame(app.data[[1]]$covariates))
xtable(t(aggregate(.~treatment,fullappdata,mean)))

#################################################################################

## Figure 5

load("~/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Code for application/BPSA.app.R")
load("~/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Code for application/PSA.app.R")

##add lower and upper intervals
PSA.sim.results = data.frame(lapply(psa.result,unlist))
#names(PSA.sim.results) = rownames(psa.result)
BPSA.sim.results = data.frame(bpsa.result)

PSA.sim.results$interval.lower = PSA.sim.results$ATE - 1.96*PSA.sim.results$se
PSA.sim.results$interval.upper = PSA.sim.results$ATE + 1.96*PSA.sim.results$se
BPSA.sim.results$interval.lower = BPSA.sim.results$ATE - 1.96*BPSA.sim.results$rubinse
BPSA.sim.results$interval.upper = BPSA.sim.results$ATE + 1.96*BPSA.sim.results$rubinse


##make combined data.frame

BPSA.data = data.frame(ATE = BPSA.sim.results$ATE,interval.lower = BPSA.sim.results$interval.lower,
                       interval.upper = BPSA.sim.results$interval.upper,method = BPSA.sim.results$method,
                       analysis = "BPSA",pu = BPSA.sim.results$percused)

PSA.data = data.frame(ATE = PSA.sim.results$ATE,interval.lower = PSA.sim.results$interval.lower,
                       interval.upper = PSA.sim.results$interval.upper,method = PSA.sim.results$method,
                       analysis = "PSA",pu=PSA.sim.results$percused)

alldata = rbind(BPSA.data,PSA.data)

#add ate and interval from unadjusted analysis
data = data.frame(treatment = app.data[[1]]$treatment,outcome = app.data[[1]]$outcome,app.data[[1]]$covariates)
unadj.anal = lm(outcome~treatment,data)
adj.anal = lm(outcome~.,data)

ua.row = data.frame(ATE = coef(unadj.anal)[2],interval.lower=coef(unadj.anal)[2]-1.96*summary(unadj.anal)$coef[2,2],
           interval.upper=coef(unadj.anal)[2]+1.96*summary(unadj.anal)$coef[2,2],
           method="regression",analysis="UNADJ",pu=1)

aa.row = data.frame(ATE = coef(adj.anal)[2],interval.lower=coef(adj.anal)[2]-1.96*summary(adj.anal)$coef[2,2],
           interval.upper=coef(adj.anal)[2]+1.96*summary(adj.anal)$coef[2,2],
           method="regression",analysis="ADJ",pu=1)

alldata = rbind(alldata,ua.row,aa.row)

##plot all 
alldata = alldata[!(alldata$method=="caliper.matching.5"),]
imp.names = c('caliper.matching' = "Caliper", "dr" = "DR", "ipw" = "IPW",
              "nn.matching" = "NN", "strat" = "Strat",'regression'="Reg")
alldata$method = factor(alldata$method,levels = c("caliper.matching","nn.matching","strat","ipw","dr","regression"))


tiff(file = "~/Dropbox/Shirley's Dissertation/Paper1/Figures/app_results.tiff")
#png(file = "~/Dropbox/Shirley's Dissertation/Paper1/Figures/app_results.png")

ggplot(data=alldata,aes(x=analysis, y=ATE)) + 
  geom_point(size=0.5) + 
  geom_errorbar(data=alldata,aes(x=analysis,ymin = interval.lower, ymax = interval.upper))+
  facet_grid(.~method,scales="free",labeller = as_labeller(imp.names)) +
  labs(title="ATE and 95% CI", 
                                  x = "Analysis method", y = "Difference in micro grams/cubic meter of PM") + 
  theme(legend.position="top") +
  theme(axis.text=element_text(size=12),
        #                                axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        strip.text.x = element_text(size = 10))+ 
  scale_linetype_discrete(name = "SE", 
                          labels = c("Conditional","Empirical","Marginal")) +
  scale_color_discrete(name = "SE", 
                       labels = c("Conditional","Empirical","Marginal"))

dev.off()
