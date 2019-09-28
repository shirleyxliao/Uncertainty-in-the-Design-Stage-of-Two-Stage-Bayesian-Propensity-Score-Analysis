################### final simulation graphs
# function to restructure data

data.fxn = function(data){
  ncolumns = dim(data)[2]
  
  for(c in 1:ncolumns){
    if(c==1){
      returnframe = data.frame(data[,c])
      returnframe$rep = 1
    } else{
      xx = data.frame(data[,c])
      xx$rep = c
      returnframe = rbind(returnframe,xx)
    }
  }
  returnframe
}

#################################
# read in results and restructure them

load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Results/BPSA.confound.R")
all.methods = data.frame(data.type = "confound",data.fxn(bpsa.result))

load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Results/BPSA.instru.R")
all.methods = rbind(all.methods,data.frame(data.type = "instru",data.fxn(bpsa.result)))

load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Results/BPSA.instru_prog.R")
all.methods = rbind(all.methods,data.frame(data.type = "instru_prog",data.fxn(bpsa.result)))

load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Results/BPSA.instru_prog_noise.R")
all.methods = rbind(all.methods,data.frame(data.type = "instru_prog_noise",data.fxn(bpsa.result)))

#all.methods$sandwich_se = sqrt(all.methods$between.design.SE^2 + all.methods$within.design.sandwich.SE^2)

all.methods$PROP = (all.methods$between.design.SE)^2/(all.methods$rubinse)^2
all.methods$interval.width = all.methods$interval.upper - all.methods$interval.lower
#all.methods$PROP_sandwich = (all.methods$between.design.SE)^2/(all.methods$sandwich_se)^2

#all.methods$interval.lower = all.methods$ATE - 1.96*all.methods$rubinse
#all.methods$upper = all.methods$ATE + 1.96*all.methods$rubinse

#all.methods$interval.lower_sandwich = all.methods$ATE - 1.96*all.methods$sandwich_se
#all.methods$upper_sandwich = all.methods$ATE + 1.96*all.methods$sandwich_se

all.methods$covered = as.numeric(1.5<all.methods$interval.upper & 1.5>all.methods$interval.lower)
#all.methods$covered_sandwich = as.numeric(1.5<all.methods$upper_sandwich & 1.5>all.methods$lower_sandwich)
all.methods$bias = all.methods$ATE - 1.5

save(all.methods,file="/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Test code/scenario_1_sims_result.R")

results = subset(all.methods,select = c(rep,data.type,ATE,between.design.SE,
                                          within.design.SE,method,rubinse,interval.lower,
                                          interval.upper,covered,bias))
results$interval.width = results$interval.upper - results$interval.lower
results$model = "BPSA"

load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Results/PSA.confound.R")
results.psa = data.frame(data.type = "confound",data.fxn(psa.result))

load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Results/PSA.intru.R")
results.psa = rbind(results.psa,data.frame(data.type = "instru",data.fxn(psa.result)))

load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Results/PSA.intru_prog.R")
results.psa = rbind(results.psa,data.frame(data.type = "instru_prog",data.fxn(psa.result)))

load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Results/PSA.intru_prog_noise.R")
results.psa = rbind(results.psa,data.frame(data.type = "instru_prog_noise",data.fxn(psa.result)))

xx = data.frame(rep = results.psa$rep,data.type = results.psa$data.type, ATE = results.psa$ATE, between.design.SE = 0,  
                within.design.SE = results.psa$se, method = results.psa$method,rubinse = results.psa$se,
                interval.lower = results.psa$ATE - 1.96*results.psa$se, 
                interval.upper = results.psa$ATE + 1.96*results.psa$se,
                #covered = ((interval.lower = results.psa$ATs) < 1.5 & (interval.upper) > 1.5 ),
                bias = results.psa$ATE - 1.5,
                interval.width = 2*1.96*results.psa$se,model = "PSA")
xx$covered = (xx$interval.lower < 1.5 & xx$interval.upper > 1.5)

results = rbind(results,xx)
#results = xx
#save(results,file="/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 3/Code/psa_results.R")
save(results,file="/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 3/Code/interval_widths.R")

#####################################

#load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 3/Code/interval_widths.R")


#xx = aggregate(interval.width~data.type+method+model,results,mean)

#xx[order(xx$data.type,xx$method),]




################################################################
### code for table B1

library(ggplot2)
load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 2/Test code/scenario_1_sims_result.R")
#all.methods = results
#all.methods = subset(all.methods,data.type %in% c("confound","instru"))

all.methods$covered = as.character(all.methods$covered)
#ggplot(all.methods,aes(x = data.type,y=PROP,color=covered)) + geom_jitter() + facet_wrap(~method)

all.methods$covered = as.numeric(all.methods$covered)
all.methods$Var.ATE = all.methods$rubinse^2
all.methods$Var.between = all.methods$between.design.SE^2
all.methods$Var.within = all.methods$within.design.SE^2

xx = aggregate(ATE~data.type+method,all.methods,var)

xx = merge(xx,aggregate(Var.ATE ~ method+data.type,all.methods,mean),by=c('method','data.type',
                                                                            'data.type'))
xx = merge(xx,aggregate(Var.between ~ method+data.type,all.methods,mean),by=c('method','data.type',
                                                                          'data.type'))
xx = merge(xx,aggregate(Var.within ~ method+data.type,all.methods,mean),by=c('method','data.type',
                                                                          'data.type'))
xx = merge(xx,aggregate(bias ~ method+data.type,all.methods,mean),by=c('method','data.type',
                                                                       'data.type'))
xx = merge(xx,aggregate(PROP ~ method+data.type,all.methods,mean),by=c('method','data.type',
                                                                       'data.type'))
xx = merge(xx,aggregate(covered ~ method+data.type,all.methods,mean),by=c('method','data.type',
                                                                          'data.type'))
mse = function(x){mean((x-1.5)^2)}

xx = merge(xx,aggregate(ATE ~ method+data.type,all.methods,mse),by=c('method','data.type',
                                                                          'data.type'))

xx = merge(xx,aggregate(interval.width~method+data.type,all.methods,mean),by=c('method','data.type'))
#names(xx) = c('method','data.type','empirical.var','avg.total.var','avg.between.var','avg.within.var','bias','PROP','coverage')
names(xx) = c('Implementation','PS model','Empirical var',
              'Avg total var','Avg between var','Avg within var','Bias','PROP_DU','Coverage',
              'MSE','Interval width')

#xx

library(xtable)
xtable(xx,digits=3)

#######################################################
# Table B2
load("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 3/Code/psa_results.R")

results$Var.ATE = results$within.design.SE^2
bias = 1.5

xx = aggregate(ATE~method+data.type,results,var)

xx = merge(xx,aggregate(Var.ATE ~ method+data.type,results,mean),by=c('method','data.type'))

xx = merge(xx,aggregate(bias ~ method+data.type,results,mean),by=c('method','data.type'))
xx = merge(xx,aggregate(covered ~ method+data.type,results,mean),by=c('method','data.type'))
mse = function(x,bias){mean((x-bias)^2)}
xx = merge(xx,aggregate(ATE ~ method+data.type,results,mse,bias=bias),by=c('method','data.type'))
xx = merge(xx,aggregate(interval.width~method+data.type,results,mean),by=c('method','data.type'))

names(xx) = c('Implementation','PS model','Empirical var',
              'Avg total var','Bias','Coverage',
              'MSE','Interval width')

library(xtable)
xtable(xx,digits=3)

#####################

#interval widths figures (Figures B1-B3)

load("~/Dropbox/Shirley's Dissertation/Non-Github folder/BPSA paper/Draft 3/Code/interval_widths.R")

metric = "bias"
xx = eval(parse(text = paste("aggregate(", metric, "~data.type+method+model,results,mean)", sep="")))

#tiff("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/bias_plot.tiff")
png("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/bias_plot.png")
ggplot(xx,aes(x=data.type,y=bias,color=model))+
  geom_point()+facet_wrap(.~method)+geom_hline(yintercept=0)+
  labs(title="Bias under PSA/BPSA",x="Implementation",y="Bias") + 
  theme(axis.text.x = element_text(angle = 45))
dev.off()

metric = "covered"
xx = eval(parse(text = paste("aggregate(", metric, "~data.type+method+model,results,mean)", sep="")))

png("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/coverage_plot.png")
ggplot(xx,aes(x=data.type,y=covered,color=model))+
  geom_point()+facet_wrap(.~method)+geom_hline(yintercept=0.95)+
  labs(title="Coverage of 95% CI under PSA/BPSA",x="Implementation",y="Proportion coverage")+ 
  theme(axis.text.x = element_text(angle = 45))
dev.off()

metric = "mse"
results$mse = (results$bias)^2 + results$rubinse^2
xx = eval(parse(text = paste("aggregate(", metric, "~data.type+method+model,results,mean)", sep="")))

png("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/mse_plot.png")
ggplot(xx,aes(x=data.type,y=mse,color=model))+
  geom_point()+facet_wrap(.~method)+
  labs(title="MSE under PSA/BPSA",x="Implementation",y="MSE")+ 
  theme(axis.text.x = element_text(angle = 45))
dev.off()


######################################################
# Unused figures

xx$log.empirical_var = log(xx$empirical.var)
tiff("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/emp_var.tiff")
ggplot(xx,aes(x=PROP,y=log.empirical_var,color=method)) + geom_point()+ggtitle("Empirical variance")
dev.off()

all.methods$log_var.ATE = log(all.methods$Var.ATE)
tiff("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/total_var.tiff")
ggplot(all.methods,aes(x=PROP,y=log_var.ATE,color=method)) + geom_point()+ggtitle("Total variance")
dev.off()


tiff("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/coverage.tiff")
ggplot(xx,aes(x=PROP,y=coverage,color=method)) + geom_point()+ggtitle("Coverage")
dev.off()

###################################################
# Figure 1

all.methods$log_var.between = log(all.methods$Var.between)
tiff("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/between_var.tiff")
ggplot(all.methods,aes(x=PROP,y=log_var.between,color=method)) + geom_point()+ggtitle("Between design variance")
dev.off()

##########################################################
# Figure 2

all.methods$log_var.within = log(all.methods$Var.within)
tiff("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/within_var.tiff")
ggplot(all.methods,aes(x=PROP,y=log_var.within,color=method)) + geom_point()+ggtitle("Within design variance")
dev.off()
########################################################
# Figure 3

tiff("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/bias.tiff")
ggplot(all.methods,aes(x=PROP,y=bias,color=method)) + geom_point()+ggtitle("Bias")
dev.off()
#######################