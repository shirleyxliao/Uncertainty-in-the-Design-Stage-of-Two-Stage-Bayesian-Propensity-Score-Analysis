################### final simulation graphs
# function to restructure data

data.fxn = function(data){
  ncolumns = dim(data)[2]
  
  for(c in 1:ncolumns){
    if(c==1){
      returnframe = data.frame(data[,c])
    } else{
      returnframe = rbind(returnframe,data.frame(data[,c]))
    }
  }
  returnframe
}

#################################
# read in results and restructure them

load("~/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Results/BPSA.confound.R")
all.methods = data.frame(data.type = "confound",data.fxn(bpsa.result))

load("~/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Results/BPSA.instru.R")
all.methods = rbind(all.methods,data.frame(data.type = "instru",data.fxn(bpsa.result)))

load("~/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Results/BPSA.instru_prog.R")
all.methods = rbind(all.methods,data.frame(data.type = "instru_prog",data.fxn(bpsa.result)))

load("~/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Results/BPSA.instru_prog_noise.R")
all.methods = rbind(all.methods,data.frame(data.type = "instru_prog_noise",data.fxn(bpsa.result)))

all.methods$sandwich_se = sqrt(all.methods$between.design.SE^2 + all.methods$within.design.sandwich.SE^2)

all.methods$PROP = (all.methods$between.design.SE)^2/(all.methods$rubinse)^2
all.methods$PROP_sandwich = (all.methods$between.design.SE)^2/(all.methods$sandwich_se)^2

all.methods$lower = all.methods$ATE - 1.96*all.methods$rubinse
all.methods$upper = all.methods$ATE + 1.96*all.methods$rubinse

all.methods$lower_sandwich = all.methods$ATE - 1.96*all.methods$sandwich_se
all.methods$upper_sandwich = all.methods$ATE + 1.96*all.methods$sandwich_se

all.methods$covered = as.numeric(1.5<all.methods$upper & 1.5>all.methods$lower)
all.methods$covered_sandwich = as.numeric(1.5<all.methods$upper_sandwich & 1.5>all.methods$lower_sandwich)
all.methods$bias = all.methods$ATE - 1.5

save(all.methods,file="~/Dropbox/Shirley\'s Dissertation/Paper1/DATA FILES/Test code/scenario_1_sims_result.R")

################################################################
### code for table B1

library(ggplot2)
load("~/Dropbox/Shirley\'s Dissertation/Paper1/DATA FILES/Test code/scenario_1_sims_result.R")

#all.methods = subset(all.methods,data.type %in% c("confound","instru"))

all.methods$covered = as.character(all.methods$covered)
ggplot(all.methods,aes(x = data.type,y=PROP,color=covered)) + geom_jitter() + facet_wrap(~method)

all.methods$covered = as.numeric(all.methods$covered)
all.methods$Var.ATE = all.methods$rubinse^2
all.methods$Var.between = all.methods$between.design.SE^2
all.methods$Var.within = all.methods$within.design.SE^2

xx = aggregate(ATE~data.type+method,all.methods,var)

xx = merge(xx,aggregate(Var.ATE ~ method+data.type,all.methods,mean),by=c('method',
                                                                            'data.type'))
xx = merge(xx,aggregate(Var.between ~ method+data.type,all.methods,mean),by=c('method',
                                                                          'data.type'))
xx = merge(xx,aggregate(Var.within ~ method+data.type,all.methods,mean),by=c('method',
                                                                          'data.type'))
xx = merge(xx,aggregate(bias ~ method+data.type,all.methods,mean),by=c('method',
                                                                       'data.type'))
xx = merge(xx,aggregate(PROP ~ method+data.type,all.methods,mean),by=c('method',
                                                                       'data.type'))
xx = merge(xx,aggregate(covered ~ method+data.type,all.methods,mean),by=c('method',
                                                                          'data.type'))
mse = function(x){mean((x-1.5)^2)}

xx = merge(xx,aggregate(ATE ~ method+data.type,all.methods,mse),by=c('method',
                                                                          'data.type'))
#names(xx) = c('method','data.type','empirical.var','avg.total.var','avg.between.var','avg.within.var','bias','PROP','coverage')
names(xx) = c('Implementation','PS model','Empirical var',
              'Avg total var','Avg between var','Avg within var','Bias','PROP_DU','Coverage',
              'MSE')

xx

library(xtable)
xtable(xx,digits=3)

########################################################
# Figure 3

tiff("/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper1/Figures/bias.tiff")
ggplot(all.methods,aes(x=PROP,y=bias,color=method)) + geom_point()+ggtitle("Bias")
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

#######################