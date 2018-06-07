#####run file
foldername = "DissertationJune"

#load run fxns
source(paste('/n/home10/silverfoxflower/',foldername,'/PSAfxn.R',sep=""))
source(paste('/n/home10/silverfoxflower/',foldername,'/BPSAfxn.R',sep=""))

#run application data
load(paste('/n/home10/silverfoxflower/',foldername,'/app.data.R',sep=""))
PSA.fxn(app.data,savefile='app.PSA.large',foldername)
BPSA.fxn(app.data,savefile='app.BPSA.large',foldername)
