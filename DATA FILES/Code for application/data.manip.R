#call necessary fxns
library(readr)

#import .csv file
application_data_simulated_PM25 <- read_csv("~/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Data for application/application_data_simulated_PM25.csv")
data = application_data_simulated_PM25


#complete case analysis
data = data[complete.cases(data),]

#subset data to desired regions
data = subset(data,region %in% c("Northeast","Southeast","IndustrialMidwest"))

#change exposure to binary
data$TREAT = (data$Exposure_InMAP>4)
treatment = as.numeric(data$TREAT)
outcome = data$PM25

######choose covariates
cov.names = c("latitude","longitude","smokerate2000","TotPop",
              "PctRural","PctWhite","PctBlack","PctHighSchool","MedianHHInc",
              "PctPoor","PctFemale","PctOccupied","PctMovedIn5","MedianHValue","PopPerSQM",
              "avtmpf","avrelh")       

covariates = as.matrix(data[,names(data)%in%cov.names])

#load data into format which may be run by existing code
app.data = list()
app.data[[1]] = list()
app.data[[1]][["treatment"]] = treatment
app.data[[1]][["covariates"]] = data.frame(apply(covariates,2,as.numeric))
app.data[[1]][["outcome"]] = outcome
app.data[[2]] = 1

save(app.data,file="~/Dropbox/Shirley's Dissertation/Paper1/DATA FILES/Code for application/app.data.R")

