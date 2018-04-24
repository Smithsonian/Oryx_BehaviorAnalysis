#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Project: SHO Stress and Behavior Analysis
# Date: 10 November 2016
# Author: Stephanie Cunningham, Jared Stabach, Grant Connette
# Description: Summarize/investigate behavioral changes in Scimitar-horned oryx fit with GPS collars
#               Fit data in a Bayesian framework to estimate the probability of each behavioral activity
#               Data fit based on a multinomial likelihood
#               How does each behavior change across the time periods?  Using each animal as a control.
#               Expectation is that adverse behaviors, such as head-shaking, should increase during the period animals are collared 
#               and then return to normal.


#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Clear the cache
rm(list=ls())

# Load necessary libraries
library(tidyr)

# Set working directory
#setwd("C:/Users/stabachj/Dropbox (Smithsonian)/Projects/Oryx/StressAnalysis/Behavior")
setwd("C:/Users/Jared/Dropbox (Smithsonian)/Projects/Oryx/StressAnalysis/Behavior")

# Read in file
bdata <- read.csv("Behavior.Nov3.csv")

# Fix the data/time fields
bdata$TimeStart <- as.POSIXct(strptime(paste0(bdata$Date," ",bdata$TimeStart), format="%m/%d/%Y %H:%M"))
bdata$TimeEnd <- as.POSIXct(strptime(paste0(bdata$Date," ",bdata$TimeEnd), format="%m/%d/%Y %H:%M"))

# Re-order dataframe so all the behavior data is at the end
bdata <- bdata[,c(1:15,26:27,16:25)]

# Code the Control and Treatment records.
bdata$Control <- ifelse(bdata$Treatment == "control",1,2) 
bdata.control <- bdata[which(bdata$"Treatment" == "control"),]
bdata <- bdata[which(bdata$"Treatment" != "control"),]

# Look at data quickly
# Set AdjObTime as a factor 
# Summarize Standing Head up (SHU)
bdata$AdjObTime <- as.factor(bdata$AdjObTime)
boxplot(pro.SHU~AdjObTime,data=bdata,boxwex=0.5,frame = FALSE,col=c("gray100","gray80","gray20"),main="Standing Head Up", xlab="Treatment Group", ylab="Percent of Activity") 
# This does not account for repeated measures. 

# Variables pro.walk and pro.oov both have NAs.  
# Remove or the variable cannon be included in analysis.
bdata$RSums <- rowSums(bdata[18:26], na.rm=TRUE)

summary(bdata$RSums)
# Some of the rows are < 1.  Set columns to 0.
bdata$pro.walk[is.na(bdata$pro.walk)] <- 0
bdata$pro.oov[is.na(bdata$pro.oov)] <- 0
summary(bdata)

# ***********************************************************************
# ***********************************************************************

# Load library
library(R2jags)

# Set-up burn-in/iterations for JAGS
n.iter=10000 # Number of iterations
n.update=n.iter*0.20 # burn-in iterations (0.20 percent)
n.adapt=1000 # adaptation iterations

# Set up blank list
data.list <- vector("list")

# Reformat data and bind together
SHU <- as.integer(bdata$ModTotObs*bdata$pro.SHU)
SHD <- as.integer(bdata$ModTotObs*bdata$pro.SHD)
lay <- as.integer(bdata$ModTotObs*bdata$pro.lay)
HDSK <- as.integer(bdata$ModTotObs*bdata$pro.headshake)
WALK <- as.integer(bdata$ModTotObs*bdata$pro.walk)
FHU <- as.integer(bdata$ModTotObs*bdata$pro.FHU)
FHD <- as.integer(bdata$ModTotObs*bdata$pro.FHD)
SCRATCH <- as.integer(bdata$ModTotObs*bdata$pro.scratch)
SOCIAL <- as.integer(bdata$ModTotObs*bdata$pro.social)

y <- cbind(SHU,SHD,lay,HDSK,WALK,FHU,FHD,SCRATCH,SOCIAL) 
class(y)

# Setup the data list
data.list=list(
  Y = y, 
  n.outcomes = ncol(y),
  #ID = as.numeric(bdata$Animal),
  PERIOD = bdata$AdjObTime,
  N = apply(y,1,sum),
  #n.groups = length(unique(bdata$Animal)),
  n = nrow(y)
)

# Fit model
jm2=jags.model("Multinomial.R",data=data.list,n.chains=3,n.adapt=n.adapt)
update(jm2, n.iter=n.update) # Burn-in the chain
zm2=coda.samples(jm2,variable.names=c("alpha","beta","PROBS"), n.iter=n.iter, n.thin=1) # generate the coda object

# Deviance Information Criteria
zdic=dic.samples(jm2,n.iter=n.iter)
zdic

# Summarize object
print("*********************************************************************")
print(summary(zm2))

# Run convergence diagnostics
gelman.diag(zm2, multivariate=FALSE)

# Need to create the dataframes from the code objects
# Plot the histograms and trace plots to examine and make sure that the parameter space has been explored
df1 = as.data.frame(rbind(zm2[[1]]))
df2 = as.data.frame(rbind(zm2[[2]]))
df3 = as.data.frame(rbind(zm2[[3]]))

# Remove the burn-in period
df1 <- df1[(n.update+1):n.iter,]
df2 <- df2[(n.update+1):n.iter,]
df3 <- df3[(n.update+1):n.iter,]

# Setup variables to plot and plotting window
val.xlab <- colnames(y)

# Look at the trace plots for all the posterior distributions for the behaviors
par(mfrow=c(3,2))

for (i in 2:ncol(y)){
  # Plot histogram, eliminating burn-in
  hist(df1[,i], freq=FALSE, breaks=100, xlim=c(min(df1[,i]),max(df1[,i])), main= paste0("Posterior Distribution of ",val.xlab[i]), xlab=val.xlab[i])
  # Overlay posterior distribution
  lines(density(df1[,i],adjust=3),col="black",lwd=2)
  lines(density(df2[,i],adjust=3),col="red",lwd=2)
  lines(density(df3[,i],adjust=3),col="blue",lwd=2)
  
  # Plot trace plot
  plot(df1[,i],xlab="Iteration Number",ylab=paste0("Value of "," ",val.xlab[i]),type="l", main="Trace Plot")
  lines(df2[,i],col="red")
  lines(df3[,i],col="blue")
  abline(a=mean(df1[,i]),b=0,col="green")
}

# Summarize the activity coefficients....am just using chain 1 here (df1)
coefs.bhv <- apply(df1[,1:9],2,mean)
quant.bhv <- apply(df1[,1:9],2,quantile)

# Probability of control activities
control.probs <- exp(coefs.bhv)/sum(exp(coefs.bhv))

seq.val1 <- seq(11,35,3)
coefs.time2 <- apply(df1[,1:9],2,mean) + apply(df1[,seq.val1],2,mean)
per2.probs <- exp(coefs.time2)/sum(exp(coefs.time2))

seq.val2 <- seq(12,36,3)
coefs.time3 <- apply(df1[,1:9],2,mean) + apply(df1[,seq.val2],2,mean)
per3.probs <- exp(coefs.time3)/sum(exp(coefs.time3))

# ***********************************************************************
# ***********************************************************************

testing1 <- df1[,9] + df1[,34]
testing2 <- df1[,35]
testing3 <- df1[,36]

test <- cbind(testing2,testing3)

test <- as.matrix(test)

par(mfrow=c(1,1))
MCMCplot(test, labels=c("Treatment1","Treatment2"),xlim=c(-1.5,1.5))

# To graph the differences, need to append to a dataframe
test <- exp(df1[,1])/sum(exp(df1[,1:9]))
test <- df1[,1] - df1[,12]

# End Code

# ***********************************************************************
# ***********************************************************************