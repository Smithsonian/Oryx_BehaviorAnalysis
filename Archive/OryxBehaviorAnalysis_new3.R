#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Project: SHO Stress and Behavior Analysis
# Date: 10 November 2016
# Author: Stephanie Cunningham
# Description: Plotting behavior means

# Editing Date: 06-Aug-2017
# Author: Jared Stabach
# Description: Going through code to understand what has been done to date.  
#   Some code cleanup. Need to incorporate random effects, since data on individual animals was repeatedly collected.
#   This could be conducted in frequentist perspective using a GLMM with individual as a random effect, or Bayesian

#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Clear the cache
rm(list=ls())

# Load necessary libraries
library(tidyr)

# Set working directory
setwd("C:/Users/stabachj/Dropbox (Smithsonian)/Projects/Oryx/StressAnalysis/Behavior")

# Read in file that was fixed in Excel
bdata <- read.csv("Behavior.Nov3.csv")

# Fix the data/time fields
bdata$TimeStart <- as.POSIXct(strptime(paste0(bdata$Date," ",bdata$TimeStart), format="%m/%d/%Y %H:%M"))
bdata$TimeEnd <- as.POSIXct(strptime(paste0(bdata$Date," ",bdata$TimeEnd), format="%m/%d/%Y %H:%M"))

# Re-order dataframe so all the behavior data is at the end
bdata <- bdata[,c(1:15,26:27,16:25)]

# I'm not sure that comparing the Treatment animals with the Control animals is the interesting question
# More, I think we want to know how behavior has changed across periods.  
# Essentially, each animal is its own control. 
# Our expectation is that adverse behaviors, such as head-shaking, should increase during the period animals are collared and then return to normal.
# The difficult piece is determing how long the after collaring period should be.

# Code the Control and Treatment records.
bdata$Control <- ifelse(bdata$Treatment == "control",1,2) 
bdata.control <- bdata[which(bdata$"Treatment" == "control"),]
bdata <- bdata[which(bdata$"Treatment" != "control"),]

# Comparison of behavior with Time factor (AdjObTime)
# This term has been coded for:
# 1: Prior to collar fitting
# 2: Immediately after collar fitting (0-3 days)
# 3: Greater than 3 days after collar fitting

# Set AdjObTime as a factor and look at Standing Head up (SHU)
bdata$AdjObTime <- as.factor(bdata$AdjObTime)
boxplot(pro.SHU~AdjObTime,data=bdata,boxwex=0.5,frame = FALSE,col=c("gray100","gray80","gray20"),main="Standing Head Up", xlab="Treatment Group", ylab="Percent of Activity") 
# The problem with this graph is that it does not take into account repeated measures.  

# Use lme4 to conduct linear regression analysis (glm), incorporating factor into model
# Must set the weight, since response variable is a proportion

library(lme4)
gm1 <- glmer(pro.SHU ~ AdjObTime + (1 | Animal), weights = TotalObs, 
             data = bdata, family = binomial)

#library(ggplot2)
#ggCaterpillar(ranef(gm1, postVar = TRUE))
#library(lattice)
#qqmath(ranef(gm1, postVar=TRUE)) 
#dotplot(ranef(gm1, condVar=TRUE))

# Warning message is expected.  We're tricking the model into not expected successes out of a total number of observations.

# Summarize model output
summary(gm1)

# Is this better than the null (no covariate)
model.null <- glmer(pro.SHU ~ 1 + (1 | Animal), weights = TotalObs, data=bdata, family = binomial)
summary(model.null)

anova(model.null, gm1)
# Yup...no surprise

# Now run multiple comparison test to identify where distances exist
library(multcomp)
posthoc <- glht(gm1, linfct = mcp(AdjObTime = "Tukey"))
summary(posthoc)

# For this variable, there is a decline in SHU in the period following collaring (both 3 days post collaring and > 3 days post collaring)
# Pre > 3-days = Post
# Observed when accounting for random effects

# Loop through all variabiles and place in a list
head(bdata)
str(bdata)
summary(bdata)

# The variables pro.walk and pro.oov both have NAs.  
# These need to be either removed or the variable cannon be included in analysis.
bdata$RSums <- rowSums(bdata[18:26], na.rm=TRUE)

summary(bdata$RSums)
# Some of the rows are < 1.  Set columns to 0.
bdata$pro.walk[is.na(bdata$pro.walk)] <- 0
bdata$pro.oov[is.na(bdata$pro.oov)] <- 0
summary(bdata)

# Create a null list to store the results
mod.list <- list()
anova.list <- list()
diff.list <- list()

# Reorder the dataframe
#bdata <- bdata[,c(1:19,21:27,20,28:29)]

# pro.lay does not converge....so removed
for(i in c(18:26)){
  print(paste0("Working On Variable: ",names(bdata[i])))
  name.Val <- names(bdata[i])
  gm1 <- glmer(bdata[,i] ~ AdjObTime + (1 | Animal), weights = TotalObs, 
               data = bdata, family = binomial)
  mod.list[[i-17]] <- summary(gm1)
  names(mod.list)[i-17] <- name.Val

  # Significant different than the null?
  anova.list[[i-17]] <- anova(model.null, gm1)
  names(anova.list)[i-17] <- name.Val
  
  # Conduct post-hoc tests to identify difference across treatments
  posthoc <- glht(gm1, linfct = mcp(AdjObTime = "Tukey"))
  diff.list[[i-17]] <-summary(posthoc)
  names(diff.list)[i-17] <- name.Val
}

# ***************************************
# ***************************************
# Model for Laying does not converge (Column 20)
# Model is unidentifiable.  Must check.

# But, rest of results are summarized in list.

# ***************************************
# ***************************************

# Better is to fit all models in a Bayesian framework so that can estimate the probabilities for each behavior
# And compared the posterior distributions
# Load library
library(R2jags)

# Set-up burn-in/iterations for JAGS
n.iter=10000 # Number of iterations
n.update=n.iter*0.20 # burn-in iterations (0.20 percent)
n.adapt=1000 # adaptation iterations
  
# Set up blank list
data.list <- vector("list")

bdata$Animal <- droplevels(bdata$Animal)
n.groups <- length(unique(bdata$Animal))

# Calculate the differenc in time from 7 am each day
bdata$Time <- difftime(bdata$TimeStart, paste0(strptime(bdata$TimeStart, format = "%Y-%m-%d")," ","7:00:00"),units="mins")
Time <-as.numeric(scale(bdata$Time))
bdata$Time <- Time

data.list=list(
  # There is both a TotalObs field and a ModTotObs field
  #C = round(bdata$TotalObs*bdata$pro.SHU), 
  C = bdata$ModTotObs*bdata$pro.SHU, 
  #N = bdata$TotalObs, 
  N = bdata$ModTotObs, 
  ID = as.numeric(bdata$Animal),
  PERIOD = bdata$AdjObTime,
  n.groups = length(unique(bdata$Animal)),
  n = length(bdata$Date),
  TimePast = bdata$Time
)

# Setup initial values
inits=list(
  list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1),mu.int=rnorm(1,0,1)), # Chain 1
  list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1),mu.int=rnorm(1,0,1)), # Chain 2
  list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1),mu.int=rnorm(1,0,1)) # Chain 3
)

# Fit model
jm2=jags.model("Binomial_GLMM.R",data=data.list,inits=inits,n.chains=length(inits),n.adapt=n.adapt)
update(jm2, n.iter=n.update) # Burn-in the chain
zm2=coda.samples(jm2,variable.names=c("alpha","beta1","beta2","beta3", "mu.int","sigma.int","Contrast1v2","Contrast1v3","Contrast2v3"), n.iter=n.iter, n.thin=1) # generate the coda object

#beta1 is the contrast with period 1 and 2
#beta2 is the contrast with period 1 and 3
#alphas are the random effects for each individual

# Deviance Information Criteria
zdic=dic.samples(jm2,n.iter=n.iter)
# Return DIC and deviaance
zdic

# Summarize object
print("*********************************************************************")
print(summary(zm2))

# Run convergence diagnostics
gelman.diag(zm2, multivariate=FALSE)
# Okay, values are 1.

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
Val.2.Plot <- c(16,1,2,3)
val.xlab <- c("Intercept","Beta 1","Beta 2","Beta 3")
val.header <- c("Intercept","Contrast - 1v2","Contrast - 1v3","Contrast - 2v3")

par(mfrow=c(length(Val.2.Plot),2))

for (i in 1:length(Val.2.Plot)){
  # Plot histogram, eliminating burn-in
  hist(df1[,Val.2.Plot[i]], freq=FALSE, breaks=100, xlim=c(min(df1[,Val.2.Plot[i]]),max(df1[,Val.2.Plot[i]])), main= paste0("Posterior Distribution of ",val.header[i]), xlab=val.xlab[i])
  # Overlay posterior distribution
  lines(density(df1[,Val.2.Plot[i]],adjust=3),col="black",lwd=2)
  lines(density(df2[,Val.2.Plot[i]],adjust=3),col="red",lwd=2)
  lines(density(df3[,Val.2.Plot[i]],adjust=3),col="blue",lwd=2)

  # Plot trace plot
  plot(df1[,Val.2.Plot[i]],xlab="Iteration Number",ylab=paste0("Value of "," ",val.header[i]),type="l", main="Trace Plot")
  lines(df2[,Val.2.Plot[i]],col="red")
  lines(df3[,Val.2.Plot[i]],col="blue")
  abline(a=mean(df1[,Val.2.Plot[i]]),b=0,col="green")
}

# **********************************************************
# **********************************************************
# Now loop through and run for all behaviors

# Create a null list to store the results
mod.summary <- list()
gelman.list <- list()
result.list <- list()

for(i in c(18:26)){
  print(paste0("Working On Variable: ",names(bdata[i])))
  name.Val <- names(bdata[i])
  
  data.list=list(
    #C = round(bdata$TotalObs*bdata[,i]), 
    #N = bdata$TotalObs, 
    C = bdata$ModTotObs*bdata[,i], 
    N = bdata$ModTotObs, 
    ID = as.numeric(bdata$Animal),
    PERIOD = bdata$AdjObTime,
    n.groups = length(unique(bdata$Animal)),
    n = length(bdata$Date),
    TimePast = bdata$Time
  )
  
  # Setup initial values
  inits=list(
    list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1),mu.int=rnorm(1,0,1)), # Chain 1
    list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1),mu.int=rnorm(1,0,1)), # Chain 2
    list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1),mu.int=rnorm(1,0,1)) # Chain 3
  )
  
  # Fit model
  jm2=jags.model("Binomial_GLMM.R",data=data.list,inits=inits,n.chains=length(inits),n.adapt=n.adapt)
  update(jm2, n.iter=n.update) # Burn-in the chain
  zm2=coda.samples(jm2,variable.names=c("alpha","beta1","beta2","beta3","mu.int","sigma.int","Contrast1v2","Contrast1v3","Contrast2v3"), n.iter=n.iter, n.thin=1) # generate the coda object
  
  # Deviance Information Criteria
  zdic=dic.samples(jm2,n.iter=n.iter)
  # Return DIC and deviaance
  zdic
  
  # Summarize object and put in a list
  mod.summary[[i-17]] <- summary(zm2)
  names(mod.summary)[i-17] <- name.Val
  
  # Run convergence diagnostics
  gelman.list[[i-17]] <- gelman.diag(zm2, multivariate=FALSE)
  names(gelman.list)[i-17] <- name.Val
  
  # create dataframes from coda objects
  # Plot the histograms and trace plots to examine and make sure that the parameter space has been explored
  df1 = as.data.frame(rbind(zm2[[1]]))
  df2 = as.data.frame(rbind(zm2[[2]]))
  df3 = as.data.frame(rbind(zm2[[3]]))
  
  # Remove the burn-in period
  df1 <- df1[(n.update+1):n.iter,]
  df2 <- df2[(n.update+1):n.iter,]
  df3 <- df3[(n.update+1):n.iter,]
  
  # Setup variables to plot and plotting window
  Val.2.Plot <- c(16,1,2,3)
  val.xlab <- c("Intercept","Beta 1","Beta 2","Beta 3")
  val.header <- c("Intercept","Contrast - 1v2","Contrast - 1v3","Contrast - 2v3")
  
  for (j in 1:length(Val.2.Plot)){
    # Plot histogram, eliminating burn-in
    hist(df1[,Val.2.Plot[j]], freq=FALSE, breaks=100, xlim=c(min(df1[,Val.2.Plot[j]]),max(df1[,Val.2.Plot[j]])), main= paste0("Posterior Distribution of ",val.header[j]), xlab=val.xlab[j])
    # Overlay posterior distribution
    lines(density(df1[,Val.2.Plot[j]],adjust=3),col="black",lwd=2)
    lines(density(df2[,Val.2.Plot[j]],adjust=3),col="red",lwd=2)
    lines(density(df3[,Val.2.Plot[j]],adjust=3),col="blue",lwd=2)
    
    # Plot trace plot
    plot(df1[,Val.2.Plot[j]],xlab="Iteration Number",ylab=paste0("Value of "," ",val.header[j]),type="l", main="Trace Plot")
    lines(df2[,Val.2.Plot[j]],col="red")
    lines(df3[,Val.2.Plot[j]],col="blue")
    abline(a=mean(df1[,Val.2.Plot[j]]),b=0,col="green")
  }
  
  # Bind together for future plotting/summary
  df1 <- rbind(df1,df2,df3)
  result.list[[i-17]] <- df1
  names(result.list)[i-17] <- name.Val
}

# **********************************************************
# **********************************************************
# Look at results
# **********************************************************
# **********************************************************

head(result.list$pro.SHU)
str(result.list)

# Look at convergence....all look good
gelman.list

# Can also look at the model summaries
mod.summary

# Bind the results together to plot together using the MCMCvis package
library(MCMCvis)

# Remove all of the columns except the first through contrast evaluations
# Then bind together in a dataframe
test <- lapply(result.list, function(x){x[,1:3]})
test.names <- names(test)
test <- do.call(cbind.data.frame, test)

# Convert to a matrix and plot
test <- as.matrix(test)

t1 <- rep(test.names,each=3)
t2 <- rep(c("- 1v2","- 1v3","- 2v3"),times=length(test.names))
t3 <- sprintf(paste0(t1,t2))
  
par(mfrow=c(1,1))
MCMCplot(test, labels=c(t3))

# **********************************************
MCMCplot(test, labels=c(t3), xlim=c(-5,5))

par(mfrow=c(2,2))
# Standing - Head up and Head down
stand <- test[,1:6]
tstand <- t3[1:6]
MCMCplot(stand, labels=c(tstand), xlim=c(-2,2))

# Feeding - Head up and Head down
feed <- test[,16:21]
tfeed <- t3[16:21]
MCMCplot(feed, labels=c(tfeed), xlim=c(-2,2))

# Headshake + scratch
head <- test[,c(10:12,22:24)]
thead <- t3[c(10:12,22:24)]
MCMCplot(head, labels=c(thead), xlim=c(-2,2))

# Walk + social
walk <- test[,c(13:15,25:27)]
twalk <- t3[c(13:15,25:27)]
MCMCplot(walk, labels=c(twalk), xlim=c(-2,2))

# Laying Down
lay <- test[,7:9]
tlay <- t3[7:9]
MCMCplot(lay, labels=c(tlay), xlim=c(-4,4))

# Replace the contrast between Treatment 2 and 3 with a dummy value for plotting
# This is a complete hack, but it help the plots look better
stand[,c(3,6)] <- -5
feed[,c(3,6)] <- -5
head[,c(3,6)] <- -5
walk[,c(3,6)] <- -5
lay[,3] <- -5

# Now move around
stand <- stand[,c(3,1:2,6,4:5)]
feed <- feed[,c(3,1:2,6,4:5)]
head <- head[,c(3,1:2,6,4:5)]
walk <- walk[,c(3,1:2,6,4:5)]
lay <- lay[,c(3,1:2)]

# Check
MCMCplot(stand, labels=c(tstand), xlim=c(-2,2))
MCMCplot(feed, labels=c(tfeed), xlim=c(-2,2))
MCMCplot(head, labels=c(thead), xlim=c(-2,2))
MCMCplot(walk, labels=c(twalk), xlim=c(-2,2))
MCMCplot(lay, labels=c(tlay), xlim=c(-4,4))

# *************************************************
# *************************************************

par(mfrow=c(2,2))

MCMCplot(head, labels=NULL,xlim=c(-2,2),labels_sz=1, 
         med_sz=1, thick_sz = 2.5, thin_sz = 1.25, ax_sz = 2.75, 
         x_axis_text_sz = 1, x_tick_text_sz = 1, mar = c(5,9,5,8), ref=0)
abline(h=3.5, lty=2, col=1)

mtext("(A)",side=2,at=7,adj=2,las=2,cex=1)

mtext("Headshake",side=2,at=6,adj=0.75,las=2,cex=1)
mtext("< Control",side=2,at=5,adj=0.5,las=2,cex=0.75)
mtext("< Control",side=2,at=4,adj=0.5,las=2,cex=0.75)

mtext("Scratch",side=2,at=3,adj=0.75,las=2,cex=1)
mtext("< Control",side=2,at=2,adj=0.5,las=2,cex=0.75)
mtext("< Control",side=2,at=1,adj=0.5,las=2,cex=0.75)

mtext("Trtmt 1 >",side=4,at=5,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 1 >",side=4,at=2,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 2 >",side=4,at=1,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 2 >",side=4,at=4,adj=2.25,las=2,cex=0.75)

# *************************************************

MCMCplot(stand, labels=NULL,xlim=c(-2,2),labels_sz=1, 
         med_sz=1, thick_sz = 2.5, thin_sz = 1.25, ax_sz = 2.75, 
         x_axis_text_sz = 1, x_tick_text_sz = 1, mar = c(5,9,5,8), ref=0)
abline(h=3.5, lty=2, col=1)

mtext("(B)",side=4,at=7,adj=1.5,las=2,cex=1)

mtext("Standing (HU)",side=2,at=6,adj=0.75,las=2,cex=1)
mtext("< Control",side=2,at=5,adj=0.5,las=2,cex=0.75)
mtext("< Control",side=2,at=4,adj=0.5,las=2,cex=0.75)

mtext("Standing (HD)",side=2,at=3,adj=0.75,las=2,cex=1)
mtext("< Control",side=2,at=2,adj=0.5,las=2,cex=0.75)
mtext("< Control",side=2,at=1,adj=0.5,las=2,cex=0.75)

mtext("Trtmt 1 >",side=4,at=5,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 1 >",side=4,at=2,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 2 >",side=4,at=1,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 2 >",side=4,at=4,adj=2.25,las=2,cex=0.75)

# *************************************************

MCMCplot(walk, labels=NULL,xlim=c(-2,2),labels_sz=1, 
         med_sz=1, thick_sz = 2.5, thin_sz = 1.25, ax_sz = 2.75, 
         x_axis_text_sz = 1, x_tick_text_sz = 1, mar = c(5,9,5,8), ref=0)
abline(h=3.5, lty=2, col=1)

mtext("(C)",side=2,at=7,adj=2,las=2,cex=1)

mtext("Walking",side=2,at=6,adj=0.75,las=2,cex=1)
mtext("< Control",side=2,at=5,adj=0.5,las=2,cex=0.75)
mtext("< Control",side=2,at=4,adj=0.5,las=2,cex=0.75)

mtext("Social",side=2,at=3,adj=0.75,las=2,cex=1)
mtext("< Control",side=2,at=2,adj=0.5,las=2,cex=0.75)
mtext("< Control",side=2,at=1,adj=0.5,las=2,cex=0.75)

mtext("Trtmt 1 >",side=4,at=5,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 1 >",side=4,at=2,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 2 >",side=4,at=1,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 2 >",side=4,at=4,adj=2.25,las=2,cex=0.75)

# *************************************************

MCMCplot(feed, labels=NULL,xlim=c(-2,2),labels_sz=1, 
         med_sz=1, thick_sz = 2.5, thin_sz = 1.25, ax_sz = 2.75, 
         x_axis_text_sz = 1, x_tick_text_sz = 1, mar = c(5,9,5,8), ref=0)
abline(h=3.5, lty=2, col=1)

mtext("(D)",side=4,at=7,adj=1.5,las=2,cex=1)

mtext("Feeding (HU)",side=2,at=6,adj=0.75,las=2,cex=1)
mtext("< Control",side=2,at=5,adj=0.5,las=2,cex=0.75)
mtext("< Control",side=2,at=4,adj=0.5,las=2,cex=0.75)

mtext("Feeding (HD)",side=2,at=3,adj=0.75,las=2,cex=1)
mtext("< Control",side=2,at=2,adj=0.5,las=2,cex=0.75)
mtext("< Control",side=2,at=1,adj=0.5,las=2,cex=0.75)

mtext("Trtmt 1 >",side=4,at=5,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 1 >",side=4,at=2,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 2 >",side=4,at=1,adj=2.25,las=2,cex=0.75)
mtext("Trtmt 2 >",side=4,at=4,adj=2.25,las=2,cex=0.75)

# *************************************************
# *************************************************

# Check into Principal Component Analysis
behavior.pca <- prcomp(bdata[,18:26])
behavior.pca

# Plot the standard deviation of the PCs
plot(behavior.pca, type="l")

# Summary PC importance
summary(behavior.pca)

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(behavior.pca, obs.scale = 1, var.scale = 1, 
              groups = bdata$AdjObTime, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

require(caret)
trans = preProcess(bdata[,18:19], 
                   method=c("BoxCox", "center", 
                            "scale", "pca"))
PC = predict(trans, bdata[,18:19])









# Could also do this with the MCMCvis library
library(MCMCvis)

# Plot all trace plots
MCMCtrace(zm2)

# Trace Plots all look good.
# Put the dfs together into 1
df1 <- rbind(df1,df2,df3)

# Now can start to address some cool items
# What is the probability that an animal is standing head up in an of the time periods
par(mfrow=c(1,1))
hist(plogis(df1$mu.int),freq=FALSE,breaks=100, main = "Standing Head-Up (SHU)", xlab="Probability", ylim=c(0,12), col="light grey")
lines(density(plogis(df1$mu.int)),lwd=2,col="black")
#plot(density(plogis(df1$mu.int)),lwd=2,col="black",ylim=c(0,12))
lines(density(plogis(df1$mu.int+df1$beta1)),lwd=2,col="blue") # This is the population in categorie 2
lines(density(plogis(df1$mu.int+df1$beta2)),lwd=2,col="red") # Categorie 3

# You could also calculate this probability using plogis and the model coefficient
Res.Df <- summary(zm2)
plogis(Res.Df$statistics[16])
quantile(plogis(df1$mu.int), c(0.025,0.975))

# Or, just grab them from the output table....might be slighly different because the table includes the burn-in period
plogis(Res.Df$quantiles[16,1])
plogis(Res.Df$quantiles[16,5])

# And what about in the secondary or tertiary periods
plogis(Res.Df$statistics[16,1]+Res.Df$statistics[14,1])
quantile(plogis(df1$mu.int+df1$beta1), c(0.025, 0.975))

plogis(Res.Df$statistics[16,1]+Res.Df$statistics[15,1])
quantile(plogis(df1$mu.int+df1$beta2), c(0.025, 0.975))

# You could also look at any of the individual animals
plogis(Res.Df$statistics[11,1])
# And in the second and third periods?
plogis(Res.Df$statistics[11,1]+Res.Df$statistics[14,1])
plogis(Res.Df$statistics[11,1]+Res.Df$statistics[15,1])

# Might also want to calculate to see if there is a difference between the time periods
# This is essentially a Tukey's test, like above:
hist(df1$beta2-df1$beta1,freq=FALSE,breaks=100, ylim=c(0,10),xlim=c(-1,0.5), main="Comparison Across Time Periods", xlab="Difference")
lines(density(df1$beta2-df1$beta1), lwd=2, col="black")
# This is the same as doing this...which I created within the model:
# lines(density(df1$Contrast2v3),lwd=2,col="red")
abline(v=0,col="red",lty=1,lwd=2)

# No difference between periods 2 and 3
# Is there a change from initial period to time 1?
#hist(df1$beta1)
lines(density(df1$beta1),lwd=2,col="blue")
#abline(v=0,col="red",lty=2,lwd=2)

# Is there a change from initial period to time 2?
# Yes, does not overlap 0
#hist(df1$beta2)
lines(density(df1$beta2),lwd=2,col="green")
#abline(v=0,col="red",lty=2,lwd=2)

# Then, using the MCMCvis packages, could also summarize
MCMCsummary(zm2)
# Same as summary(zm2)

# And could plot using a caterpillar plot
MCMCplot(zm2, params = "Contrast")

# Look at ggplot for creating caterpillar plot
# Dump results into table
# Append results for all activities together into 1 dataframe to graph together...showing differences.

#library(xtable)
#xtable(summary(df1[,1:3]))

#library(mcmcplots)
#caterplot(zm2,c("Contrast1v2","Contrast1v3","Contrast2v3"),collapse=TRUE,reorder=FALSE,las=3)

# Can convert the code mcmc.list to a matrix and plot
# This will allow multiple mcmc.list to be bounded together

# Plot the mcmc.list
MCMCplot(zm2)

# Or convert to a matrix and plot
test <- as.matrix(df1)
MCMCplot(test)

# See plot command for all the parameters to control the plot and what the open and closed circles mean
MCMCplot(test, med_sz=1)

test <- test[,1:3]
MCMCplot(test, med_sz=1)

test <- cbind(test,test)
MCMCplot(test, labels=c("what","what","what","what","what","what"))

# ***********************************************************************
# ***********************************************************************

# Set-up burn-in/iterations for JAGS
n.iter=10000 # Number of iterations
n.update=n.iter*0.20 # burn-in iterations (0.20 percent)
n.adapt=1000 # adaptation iterations

# Set up blank list
data.list <- vector("list")

# Reformat data
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

data.list=list(
  Y = y, 
  n.outcomes = ncol(y),
  #N = bdata$ModTotObs, 
  ID = as.numeric(bdata$Animal),
  PERIOD = bdata$AdjObTime,
  N = apply(y,1,sum),
  n.groups = length(unique(bdata$Animal)),
  n = nrow(y)
)

# Setup initial values
#inits=list(
 # list(beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1)), # Chain 1
 # list(beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1)), # Chain 2
 # list(beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1)) # Chain 3
#)

# Fit model
jm2=jags.model("Multinomial.R",data=data.list,n.chains=3,n.adapt=n.adapt)
update(jm2, n.iter=n.update) # Burn-in the chain
zm2=coda.samples(jm2,variable.names=c("alpha","beta"), n.iter=n.iter, n.thin=1) # generate the coda object

#beta1 is the contrast with period 1 and 2
#beta2 is the contrast with period 1 and 3
#alphas are the random effects for each individual

# Deviance Information Criteria
zdic=dic.samples(jm2,n.iter=n.iter)
# Return DIC and deviaance
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

test <- apply(df1[,1:9],2,mean)
# Probability of control activities
control.probs <- exp(test)/sum(exp(test))

test2 <- apply(df1[,1:9],2,mean) + apply(df1[,c(11,14,17,20,23,26,29,32,35)],2,mean)
per2.probs <- exp(test2)/sum(exp(test2))

seq.val <- seq(12,36,3)
test3 <- apply(df1[,1:9],2,mean) + apply(df1[,seq.val],2,mean)
per3.probs <- exp(test3)/sum(exp(test3))

# ***********************************************************************
# ***********************************************************************

# End Code

# ***********************************************************************
# ***********************************************************************