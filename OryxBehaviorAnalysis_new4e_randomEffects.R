#*********************************************************************************************************
#*********************************************************************************************************

# Project: SHO Stress and Behavior Analysis
# Date: 10 November 2016
# Author: Stephanie Cunningham, Jared Stabach, Grant Connette
# Description: Summarize/investigate behavioral changes in Scimitar-horned oryx fit with GPS collars
#               Fit data in a Bayesian framework to estimate the probability of each behavioral activity
#               Data fit based on a multinomial likelihood
#               How does each behavior change across the time periods?  Using each animal as a control.
#               Expectation is that adverse behaviors, such as head-shaking, should increase during the period animals are collared and then return to normal.
# New comment

#*********************************************************************************************************
#*********************************************************************************************************

# Clear the cache
rm(list=ls())

# Load necessary libraries
library(tidyr)

# Read in file
bdata <- read.csv("Behavior.Nov3.csv")

# Let's combine the proportions for HU (SHU + FHU: Standing/Feeding Head-Up) and HD (SHD + FHD: Standing/Feeding Head-Down).
bdata$pro.HU <- bdata$pro.SHU + bdata$pro.FHU
bdata$pro.HD <- bdata$pro.SHD + bdata$pro.FHD

bdata$pro.walk[is.na(bdata$pro.walk)] <- 0
bdata$pro.Loco <- bdata$pro.walk + bdata$pro.social

# Delete and reorganize
bdata <- bdata[,-c(16,17,20,21,22,24)]
# Reorganize
bdata <- bdata[,c(1:15,22:24,16:21)]

# Fix the data/time fields
bdata$TimeStart <- as.POSIXct(strptime(paste0(bdata$Date," ",bdata$TimeStart), format="%m/%d/%Y %H:%M"))
bdata$TimeEnd <- as.POSIXct(strptime(paste0(bdata$Date," ",bdata$TimeEnd), format="%m/%d/%Y %H:%M"))

# Re-order dataframe so all the behavior data is at the end
#bdata <- bdata[,c(1:15,26:27,16:25)]
bdata <- bdata[,c(1:15,23:24,16:22)]

# Remove any rows where variables of interest (columns 18:24) have NAs
# These are all the oov (Out of View) records, so doesn't really matter, since not included in analysis.
# Set to Zero
bdata[18:24][is.na(bdata[18:24])] <- 0

bdata$RSums <- rowSums(bdata[18:24])
summary(bdata$RSums)

# Reformat data and bind together
bdata$HU <- as.integer(bdata$ModTotObs*bdata$pro.HU)
bdata$HD <- as.integer(bdata$ModTotObs*bdata$pro.HD)
bdata$LAY <- as.integer(bdata$ModTotObs*bdata$pro.lay)
bdata$HDSK <- as.integer(bdata$ModTotObs*bdata$pro.headshake)
#WALK <- as.integer(bdata$ModTotObs*bdata$pro.walk)
bdata$LOCO <- as.integer(bdata$ModTotObs*bdata$pro.Loco)
#FHU <- as.integer(bdata$ModTotObs*bdata$pro.FHU)
#FHD <- as.integer(bdata$ModTotObs*bdata$pro.FHD)
bdata$SCRATCH <- as.integer(bdata$ModTotObs*bdata$pro.scratch)
bdata$OOV <- as.integer(bdata$ModTotObs*bdata$pro.oov)
#SOCIAL <- as.integer(bdata$ModTotObs*bdata$pro.social)

# Re-Order
bdata <- bdata[,c(1:17,26:32,25)]

# ************************************
# ************************************

# Write to a new file.  This will be the file that is shared
write.csv(bdata, file = "bdata.csv")

# ***********************************************************************
# ***********************************************************************

# Read in file
# Some of the factors (Feeding Head-Up + Standing Head-Up have already been collapsed (Head-Up))
bdata <- read.csv("bdata.csv", header=T, sep=",", row.names=1)

# View data
head(bdata)

# Set/Update the data/time fields
bdata$TimeStart <- as.POSIXct(bdata$TimeStart, format="%Y-%m-%d %H:%M")
bdata$TimeEnd <- as.POSIXct(bdata$TimeEnd, format="%Y-%m-%d %H:%M")

# Code the Control and Treatment records
# Remove the Control, too few to be useful
bdata$Control <- ifelse(bdata$Treatment == "control",1,2) 
bdata.control <- bdata[which(bdata$Treatment == "control"),]
bdata <- bdata[which(bdata$Treatment != "control"),]

# Set AdjObTime as a factor 
bdata$AdjObTime <- as.factor(bdata$AdjObTime)

# ***********************************************************************
# ***********************************************************************

# Load library
library(jagsUI)

# Set-up burn-in/iterations for JAGS
n.iter=500000 # Number of iterations
n.update=n.iter*0.20 # burn-in iterations (0.20 percent)
#n.adapt=1000 # adaptation iterations

# Set up blank list
data.list <- vector("list")

# Reformat data and bind together
##HU <- as.integer(bdata$ModTotObs*bdata$pro.HU)
##HD <- as.integer(bdata$ModTotObs*bdata$pro.HD)
##LAY <- as.integer(bdata$ModTotObs*bdata$pro.lay)
##HDSK <- as.integer(bdata$ModTotObs*bdata$pro.headshake)
#WALK <- as.integer(bdata$ModTotObs*bdata$pro.walk)
##LOCO <- as.integer(bdata$ModTotObs*bdata$pro.Loco)
#FHU <- as.integer(bdata$ModTotObs*bdata$pro.FHU)
#FHD <- as.integer(bdata$ModTotObs*bdata$pro.FHD)
##SCRATCH <- as.integer(bdata$ModTotObs*bdata$pro.scratch)
#SOCIAL <- as.integer(bdata$ModTotObs*bdata$pro.social)

y <- cbind(bdata$HU,bdata$HD,bdata$LAY,bdata$HDSK,bdata$LOCO,bdata$SCRATCH) 
class(y)

# Create matrix for inverse Wishart prior on individual random effects
R <- matrix(0,nrow=6,ncol=6)
for (i in 1:6){
  R[i,i] <- 0.1
}

# Setup the data list
data.list=list(
  Y = y, 
  n.outcomes = ncol(y),
  #ID = as.numeric(bdata$Animal),
  PERIOD = as.numeric(bdata$AdjObTime),
  N = apply(y,1,sum),
  #n.groups = length(unique(bdata$Animal)),
  n = nrow(y),
  ind = as.numeric(droplevels(bdata$Animal)),
  nind = length(unique(bdata$Animal)),
  R = R
)

# Fit model
jm2=jags(model.file = "Multinomial_withREs.R",
         data=data.list,
         n.chains=3,n.iter=n.iter,n.thin=20,parallel = F,
         parameters.to.save = c("alpha","beta","sigma","PROBS","eps"))

# Save the jags model
#save(jm2, file = "Behavior_Models.Rda")
load("Behavior_Models.Rda")

# Summarize object
print("*********************************************************************")
jm2

# Look at the variance/covariance matrix
image(jm2$mean$sigma)

# Converting variance/covariance to correlation matrix
rho <- matrix(NA,6,6)
for (i in 1:6){
  for (j in 1:6){
    rho[i,j] <- jm2$mean$sigma[i,j]/sqrt(jm2$mean$sigma[i,i]*jm2$mean$sigma[j,j])  # correlation
    rho[i,j] <- round(rho[i,j],2)
  }
}

# Specify column names
rownames(rho) <- colnames(rho) <- c("HU","HD","LAY","HDSK","LOCO","SCRATCH")

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(rho)
upper_tri

# Melt the correlation matrix
library(reshape2)
library(ggplot2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
names(melted_cormat) <- c("Behavior1", "Behavior2", "value")

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
rho <- reorder_cormat(rho)
upper_tri <- get_upper_tri(rho)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

# Print the heatmap
print(ggheatmap)

ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.4, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

# eps are the individual random effects
# tau.j's are the random effects for each behavior
# Likely want to graph the probabilities, or at least those that are interesting.
# Then have a table of all the probabilities, inclusive of credible intervals

# Use the HDInterval package, to more appropriately calculate the credible intervals on posterior distributions that might be skewed
library(HDInterval)
# Interesting to view the alpha's and the tau.j random effects (Precision).  
# These can be converted by divided by the inverse to get the variances.  
# Will provide strength of coefficients estimates, along with variance around estimates.

# Graph the Probabilities
#


head(jm2$sims.list$PROBS)
ncol(jm2$sims.list$PROBS)

# Look at traceplots to determine if proper exploration of the parameter space has occurred.
par(ask=FALSE)
traceplot(jm2)
traceplot(jm2, parameters = 'PROBS[1:3]')

par(mfrow=c(2,1))

whiskerplot(jm2, parameters='PROBS')

jm2$mean
jm2$summary

library(MCMCvis)
# Look at trace and density plots to assess model convergence
MCMCtrace(jm2, params = 'PROBS', ind=TRUE, pdf=FALSE)
temp <- rep('First Param', 18)

MCMCplot(jm2, params = 'PROBS', ref=NULL, xlab='ESTIMATE', main = "MCMCvis plot", labels = temp)

ex2 <- MCMCchains(jm2, params = 'alpha', excl = 'alpha\\[1\\]', mcmc.list = TRUE, ISB = FALSE)
ex2 <- MCMCchains(jm2, params = c('alpha\\[2\\]','alpha\\[3\\]','alpha\\[4\\]'), mcmc.list = TRUE, ISB = FALSE)
MCMCplot(ex2)
ex2 <- MCMCchains(jm2, params = 'beta\\[1\\]', mcmc.list = TRUE, ISB = FALSE)


ex2 <- MCMCchains(jm2, params = 'PROBS\\[1,1\\]', mcmc.list = TRUE, ISB = FALSE)
MCMCplot(jm2, params = c('alpha\\[2\\]', 'alpha\\[4\\]', 'alpha\\[3\\]'), ISB=FALSE, main="Cntl v Trmt 2",med_sz=0, thin_sz = 1, thick_sz = 3, ax_sz=1, main_text_sz=1)


ex2 <- MCMCchains(jm2, params = 'alpha\\[2\\]', mcmc.list = TRUE, ISB = FALSE)
MCMCplot(jm2, params = c('alpha\\[2\\]', 'alpha\\[4\\]', 'alpha\\[3\\]'), ISB=FALSE, labels = c('whatever 1','whatever 2','whatever 3'))
         

MCMCplot(jm2, params = c('PROBS\\[1,6\\]', 'PROBS\\[2,6\\]', 'PROBS\\[3,6\\]'), ISB=FALSE, ref_ovl = FALSE)

MCMCplot(jm2, params = 'PROBS')
MCMCsummary(ex2)

MCMCplot(jm2, params = c('PROBS\\[1,1\\]', 'PROBS\\[2,1\\]', 'PROBS\\[3,1\\]'), ISB=FALSE, ref_ovl = FALSE)


par(mfrow=c(1,2))




# Setup variables to plot and plotting window
val.xlab <- colnames(y)

# Look at the trace plots for all the posterior distributions for the behaviors
par(mfrow=c(3,2))

# Plot the control probabilities for each behavior to investigate proper exploration of the parameter space
val.2.plot <- 22:28
#plot.seq <- seq(1,25,3)
#plot.seq.Trmt1 <- seq(2,26,3)
#plot.seq.Trmt2 <- seq(3,27,3)

for (i in 2:length(val.xlab)){
  # Plot histogram, eliminating burn-in
  #print(i)
  hist(df1[,val.2.plot[i]], freq=FALSE, breaks=100, xlim=c(min(df1[,val.2.plot[i]]),max(df1[,val.2.plot[i]])), main= paste0("Posterior Distribution of ",val.xlab[i]), xlab=val.xlab[i])
  # Overlay posterior distribution
  lines(density(df1[,val.2.plot[i]],adjust=3),col="black",lwd=2)
  lines(density(df2[,val.2.plot[i]],adjust=3),col="red",lwd=2)
  lines(density(df3[,val.2.plot[i]],adjust=3),col="blue",lwd=2)
  
  # Plot trace plot
  plot(df1[,val.2.plot[i]],xlab="Iteration Number",ylab=paste0("Value of "," ",val.xlab[i]),type="l", main="Trace Plot")
  lines(df2[,val.2.plot[i]],col="red")
  lines(df3[,val.2.plot[i]],col="blue")
  abline(a=mean(df1[,val.2.plot[i]]),b=0,col="green")
}

# Summarize the activity probabilities....am just using chain 1 here (df1)
# plot.seq... are indexes to calculate the probability summaries for each of the treatment groups, including the control.
# The summed probabilities should all sum to 1 across each activity.  This was the whole reason for going to moving to a multinomial regression
coefs.bhv <- apply(df1[,val.2.plot],2,mean)
quant.bhv <- apply(df1[,val.2.plot],2,quantile)

# Probability of control activities
(control.probs <- exp(coefs.bhv)/sum(exp(coefs.bhv)))
sum(control.probs)
# Which should be the same as
control.seq <- seq(1,19,3)
(coefs.bhv.test <- apply(df1[,control.seq],2,mean))

# Now do the same for Trmt 1
#seq.val1 <- seq(11,35,3)
plot.seq.Trmt1 <- seq(30,49,3)
coefs.time2 <- apply(df1[,22:28],2,mean) + apply(df1[,plot.seq.Trmt1],2,mean)
(per2.probs <- exp(coefs.time2)/sum(exp(coefs.time2)))

# Which should be the same as
control.seq1 <- seq(2,20,3)
(coefs.bhv.test <- apply(df1[,control.seq1],2,mean))

#seq.val2 <- seq(12,36,3)
plot.seq.Trmt2 <- seq(31,49,3)
coefs.time3 <- apply(df1[,22:28],2,mean) + apply(df1[,plot.seq.Trmt2],2,mean)
(per3.probs <- exp(coefs.time3)/sum(exp(coefs.time3)))

# Which should be the same as
control.seq <- seq(3,21,3)
(coefs.bhv.test <- apply(df1[,control.seq],2,mean))

# Look at the probabilities for each time period
colnames(y)
control.probs
per2.probs
per3.probs

# Probabilities should add up to 1
sum(control.probs)
sum(per2.probs)
sum(per3.probs)

# ***********************************************************************
# ***********************************************************************

# Separate out the probabilities
# From this, could graph the probability of doing each activity or across each treatment.
df.prob <- df1[,1:21]

# Separate the alpha and beta coefficients, to compare the effects
df.test <- df1[,-1:-21]
testing1 <- df.test[,7] + df.test[,26]
testing2 <- df.test[,27]
testing3 <- df.test[,28]

test <- cbind(testing2,testing3)

test <- as.matrix(test)

library(MCMCvis)

par(mfrow=c(1,1))
MCMCplot(test, labels=c("Treatment1","Treatment2"),xlim=c(-5,5))

# Loop through all the behaviors, creating a graph for each
Trt1 <- seq(9,28,3)
Trt2 <- seq(10,28,3)
val.xlab

par(mfrow=c(1,2))

# Extract the values to plot
testing2 <- df.test[,Trt1]
testing3 <- df.test[,Trt2]

# Convert to a matrix
testing2 <- as.matrix(testing2)
testing3 <- as.matrix(testing3)

# Plot the results
MCMCplot(testing2, labels=val.xlab[1:7],xlim=c(-5,5),main="Cntl v Trmt 1", med_sz=0, thin_sz = 1, thick_sz = 3, ax_sz=1, x_axis_text_sz=1, x_tick_text_sz=1, main_text_sz=1)
MCMCplot(testing3, labels = val.xlab[1:7], xlim=c(-5,5),main="Cntl v Trmt 2",med_sz=0, thin_sz = 1, thick_sz = 3, ax_sz=1, x_axis_text_sz=1, x_tick_text_sz=1, main_text_sz=1)




# To graph the differences, need to append to a dataframe
test <- exp(df.test[,1])/sum(exp(df.test[,1:9]))
test <- df.test[,1] - df.test[,12]

# ***********************************************************************
# ***********************************************************************

# End Code

# ***********************************************************************
# ***********************************************************************
