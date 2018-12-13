#*********************************************************************************************************
#*********************************************************************************************************

# Project: Scimitar-horned oryx Behavior Analysis
# Date: 12 December 2018
# Author: Grant Connette, Jared Stabach, and Stephanie Cunningham
# Contact: grmcco@gmail.com
# Description: Investigate behavioral changes in Scimitar-horned oryx fit with GPS collars. Data fit in a Bayesian framework to estimate the probability of each behavioral activity, based on a multinomial ilkelihood.
# Each animal used as their own control to assess how each behavior changed across time periods. Expectation is that adverse behaviors, such as headshaking, should increase during the period animals are collared (treatment) and then return to normal (post-treatment).

#*********************************************************************************************************
#*********************************************************************************************************

# Clear the cache
rm(list=ls())

# Load necessary libraries
library(tidyr)
library(reshape2)
library(ggplot2)
library(jagsUI)
library(MCMCvis)

# Read in file
bdata <- read.csv("bdata.csv", header=T, sep=",", row.names=1)

# View data
head(bdata)

# Set/Update the data/time fields
bdata$TimeStart <- as.POSIXct(bdata$TimeStart, format="%Y-%m-%d %H:%M")
bdata$TimeEnd <- as.POSIXct(bdata$TimeEnd, format="%Y-%m-%d %H:%M")

# Code the Control and Treatment records
# Remove the Control, too few to be useful
# Code the Control and Treatment records.
bdata$Control <- ifelse(bdata$Treatment == "control",1,2) 
bdata.control <- bdata[which(bdata$Treatment == "control"),]
bdata <- bdata[which(bdata$Treatment != "control"),]

# Set AdjObTime as a factor 
bdata$AdjObTime <- as.factor(bdata$AdjObTime)

# ***********************************************************************
# ***********************************************************************

# Set-up burn-in/iterations for JAGS
n.iter=500000 # Number of iterations
n.update=n.iter*0.20 # burn-in iterations (0.20 percent)

# Set up blank list
data.list <- vector("list")

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
  PERIOD = as.numeric(bdata$AdjObTime),
  N = apply(y,1,sum),
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
# *********************************
# *********************************

#save(jm2, file = "Behavior_Models.Rda")
load("Behavior_Models.Rda")

# Summarize object
print("*********************************************************************")
jm2

# *********************************
# *********************************
# Create the variance/covariance matrix (Fig. 3)
# *********************************
# *********************************

# Raw image
image(jm2$mean$sigma)

# Format to make it look nice
rho <- matrix(NA,6,6)
for (i in 1:6){
  for (j in 1:6){
    rho[i,j] <- jm2$mean$sigma[i,j]/sqrt(jm2$mean$sigma[i,i]*jm2$mean$sigma[j,j])  # correlation
    rho[i,j] <- round(rho[i,j],2)
  }
}

# Specify column names (i.e., behaviors)
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

# Reorder the correlation matrix
rho <- reorder_cormat(rho)
upper_tri <- get_upper_tri(rho)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Print Again with values and a few items cleaned up
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

# *********************************
# *********************************

# eps are the individual random effects
# tau.j's are the random effects for each behavior

# Graph the Probabilities
# *********************************
# *********************************

library(HDInterval)

# Investigate values in output
jm2$mean

# Look at trace and density plots to assess model convergence
MCMCtrace(jm2, params = 'PROBS', ind=TRUE, pdf=FALSE)

# Summarize values from the output, include median and highest posterior density intervals
Post.Summary <- MCMCsummary(jm2, 
                            params = 'PROBS',
                            Rhat = TRUE,
                            n.eff = TRUE,
                            func = function(x) c(median(x), hdi(x,credMass = 0.95)),
                            func_name = c('median','hdi_low','hdi_high'))

# View
Post.Summary

# Export file
write.csv(Post.Summary, "jm2_Output_Summary.csv")

# View the differences between the distributions
hist(jm2$sims.list$PROBS[,3,2]-jm2$sims.list$PROBS[,3,1], main = "Distribution Differences", xlab="Difference")

# Summarize Head-shaking
# What's the probability that Period 2 is greater than period 1.  Answer is 0.99 Probability
length(which(jm2$sims.list$PROBS[,2,4]>jm2$sims.list$PROBS[,1,4]))/length(jm2$sims.list$PROBS[,1,4])
# What about Period 3.  Never
length(which(jm2$sims.list$PROBS[,3,4]<jm2$sims.list$PROBS[,1,4]))/length(jm2$sims.list$PROBS[,1,4])
# Smaller.  Almost always (0.9999)

# Summarizes all, compares period 2 with period 1
apply(jm2$sims.list$PROBS[,c(1,2),],3,function(x) length(which((x[,2]-x[,1])>0))/75000)

# *************************************
# *************************************
# Plot the caterpillar plots from the MCMC output
MCMCplot(jm2, params = 'PROBS')

# Set all the Labels
main.label <- c("Head-Up", "Head-Down", "Laying", "Headshaking", "Locomotion", "Scratching")

# To plot all the parameters, while cleaning up the plot a bit
png('All_variables.png')
par(mfrow=c(2,3))

MCMCplot(jm2, params = c('PROBS\\[1,1\\]', 'PROBS\\[2,1\\]', 'PROBS\\[3,1\\]'), ref = Post.Summary[1,8], ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[1],
         med_sz=1, thin_sz = 1, thick_sz = 3, ax_sz=1, main_text_sz=1,
         labels=c('Pre-Trmt','Trmt','Post-Trmt'), xlab="Probability")
MCMCplot(jm2, params = c('PROBS\\[1,2\\]', 'PROBS\\[2,2\\]', 'PROBS\\[3,2\\]'), ref = Post.Summary[4,8], ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[2],
         med_sz=1, thin_sz = 1, thick_sz = 3, ax_sz=1, main_text_sz=1,
         labels=NULL, xlab="Probability")
MCMCplot(jm2, params = c('PROBS\\[1,3\\]', 'PROBS\\[2,3\\]', 'PROBS\\[3,3\\]'), ref = Post.Summary[7,8], ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[3],
         med_sz=1, thin_sz = 1, thick_sz = 3, ax_sz=1, main_text_sz=1,
         labels=NULL, xlab="Probability")

MCMCplot(jm2, params = c('PROBS\\[1,4\\]', 'PROBS\\[2,4\\]', 'PROBS\\[3,4\\]'), ref = Post.Summary[10,8], ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[4],
         med_sz=1, thin_sz = 1, thick_sz = 3, ax_sz=1, main_text_sz=1,
         labels=c('Pre-Trmt','Trmt','Post-Trmt'), xlab="Probability")
MCMCplot(jm2, params = c('PROBS\\[1,5\\]', 'PROBS\\[2,5\\]', 'PROBS\\[3,5\\]'), ref = Post.Summary[13,8], ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[5],
         med_sz=1, thin_sz = 1, thick_sz = 3, ax_sz=1, main_text_sz=1,
         labels=NULL, xlab="Probability")
MCMCplot(jm2, params = c('PROBS\\[1,6\\]', 'PROBS\\[2,6\\]', 'PROBS\\[3,6\\]'), ref = Post.Summary[16,8], ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[6],
         med_sz=1, thin_sz = 1, thick_sz = 3, ax_sz=1, main_text_sz=1,
         labels=NULL, xlab="Probability")
dev.off()

# And to only show the graphs where a change occurred
png(file = 'PROBS_variables.png',width=15, height=5, units = 'in', res=500)
layout(matrix(c(1,2,3), 1, 3, byrow = FALSE), widths=1, heights=1)

MCMCplot(jm2, params = c('PROBS\\[1,3\\]', 'PROBS\\[2,3\\]', 'PROBS\\[3,3\\]'), ref = Post.Summary[7,8], 
         ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[3],
         med_sz=1.5, thin_sz = 1, thick_sz = 3, ax_sz=2, main_text_sz=2,axis_text_sz=1.5,tick_text_sz = 1.5, 
         labels_sz = 2,
         labels=c('Pre-Trmt','Trmt','Post-Trmt'), xlab="Probability",
         mar = c(5.1, 6.1, 4.1, 2.1))
MCMCplot(jm2, params = c('PROBS\\[1,4\\]', 'PROBS\\[2,4\\]', 'PROBS\\[3,4\\]'), ref = Post.Summary[10,8], 
         ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[4],
         med_sz=1.5, thin_sz = 1, thick_sz = 3, ax_sz=2, main_text_sz=2,axis_text_sz=1.5,tick_text_sz = 1.5, 
         labels=NULL, xlab="Probability",
         mar = c(5.1, 6.1, 4.1, 2.1))
MCMCplot(jm2, params = c('PROBS\\[1,6\\]', 'PROBS\\[2,6\\]', 'PROBS\\[3,6\\]'), ref = Post.Summary[16,8], 
         ref_ovl = TRUE, ISB=FALSE, 
         main=main.label[6],
         med_sz=1.5, thin_sz = 1, thick_sz = 3, ax_sz=2, main_text_sz=2,axis_text_sz=1.5,tick_text_sz = 1.5, 
         labels=NULL, xlab="Probability")
dev.off()

# ***********************************************************************
# ***********************************************************************

# End Code

# ***********************************************************************
# ***********************************************************************