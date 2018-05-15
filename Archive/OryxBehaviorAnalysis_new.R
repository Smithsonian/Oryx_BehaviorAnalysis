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
#setwd("C:/Users/stabachj/Dropbox (Smithsonian)/Projects/Oryx/StressAnalysis/Behavior")
setwd("C:/Users/Jared/Dropbox (Smithsonian)/Projects/Oryx/StressAnalysis/Behavior")

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
# The variables pro.walk and pro.oov both have NAs.  These need to be either removed or the variable cannon be included in analysis.
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

bdata <- bdata[,c(1:19,21:27,20,28:29)]

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

data.list=list(
  C = round(bdata$TotalObs*bdata$pro.SHU), 
  N = bdata$TotalObs, 
  ID = as.numeric(bdata$Animal),
  PERIOD = bdata$AdjObTime,
  n.groups = length(unique(bdata$Animal)),
  n = length(bdata$Date)
)

# Setup initial values
inits=list(
  list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),mu.int=rnorm(1,0,1)), # Chain 1
  list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),mu.int=rnorm(1,0,1)), # Chain 2
  list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),mu.int=rnorm(1,0,1)) # Chain 3
)

# Fit model
jm2=jags.model("Binomial_GLMM.R",data=data.list,inits=inits,n.chains=length(inits),n.adapt=n.adapt)
update(jm2, n.iter=n.update) # Burn-in the chain
zm2=coda.samples(jm2,variable.names=c("alpha","beta1","beta2","mu.int","sigma.int","Contrast1v2","Contrast1v3","Contrast2v3"), n.iter=n.iter, n.thin=1) # generate the coda object

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
    C = round(bdata$TotalObs*bdata[,i]), 
    N = bdata$TotalObs, 
    ID = as.numeric(bdata$Animal),
    PERIOD = bdata$AdjObTime,
    n.groups = length(unique(bdata$Animal)),
    n = length(bdata$Date)
  )
  
  # Setup initial values
  inits=list(
    list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),mu.int=rnorm(1,0,1)), # Chain 1
    list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),mu.int=rnorm(1,0,1)), # Chain 2
    list(alpha=rnorm(n.groups,0,2),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),mu.int=rnorm(1,0,1)) # Chain 3
  )
  
  # Fit model
  jm2=jags.model("Binomial_GLMM.R",data=data.list,inits=inits,n.chains=length(inits),n.adapt=n.adapt)
  update(jm2, n.iter=n.update) # Burn-in the chain
  zm2=coda.samples(jm2,variable.names=c("alpha","beta1","beta2","mu.int","sigma.int","Contrast1v2","Contrast1v3","Contrast2v3"), n.iter=n.iter, n.thin=1) # generate the coda object
  
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








# See plot command for all the parameters to control the plot and what the open and closed circles mean
MCMCplot(test, med_sz=1)

test <- test[,1:3]
MCMCplot(test, med_sz=1)

test <- cbind(test,test)
MCMCplot(test, labels=c("what","what","what","what","what","what"))









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


# ***************************************
# ***************************************.

# Stephanie's code from here, edited slightly.

# ***************************************
# ***************************************

# Collared Females
f.clr1 <- subset(bdata, Sex=="Female" & (Treatment != "control") & AdjObTime==1)
f.clr2 <- subset(bdata, Sex=="Female" & (Treatment != "control") & AdjObTime==2)
f.clr3 <- subset(bdata, Sex=="Female" & (Treatment != "control") & AdjObTime==3)

# Summarize the important behavior variables 
fclr.mns1 <- colMeans(f.clr1[18:26],na.rm=TRUE)
fclr.mns2 <- colMeans(f.clr2[18:26],na.rm=TRUE)
fclr.mns3 <- colMeans(f.clr3[18:26],na.rm=TRUE)
fclr.mns <- rbind(fclr.mns1,fclr.mns2,fclr.mns3)
fclr.mns <- as.data.frame(fclr.mns, row.names = FALSE)

# Assign Column headers
col.names <- c("SHU","SHD","Laying","Headshaking","Walking","FHU","FHD","Scratch","Social")
colnames(fclr.mns) <- col.names

# Do same for sd's
# *************************************************************
fcl.sd1 <- apply(f.clr1[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
fcl.sd2 <- apply(f.clr2[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
fcl.sd3 <- apply(f.clr3[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
fcl.sd <- rbind(fcl.sd1,fcl.sd2,fcl.sd3)
fcl.sd <- as.data.frame(fcl.sd, row.names = FALSE)

# Column headers
colnames(fcl.sd) <- col.names

# Add in character field and re-order
fclr.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fclr.mns <- fclr.mns[,c(10,1:9)] 
fcl.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fcl.sd <- fcl.sd[,c(10,1:9)]

# Use gather command from tidyr to re-organize the variables
fcl.sd <- gather(fcl.sd, "Behavior", "StdDev", 2:10)

# Do same for se's
# *************************************************************

se.func <- function(data){
  
  sd <- apply(data,MARGIN = 2, FUN = sd, na.rm=TRUE)
  se <- sd / sqrt(nrow(data))
  return(se)
}

se1 <- se.func(f.clr1[18:26])
se2 <- se.func(f.clr2[18:26])
se3 <- se.func(f.clr3[18:26])
fcl.se <- rbind(se1,se2,se3)
fcl.se <- as.data.frame(fcl.se, row.names = FALSE)

# Column headers
colnames(fcl.se) <- col.names

# Add in character field and re-order
fcl.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fcl.se <- fcl.se[,c(10,1:9)] 

# Do as before
fcl.se <- gather(fcl.se, "Behavior", "SE", 2:10)

# Now calculate the same on the means and then combine
fclr.mns <- gather(fclr.mns, "Behavior", "Mean", 2:10)
fclr.mns <- cbind(fclr.mns,fcl.sd[,3],fcl.se[,3]) # Combining all summaries together?

# Rename columns
colnames(fclr.mns)[4] <- "StdDev"
colnames(fclr.mns)[5] <- "SE"

# Order the frame
N.order <- order(fclr.mns$Behavior, decreasing=FALSE)
fclr.mns <- fclr.mns[N.order,]
write.csv(fclr.mns, "f_clr_mns.csv")

# Some more subsetting of the output
fcl.pre <- subset(fclr.mns, Phase=="Pre")
fcl.col <- subset(fclr.mns, Phase=="Post (0-3)")
fcl.post <- subset(fclr.mns, Phase=="Post (4+)")

# What are these hard-coded values?
xax <- c(1,1.8,2.6,4,4.8,5.6,7,7.8,8.6,10,10.8,11.6,13,13.8,14.6,16,16.8,17.6,19,19.8,20.6,22,22.8,23.6,25,25.8,26.6)
xpre <- c(1,4,7,10,13,16,19,22,25)
xcol <- c(1.8,4.8,7.8,10.8,13.8,16.8,19.8,22.8,25.8)
xpost <- c(2.6,5.6,8.6,11.6,14.6,17.6,20.6,23.6,26.6)

# ? Averages across each of the different behaviors?
avgpre6 <- fcl.pre$Mean  
sepre6 <- fcl.pre$SE
avgcol6 <- fcl.col$Mean  
secol6 <- fcl.col$SE
avgpost6 <- fcl.post$Mean  
sepost6 <- fcl.post$SE

# ****************************************************************************
# ****************************************************************************
# Now do the same for the combined controls
cmb.ct1 <- subset(bdata, Treatment=="control" & AdjObTime==1)
cmb.ct2 <- subset(bdata, Treatment=="control" & AdjObTime==2)
cmb.ct3 <- subset(bdata, Treatment=="control" & AdjObTime==3)

# Summarize means
# *************************************************************
mn1 <- apply(cmb.ct1[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn2 <- apply(cmb.ct2[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn3 <- apply(cmb.ct3[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
cmb.ct.mns <- rbind(mn1,mn2,mn3)
cmb.ct.mns <- as.data.frame(cmb.ct.mns, row.names = FALSE)
colnames(cmb.ct.mns) <- col.names
cmb.ct.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.ct.mns <- cmb.ct.mns[,c(10,1:9)] 

# Summarize sd
# *************************************************************
sd1 <- apply(cmb.ct1[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd2 <- apply(cmb.ct2[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd3 <- apply(cmb.ct3[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
cmb.ct.sd <- rbind(sd1,sd2,sd3)
cmb.ct.sd <- as.data.frame(cmb.ct.sd, row.names = FALSE)
colnames(cmb.ct.sd) <- col.names
cmb.ct.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.ct.sd <- cmb.ct.sd[,c(10,1:9)] 

# Summarize sd
# *************************************************************
se1 <- se.func(cmb.ct1[18:26])
se2 <- se.func(cmb.ct2[18:26])
se3 <- se.func(cmb.ct3[18:26])
cmb.ct.se <- rbind(se1,se2,se3)
cmb.ct.se <- as.data.frame(cmb.ct.se, row.names = FALSE)
colnames(cmb.ct.se) <- col.names
cmb.ct.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.ct.se <- cmb.ct.se[,c(10,1:9)] 

# Use the gather command
cmb.ct.se <- gather(cmb.ct.se, "Behavior", "SE", 2:10)
cmb.ct.mns <- gather(cmb.ct.mns, "Behavior", "Mean", 2:10)
cmb.ct.sd <- gather(cmb.ct.sd, "Behavior", "StdDev", 2:10)

cmb.ct.mns <- cbind(cmb.ct.mns, cmb.ct.sd[,3], cmb.ct.se[,3])
colnames(cmb.ct.mns)[4] <- "StdDev"
colnames(cmb.ct.mns)[5] <- "SE"

N.order <- order(cmb.ct.mns$Behavior, decreasing=FALSE)
cmb.ct.mns <- cmb.ct.mns[N.order,]

cmb.ct.pre <- subset(cmb.ct.mns, Phase=="Pre")
cmb.ct.col <- subset(cmb.ct.mns, Phase=="Post (0-3)")
cmb.ct.post <- subset(cmb.ct.mns, Phase=="Post (4+)")

avgpre5 <- cmb.ct.pre$Mean  
sepre5 <- cmb.ct.pre$SE
avgcol5 <- cmb.ct.col$Mean  
secol5 <- cmb.ct.col$SE
avgpost5 <- cmb.ct.post$Mean
sepost5 <- cmb.ct.post$SE

# ****************************************************************************
# ****************************************************************************
# And the control Males.....................................................................................................................................................
m.ctl1 <- subset(bdata, Sex=="Male" & Treatment=="control" & AdjObTime==1)
m.ctl2 <- subset(bdata, Sex=="Male" & Treatment=="control" & AdjObTime==2)
m.ctl3 <- subset(bdata, Sex=="Male" & Treatment=="control" & AdjObTime==3)

# Summarize means
# *************************************************************
mn1 <- apply(m.ctl1[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn2 <- apply(m.ctl2[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn3 <- apply(m.ctl3[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
m.ctl.mns <- rbind(mn1,mn2,mn3)
m.ctl.mns <- as.data.frame(m.ctl.mns, row.names = FALSE)
colnames(m.ctl.mns) <- col.names
m.ctl.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.ctl.mns <- m.ctl.mns[,c(10,1:9)] 

# Summarize sd
# *************************************************************
sd1 <- apply(m.ctl1[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd2 <- apply(m.ctl2[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd3 <- apply(m.ctl3[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
m.ctl.sd <- rbind(sd1,sd2,sd3)
m.ctl.sd <- as.data.frame(m.ctl.sd, row.names = FALSE)
colnames(m.ctl.sd) <- col.names
m.ctl.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.ctl.sd <- m.ctl.sd[,c(10,1:9)] 

# Summarize se
# *************************************************************
se1 <- se.func(m.ctl1[18:26])
se2 <- se.func(m.ctl2[18:26])
se3 <- se.func(m.ctl3[18:26])
m.ctl.se <- rbind(se1,se2,se3)
m.ctl.se <- as.data.frame(m.ctl.se, row.names = FALSE)
colnames(m.ctl.se) <- col.names
m.ctl.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.ctl.se <- m.ctl.se[,c(10,1:9)] 

m.ctl.sd <- gather(m.ctl.sd, "Behavior", "StdDev", 2:10)
m.ctl.se <- gather(m.ctl.se, "Behavior", "SE", 2:10)
m.ctl.mns <- gather(m.ctl.mns, "Behavior", "Mean", 2:10)
m.ctl.mns <- cbind(m.ctl.mns, m.ctl.sd[,3], m.ctl.se[,3])
colnames(m.ctl.mns)[4] <- "StdDev"
colnames(m.ctl.mns)[5] <- "SE"

N.order <- order(m.ctl.mns$Behavior, decreasing=FALSE)
m.ctl.mns <- m.ctl.mns[N.order,]

# Why not use aggregate here?  Again, no idea of the purpose here.
m.ctl.pre <- subset(m.ctl.mns, Phase=="Pre")
m.ctl.col <- subset(m.ctl.mns, Phase=="Post (0-3)")
m.ctl.post <- subset(m.ctl.mns, Phase=="Post (4+)")

avgpre4 <- m.ctl.pre$Mean 
sepre4 <- m.ctl.pre$SE
avgcol4 <- m.ctl.col$Mean  
secol4 <- m.ctl.col$SE
avgpost4 <- m.ctl.post$Mean  
sepost4 <- m.ctl.post$SE

# ****************************************************************************
# ****************************************************************************
# Combined Collars......................................................................................................................
cmb.cl1 <- subset(bdata, Treatment !="control" & AdjObTime==1)
cmb.cl2 <- subset(bdata, Treatment !="control" & AdjObTime==2)
cmb.cl3 <- subset(bdata, Treatment !="control" & AdjObTime==3)

# Summarize means
# *************************************************************
mn1 <- apply(cmb.cl1[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn2 <- apply(cmb.cl2[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn3 <- apply(cmb.cl3[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
cmb.cl.mns <- rbind(mn1,mn2,mn3)
cmb.cl.mns <- as.data.frame(cmb.cl.mns, row.names = FALSE)
  colnames(cmb.cl.mns) <- col.names
cmb.cl.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.cl.mns <- cmb.cl.mns[,c(10,1:9)] 

# Summarize sd
# *************************************************************
sd1 <- apply(cmb.cl1[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd2 <- apply(cmb.cl2[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd3 <- apply(cmb.cl3[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
cmb.cl.sd <- rbind(sd1,sd2,sd3)
cmb.cl.sd <- as.data.frame(cmb.cl.sd, row.names = FALSE)
  colnames(cmb.cl.sd) <- col.names
cmb.cl.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.cl.sd <- cmb.cl.sd[,c(10,1:9)] 

# Summarize se
# *************************************************************
se1 <- se.func(cmb.cl1[18:26])
se2 <- se.func(cmb.cl2[18:26])
se3 <- se.func(cmb.cl3[18:26])
cmb.cl.se <- rbind(se1,se2,se3)
cmb.cl.se <- as.data.frame(cmb.cl.se, row.names = FALSE)
  colnames(cmb.cl.se) <- col.names
cmb.cl.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.cl.se <- cmb.cl.se[,c(10,1:9)] 

cmb.cl.sd <- gather(cmb.cl.sd, "Behavior", "StdDev", 2:10)
cmb.cl.se <- gather(cmb.cl.se, "Behavior", "SE", 2:10)
cmb.cl.mns <- gather(cmb.cl.mns, "Behavior", "Mean", 2:10)
cmb.cl.mns <- cbind(cmb.cl.mns, cmb.cl.sd[,3], cmb.cl.se[,3])
  colnames(cmb.cl.mns)[4] <- "StdDev"
  colnames(cmb.cl.mns)[5] <- "SE"

N.order <- order(cmb.cl.mns$Behavior, decreasing=FALSE)
cmb.cl.mns <- cmb.cl.mns[N.order,]

cmb.cl.pre <- subset(cmb.cl.mns, Phase=="Pre")
cmb.cl.col <- subset(cmb.cl.mns, Phase=="Post (0-3)")
cmb.cl.post <- subset(cmb.cl.mns, Phase=="Post (4+)")

avgpre3 <- cmb.cl.pre$Mean  
sepre3 <- cmb.cl.pre$SE
avgcol3 <- cmb.cl.col$Mean  
secol3 <- cmb.cl.col$SE
avgpost3 <- cmb.cl.post$Mean  
sepost3 <- cmb.cl.post$SE

# ****************************************************************************
# ****************************************************************************
# Collared Males......................................................................................................................
m.clr1 <- subset(bdata, Sex=="Male" & Treatment=="Vectronic" & AdjObTime==1)
m.clr2 <- subset(bdata, Sex=="Male" & Treatment=="Vectronic" & AdjObTime==2)
m.clr3 <- subset(bdata, Sex=="Male" & Treatment=="Vectronic" & AdjObTime==3)

# Summarize means
# *************************************************************
mn1 <- apply(m.clr1[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn2 <- apply(m.clr2[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn3 <- apply(m.clr3[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
m.clr.mns <- rbind(mn1,mn2,mn3)
m.clr.mns <- as.data.frame(m.clr.mns, row.names = FALSE)
colnames(m.clr.mns) <- col.names
m.clr.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.clr.mns <- m.clr.mns[,c(10,1:9)] 

# Summarize sd
# *************************************************************
sd1 <- apply(m.clr1[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd2 <- apply(m.clr2[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd3 <- apply(m.clr3[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
m.clr.sd <- rbind(sd1,sd2,sd3)
m.clr.sd <- as.data.frame(m.clr.sd, row.names = FALSE)
colnames(m.clr.sd) <- col.names
m.clr.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.clr.sd <- m.clr.sd[,c(10,1:9)] 

# Summarize se
# *************************************************************
se1 <- se.func(m.clr1[18:26])
se2 <- se.func(m.clr2[18:26])
se3 <- se.func(m.clr3[18:26])
m.clr.se <- rbind(se1,se2,se3)
m.clr.se <- as.data.frame(m.clr.se, row.names = FALSE)
colnames(m.clr.se) <- col.names
m.clr.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.clr.se <- m.clr.se[,c(10,1:9)] 

m.clr.sd <- gather(m.clr.sd, "Behavior", "StdDev", 2:10)
m.clr.se <- gather(m.clr.se, "Behavior", "SE", 2:10)
m.clr.mns <- gather(m.clr.mns, "Behavior", "Mean", 2:10)
m.clr.mns <- cbind(m.clr.mns, m.clr.sd[,3], m.clr.se[,3])
colnames(m.clr.mns)[4] <- "StdDev"
colnames(m.clr.mns)[5] <- "SE"

N.order <- order(m.clr.mns$Behavior, decreasing=FALSE)
m.clr.mns <- m.clr.mns[N.order,]

m.clr.pre <- subset(m.clr.mns, Phase=="Pre")
m.clr.col <- subset(m.clr.mns, Phase=="Post (0-3)")
m.clr.post <- subset(m.clr.mns, Phase=="Post (4+)")

avgpre2 <- m.clr.pre$Mean  
sepre2 <- m.clr.pre$SE
avgcol2 <- m.clr.col$Mean  
secol2 <- m.clr.col$SE
avgpost2 <- m.clr.post$Mean  
sepost2 <- m.clr.post$SE

# ****************************************************************************
# ****************************************************************************
# Control Females.......................................................................................................................
f.ctl1 <- subset(bdata, Sex=="Female" & Treatment=="control" & AdjObTime==1)
f.ctl2 <- subset(bdata, Sex=="Female" & Treatment=="control" & AdjObTime==2)
f.ctl3 <- subset(bdata, Sex=="Female" & Treatment=="control" & AdjObTime==3)

# Summarize means
# *************************************************************
mn1 <- apply(f.ctl1[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn2 <- apply(f.ctl2[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
mn3 <- apply(f.ctl3[18:26],MARGIN = 2, FUN = mean, na.rm=TRUE)
fct.mns <- rbind(mn1,mn2,mn3)
fct.mns <- as.data.frame(fct.mns, row.names = FALSE)
colnames(fct.mns) <- col.names
fct.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fct.mns <- fct.mns[,c(10,1:9)] 

# Summarize sd
# *************************************************************
sd1 <- apply(f.ctl1[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd2 <- apply(f.ctl2[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
sd3 <- apply(f.ctl3[18:26],MARGIN = 2, FUN = sd, na.rm=TRUE)
fct.sd <- rbind(sd1,sd2,sd3)
fct.sd <- as.data.frame(fct.sd, row.names = FALSE)
colnames(fct.sd) <- col.names
fct.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fct.sd <- fct.sd[,c(10,1:9)] 

# Summarize se
# *************************************************************
se1 <- se.func(f.ctl1[18:26])
se2 <- se.func(f.ctl2[18:26])
se3 <- se.func(f.ctl3[18:26])
fct.se <- rbind(se1,se2,se3)
fct.se <- as.data.frame(fct.se, row.names = FALSE)
colnames(fct.se) <- col.names
fct.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fct.se <- fct.se[,c(10,1:9)] 

fct.sd <- gather(fct.sd, "Behavior", "StdDev", 2:10)
fct.se <- gather(fct.se, "Behavior", "SE", 2:10)
fct.mns <- gather(fct.mns, "Behavior", "Mean", 2:10)
fct.mns <- cbind(fct.mns, fct.sd[,3], fct.se[,3])
colnames(fct.mns)[4] <- "StdDev"
colnames(fct.mns)[5] <- "SE"

N.order <- order(fct.mns$Behavior, decreasing=FALSE)
fct.mns <- fct.mns[N.order,]

fct.pre <- subset(fct.mns, Phase=="Pre")
fct.col <- subset(fct.mns, Phase=="Post (0-3)")
fct.post <- subset(fct.mns, Phase=="Post (4+)")

avgpre1 <- fct.pre$Mean  
sepre1 <- fct.pre$SE
avgcol1 <- fct.col$Mean  
secol1 <- fct.col$SE
avgpost1 <- fct.post$Mean  
sepost1 <- fct.post$SE

#################################################################################################################################################################################################
#################################################################################################################################################################################################
# Put the plots together so that they are all in the same window.

par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(4,4.5,2,0))

# Create new column names files...to match how Stephanie plots
col.names1 <- c("FHD","FHU","HS","LAY","SCR","SHD","SHU","SOC","WLK")

# Collared Males
plot(xpre, m.clr.pre$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray60", axes=FALSE, xlab=" ", ylab="Proportion oF Observation Period", main="Collars", cex.main=1.5)
axis(1, at=seq(1.8,25.8,3), labels=col.names1, cex.axis=1.3)
axis(2, at=seq(0,0.75,0.05), cex.axis=1.3)
arrows(xpre, avgpre2-sepre2, xpre, avgpre2+sepre2, length=0.05, angle=90, code=3, lwd=2)

#par(new=TRUE)
#dev.new()
#plot(xcol, m.clr.col$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray80", axes=FALSE, ylab=" ", xlab=" ")
points(xcol,m.clr.col$Mean,pch=22,cex=2.5,col="black",bg="gray80")
arrows(xcol, avgcol2-secol2, xcol, avgcol2+secol2, length=0.05, angle=90, code=3, lwd=2)

par(new=TRUE)
plot(xpost, m.clr.post$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black", cex.lab=1.4, bg="gray100", axes=FALSE, ylab=" ", xlab=" ")
arrows(xpost, avgpost2-sepost2, xpost, avgpost2+sepost2, length=0.05, angle=90, code=3, lwd=2)
abline(v=seq(3.3,24.3,3) , lty=3, lwd=1)
abline(h=0, lty=3, lwd=1)
text(1.8,0.7, "A", cex=1.4)

par(new=FALSE)
# Control Male
plot(xpre, m.ctl.pre$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray60", axes=FALSE, xlab=" ", ylab=" ", main="Controls", cex.main=1.5)
axis(1, at=seq(1.8,25.8,3), labels=c("FHD","FHU","HS","LAY","SCR","SHD","SHU","SOC","WLK"), cex.axis=1.3)
axis(2, at=seq(0,0.75,0.05), cex.axis=1.3)
arrows(xpre, avgpre4-sepre4, xpre, avgpre4+sepre4, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xcol, m.ctl.col$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray80", axes=FALSE, ylab=" ", xlab=" ")
arrows(xcol, avgcol4-secol4, xcol, avgcol4+secol4, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xpost, m.ctl.post$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray100", axes=FALSE, ylab=" ", xlab=" ")
arrows(xpost, avgpost4-sepost4, xpost, avgpost4+sepost4, length=0.05, angle=90, code=3, lwd=2)
legend("topright", c("Pre","Collar","Post"), fill=c("gray60","gray80","gray100"), bty="n", cex=1.4, bg="white")
abline(v=seq(3.3,24.3,3) , lty=3, lwd=1)
abline(h=0, lty=3, lwd=1)
text(1.8,0.7, "B", cex=1.4)
par(new=TRUE)
plot(13.8, 0.75, type="p", pch=8, cex=1.2, axes=FALSE, ylab=" ", xlab=" ", ylim=c(0,0.75), xlim=c(1,27))
par(new=FALSE)

# Collared Females
plot(xpre, fcl.pre$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black", cex.lab=1.4, bg="gray60", axes=FALSE, xlab=" ", ylab="Proportion oF Observation Period")
axis(1, at=seq(1.8,25.8,3), labels=c("FHD","FHU","HS","LAY","SCR","SHD","SHU","SOC","WLK"), cex.axis=1.3)
axis(2, at=seq(0,0.75,0.05), cex.axis=1.3)
arrows(xpre, avgpre6-sepre6, xpre, avgpre6+sepre6, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xcol, fcl.col$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray80", axes=FALSE, ylab=" ", xlab=" ")
arrows(xcol, avgcol6-secol6, xcol, avgcol6+secol6, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xpost, fcl.post$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray100", axes=FALSE, ylab=" ", xlab=" ")
arrows(xpost, avgpost6-sepost6, xpost, avgpost6+sepost6, length=0.05, angle=90, code=3, lwd=2)
abline(v=seq(3.3,24.3,3) , lty=3, lwd=1)
abline(h=0, lty=3, lwd=1)
text(1.8,0.7, "C", cex=1.4)
par(new=TRUE)
plot(c(1.8, 7.8, 13.8, 16.8, 22.8, 25.8), rep(0.75, 6), type="p", pch=8, cex=1.2, axes=FALSE, ylab=" ", xlab=" ", ylim=c(0,0.75), xlim=c(1,27))
par(new=FALSE)

# Control Females
plot(xpre, fct.pre$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black", cex.lab=1.4, bg="gray60", axes=FALSE, xlab=" ", ylab=" ")
axis(1, at=seq(1.8,25.8,3), labels=c("FHD","FHU","HS","LAY","SCR","SHD","SHU","SOC","WLK"), cex.axis=1.3)
axis(2, at=seq(0,0.75,0.05), cex.axis=1.3)
arrows(xpre, avgpre1-sepre1, xpre, avgpre1+sepre1, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xcol, fct.col$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray80", axes=FALSE, ylab=" ", xlab=" ")
arrows(xcol, avgcol1-secol1, xcol, avgcol1+secol1, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xpost, fct.post$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray100", axes=FALSE, ylab=" ", xlab=" ")
arrows(xpost, avgpost1-sepost1, xpost, avgpost1+sepost1, length=0.05, angle=90, code=3, lwd=2)
abline(v=seq(3.3,24.3,3) , lty=3, lwd=1)
abline(h=0, lty=3, lwd=1)
text(1.8,0.7, "D", cex=1.4)
par(new=FALSE)

# Combined Collars
plot(xpre, cmb.cl.pre$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, cex.lab=1.4, col="black", bg="gray60", axes=FALSE, xlab="Behavior", ylab="Proportion of Observation Period")
axis(1, at=seq(1.8,25.8,3), labels=c("FHD","FHU","HS","LAY","SCR","SHD","SHU","SOC","WLK"), cex.axis=1.3)
axis(2, at=seq(0,0.75,0.05), cex.axis=1.3)
arrows(xpre, avgpre3-sepre3, xpre, avgpre3+sepre3, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xcol, cmb.cl.col$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, cex.lab=1.4, col="black", bg="gray80", axes=FALSE, ylab=" ", xlab=" ")
arrows(xcol, avgcol3-secol3, xcol, avgcol3+secol3, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xpost, cmb.cl.post$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5,cex.lab=1.4,  col="black", bg="gray100", axes=FALSE, ylab=" ", xlab=" ")
arrows(xpost, avgpost3-sepost3, xpost, avgpost3+sepost3, length=0.05, angle=90, code=3, lwd=2)
abline(v=seq(3.3,24.3,3) , lty=3, lwd=1)
abline(h=0, lty=3, lwd=1)
text(1.8,0.7, "E", cex=1.4)
par(new=TRUE)
plot(c(4.8, 7.8, 10.8, 13.8, 16.8, 22.8, 25.8), rep(0.75, 7), type="p", pch=8, cex=1.2, axes=FALSE, ylab=" ", xlab=" ", ylim=c(0,0.75), xlim=c(1,27))
par(new=FALSE)

# Combined Controls
plot(xpre, cmb.ct.pre$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22,cex.lab=1.4,  cex=2.5, col="black", bg="gray60", axes=FALSE, xlab="Behavior", ylab=" ")
axis(1, at=seq(1.8,25.8,3), labels=c("FHD","FHU","HS","LAY","SCR","SHD","SHU","SOC","WLK"), cex.axis=1.3)
axis(2, at=seq(0,0.75,0.05), cex.axis=1.3)
arrows(xpre, avgpre5-sepre5, xpre, avgpre5+sepre5, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xcol, cmb.ct.col$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5,cex.lab=1.4,  col="black", bg="gray80", axes=FALSE, ylab=" ", xlab=" ")
arrows(xcol, avgcol5-secol5, xcol, avgcol5+secol5, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xpost, cmb.ct.post$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5,cex.lab=1.4,  col="black", bg="gray100", axes=FALSE, ylab=" ", xlab=" ")
arrows(xpost, avgpost5-sepost5, xpost, avgpost5+sepost5, length=0.05, angle=90, code=3, lwd=2)
abline(v=seq(3.3,24.3,3) , lty=3, lwd=1)
abline(h=0, lty=3, lwd=1)
text(1.8,0.7, "F", cex=1.4)
par(new=TRUE)
plot(10.8, 0.75, type="p", pch=8, cex=1.2, axes=FALSE, ylab=" ", xlab=" ", ylim=c(0,0.75), xlim=c(1,27))
par(new=FALSE)
