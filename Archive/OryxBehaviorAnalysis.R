#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Project: SHO Stress and Behavior Analysis
# Date: 10 November 2016
# Author: Stephanie Cunningham
# Description: Plotting behavior means

#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Clear the cache
rm(list=ls())

# Set working directory
# setwd("Z:/Interns/StephanieCunningham/SHO Behavior/SCBI/BehaviorData")
# setwd("/Volumes/STEPH 1TB/SCBI Work/SHO Behavior and Stress/SCBI/BehaviorData")
setwd("C:/Users/Jared/Dropbox (Smithsonian)/Projects/Oryx/StressAnalysis/Behavior")

# Read in file that was fixed in Excel
bdata <- read.csv("Behavior.Nov3.csv")

# Fix the data/time fields
bdata$TimeStart <- as.POSIXct(strptime(paste0(bdata$Date," ",bdata$TimeStart), format="%m/%d/%Y %H:%M"))
bdata$TimeEnd <- as.POSIXct(strptime(paste0(bdata$Date," ",bdata$TimeEnd), format="%m/%d/%Y %H:%M"))

# Re-order dataframe so all the behavior data is at the end
bdata <- bdata[,c(1:15,26:27,16:25)]

# Collared Females
# I'm not sure what AdjObTime 1, 2, and 3 are.....................................................................
f.clr1 <- subset(bdata, Sex=="Female" & (Treatment=="Vectronic" | Treatment=="ATS") & AdjObTime==1)
f.clr2 <- subset(bdata, Sex=="Female" & (Treatment=="Vectronic" | Treatment=="ATS") & AdjObTime==2)
f.clr3 <- subset(bdata, Sex=="Female" & (Treatment=="Vectronic" | Treatment=="ATS") & AdjObTime==3)

# Creating a blank dataframe 
fclr.mns <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# How to add values to data frame
# Could use apply to make this a lot easier
pre <- f.clr1[18:27]
post2 <- f.clr2[18:27]
post3 <- f.clr3[18:27]
for (i in 1:length(fclr.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  fclr.mns[1,i] <- mean(behaviors.pre, na.rm=TRUE)
  fclr.mns[2,i] <- mean(behaviors.post2, na.rm=TRUE)
  fclr.mns[3,i] <- mean(behaviors.post3, na.rm=TRUE)
  #fclr.mns[1,i] <- mn1
  #fclr.mns[2,i] <- mn2
  #fclr.mns[3,i] <- mn3
}

fclr.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fclr.mns <- fclr.mns[,c(10,1:9)] 

fcl.sd <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- f.clr1[18:27]
post2 <- f.clr2[18:27]
post3 <- f.clr3[18:27]

for (i in 1:length(fclr.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  fcl.sd[1,i] <- sd(behaviors.pre, na.rm=TRUE)
  fcl.sd[2,i] <- sd(behaviors.post2, na.rm=TRUE)
  fcl.sd[3,i] <- sd(behaviors.post3, na.rm=TRUE)
  #fcl.sd[1,i] <- sd1
  #fcl.sd[2,i] <- sd2
  #fcl.sd[3,i] <- sd3
}

fcl.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fcl.sd <- fcl.sd[,c(10,1:9)] 
fcl.sd <- gather(fcl.sd, "Behavior", "StdDev", 2:10)
StdDev <- c(fcl.sd$StdDev)

fcl.se <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- f.clr1[18:27]
post2 <-f.clr2[18:27]
post3 <-f.clr3[18:27]
for (i in 1:length(fclr.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  l1 <- nrow(pre)
  l2 <- nrow(post2)
  l3 <- nrow(post3)
  se1 <- sd1/sqrt(l1)
  se2 <- sd2/sqrt(l2)
  se3 <- sd3/sqrt(l3)
  fcl.se[1,i] <- se1
  fcl.se[2,i] <- se2
  fcl.se[3,i] <- se3
}

fcl.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fcl.se <- fcl.se[,c(10,1:9)] 
fcl.se <- gather(fcl.se, "Behavior", "SE", 2:10)
SE <- c(fcl.se$SE)
fclr.mns <- gather(fclr.mns, "Behavior", "Mean", 2:10)
fclr.mns <- cbind(fclr.mns, StdDev, SE)

N.order <- order(fclr.mns$Behavior, decreasing=FALSE)
fclr.mns <- fclr.mns[N.order,]
write.csv(fclr.mns, "f_clr_mns.csv")


fcl.pre <- subset(fclr.mns, Phase=="Pre")
fcl.col <- subset(fclr.mns, Phase=="Post (0-3)")
fcl.post <- subset(fclr.mns, Phase=="Post (4+)")

xax <- c(1,1.8,2.6,4,4.8,5.6,7,7.8,8.6,10,10.8,11.6,13,13.8,14.6,16,16.8,17.6,19,19.8,20.6,22,22.8,23.6,25,25.8,26.6)
xpre <- c(1,4,7,10,13,16,19,22,25)
xcol <- c(1.8,4.8,7.8,10.8,13.8,16.8,19.8,22.8,25.8)
xpost <- c(2.6,5.6,8.6,11.6,14.6,17.6,20.6,23.6,26.6)

avgpre6 <- c(fcl.pre$Mean)  
sepre6 <- c(fcl.pre$SE)
avgcol6 <- c(fcl.col$Mean)  
secol6 <- c(fcl.col$SE)
avgpost6 <- c(fcl.post$Mean)  
sepost6 <- c(fcl.post$SE)


# Combined Controls...............................................................................................................................
cmb.ct1 <- subset(bdata, Treatment=="control" & AdjObTime==1)
cmb.ct2 <- subset(bdata, Treatment=="control" & AdjObTime==2)
cmb.ct3 <- subset(bdata, Treatment=="control" & AdjObTime==3)

cmb.ct.mns <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- cmb.ct1[18:27]
post2 <-cmb.ct2[18:27]
post3 <-cmb.ct3[18:27]
for (i in 1:length(cmb.ct.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  mn1 <- mean(behaviors.pre, na.rm=TRUE)
  mn2 <- mean(behaviors.post2, na.rm=TRUE)
  mn3 <- mean(behaviors.post3, na.rm=TRUE)
  cmb.ct.mns[1,i] <- mn1
  cmb.ct.mns[2,i] <- mn2
  cmb.ct.mns[3,i] <- mn3
}

cmb.ct.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.ct.mns <- cmb.ct.mns[,c(10,1:9)] 


cmb.ct.sd <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- cmb.ct1[18:27]
post2 <-cmb.ct2[18:27]
post3 <-cmb.ct3[18:27]
for (i in 1:length(cmb.ct.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  cmb.ct.sd[1,i] <- sd1
  cmb.ct.sd[2,i] <- sd2
  cmb.ct.sd[3,i] <- sd3
}

cmb.ct.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.ct.sd <- cmb.ct.sd[,c(10,1:9)] 

cmb.ct.se <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- cmb.ct1[18:27]
post2 <- cmb.ct2[18:27]
post3 <- cmb.ct3[18:27]
for (i in 1:length(cmb.ct.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  l1 <- nrow(pre)
  l2 <- nrow(post2)
  l3 <- nrow(post3)
  se1 <- sd1/sqrt(l1)
  se2 <- sd2/sqrt(l2)
  se3 <- sd3/sqrt(l3)
  cmb.ct.se[1,i] <- se1
  cmb.ct.se[2,i] <- se2
  cmb.ct.se[3,i] <- se3
}

cmb.ct.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.ct.se <- cmb.ct.se[,c(10,1:9)] 
cmb.ct.se <- gather(cmb.ct.se, "Behavior", "SE", 2:10)
SE <- c(cmb.ct.se$SE)

cmb.ct.mns <- gather(cmb.ct.mns, "Behavior", "Mean", 2:10)
cmb.ct.sd <- gather(cmb.ct.sd, "Behavior", "StdDev", 2:10)
StdDev <- c(cmb.ct.sd$StdDev)

cmb.ct.mns <- cbind(cmb.ct.mns, StdDev, SE)
N.order <- order(cmb.ct.mns$Behavior, decreasing=FALSE)
cmb.ct.mns <- cmb.ct.mns[N.order,]


cmb.ct.pre <- subset(cmb.ct.mns, Phase=="Pre")
cmb.ct.col <- subset(cmb.ct.mns, Phase=="Post (0-3)")
cmb.ct.post <- subset(cmb.ct.mns, Phase=="Post (4+)")

xax <- c(1,1.8,2.6,4,4.8,5.6,7,7.8,8.6,10,10.8,11.6,13,13.8,14.6,16,16.8,17.6,19,19.8,20.6,22,22.8,23.6,25,25.8,26.6)
xpre <- c(1,4,7,10,13,16,19,22,25)
xcol <- c(1.8,4.8,7.8,10.8,13.8,16.8,19.8,22.8,25.8)
xpost <- c(2.6,5.6,8.6,11.6,14.6,17.6,20.6,23.6,26.6)

avgpre5 <- c(cmb.ct.pre$Mean)  
sepre5 <- c(cmb.ct.pre$SE)
avgcol5 <- c(cmb.ct.col$Mean)  
secol5 <- c(cmb.ct.col$SE)
avgpost5 <- c(cmb.ct.post$Mean)  
sepost5 <- c(cmb.ct.post$SE)



# Control Male.....................................................................................................................................................
m.ctl1 <- subset(bdata, Sex=="Male" & Treatment=="control" & AdjObTime==1)
m.ctl2 <- subset(bdata, Sex=="Male" & Treatment=="control" & AdjObTime==2)
m.ctl3 <- subset(bdata, Sex=="Male" & Treatment=="control" & AdjObTime==3)

m.ctl.mns <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- m.ctl1[18:27]
post2 <-m.ctl2[18:27]
post3 <-m.ctl3[18:27]
for (i in 1:length(m.ctl.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  mn1 <- mean(behaviors.pre, na.rm=TRUE)
  mn2 <- mean(behaviors.post2, na.rm=TRUE)
  mn3 <- mean(behaviors.post3, na.rm=TRUE)
  m.ctl.mns[1,i] <- mn1
  m.ctl.mns[2,i] <- mn2
  m.ctl.mns[3,i] <- mn3
}

m.ctl.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.ctl.mns <- m.ctl.mns[,c(10,1:9)] 


m.ctl.sd <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- m.ctl1[18:27]
post2 <- m.ctl2[18:27]
post3 <- m.ctl3[18:27]
for (i in 1:length(m.ctl.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  m.ctl.sd[1,i] <- sd(behaviors.pre, na.rm=TRUE)
  m.ctl.sd[2,i] <- sd(behaviors.post2, na.rm=TRUE)
  m.ctl.sd[3,i] <- sd(behaviors.post3, na.rm=TRUE)
  #m.ctl.sd[1,i] <- sd1
  #m.ctl.sd[2,i] <- sd2
  #m.ctl.sd[3,i] <- sd3
}

m.ctl.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.ctl.sd <- m.ctl.sd[,c(10,1:9)] 
m.ctl.sd <- gather(m.ctl.sd, "Behavior", "StdDev", 2:10)
StdDev <- c(m.ctl.sd$StdDev)

m.ctl.se <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- m.ctl1[18:27]
post2 <- m.ctl2[18:27]
post3 <- m.ctl3[18:27]
for (i in 1:length(m.ctl.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  l1 <- nrow(pre)
  l2 <- nrow(post2)
  l3 <- nrow(post3)
  se1 <- sd1/sqrt(l1)
  se2 <- sd2/sqrt(l2)
  se3 <- sd3/sqrt(l3)
  m.ctl.se[1,i] <- se1
  m.ctl.se[2,i] <- se2
  m.ctl.se[3,i] <- se3
}

m.ctl.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.ctl.se <- m.ctl.se[,c(10,1:9)] 
m.ctl.se <- gather(m.ctl.se, "Behavior", "SE", 2:10)
SE <- c(m.ctl.se$SE)
m.ctl.mns <- gather(m.ctl.mns, "Behavior", "Mean", 2:10)
m.ctl.mns <- cbind(m.ctl.mns, StdDev, SE)
N.order <- order(m.ctl.mns$Behavior, decreasing=FALSE)
m.ctl.mns <- m.ctl.mns[N.order,]


m.ctl.pre <- subset(m.ctl.mns, Phase=="Pre")
m.ctl.col <- subset(m.ctl.mns, Phase=="Post (0-3)")
m.ctl.post <- subset(m.ctl.mns, Phase=="Post (4+)")

xax <- c(1,1.8,2.6,4,4.8,5.6,7,7.8,8.6,10,10.8,11.6,13,13.8,14.6,16,16.8,17.6,19,19.8,20.6,22,22.8,23.6,25,25.8,26.6)
xpre <- c(1,4,7,10,13,16,19,22,25)
xcol <- c(1.8,4.8,7.8,10.8,13.8,16.8,19.8,22.8,25.8)
xpost <- c(2.6,5.6,8.6,11.6,14.6,17.6,20.6,23.6,26.6)

avgpre4 <- c(m.ctl.pre$Mean)  
sepre4 <- c(m.ctl.pre$SE)
avgcol4 <- c(m.ctl.col$Mean)  
secol4 <- c(m.ctl.col$SE)
avgpost4 <- c(m.ctl.post$Mean)  
sepost4 <- c(m.ctl.post$SE)


# Combined Collars......................................................................................................................
cmb.cl1 <- subset(bdata, (Treatment=="Vectronic" | Treatment=="ATS") & AdjObTime==1)
cmb.cl2 <- subset(bdata, (Treatment=="Vectronic" | Treatment=="ATS") & AdjObTime==2)
cmb.cl3 <- subset(bdata, (Treatment=="Vectronic" | Treatment=="ATS") & AdjObTime==3)

cmb.cl.mns <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- cmb.cl1[18:27]
post2 <-cmb.cl2[18:27]
post3 <-cmb.cl3[18:27]
for (i in 1:length(cmb.cl.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  mn1 <- mean(behaviors.pre, na.rm=TRUE)
  mn2 <- mean(behaviors.post2, na.rm=TRUE)
  mn3 <- mean(behaviors.post3, na.rm=TRUE)
  cmb.cl.mns[1,i] <- mn1
  cmb.cl.mns[2,i] <- mn2
  cmb.cl.mns[3,i] <- mn3
}

cmb.cl.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.cl.mns <- cmb.cl.mns[,c(10,1:9)] 

cmb.cl.sd <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- cmb.cl1[18:27]
post2 <- cmb.cl2[18:27]
post3 <- cmb.cl3[18:27]
for (i in 1:length(cmb.cl.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  cmb.cl.sd[1,i] <- sd1
  cmb.cl.sd[2,i] <- sd2
  cmb.cl.sd[3,i] <- sd3
}

cmb.cl.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.cl.sd <- cmb.cl.sd[,c(10,1:9)] 
cmb.cl.sd <- gather(cmb.cl.sd, "Behavior", "StdDev", 2:10)
StdDev <- c(cmb.cl.sd$StdDev)

cmb.cl.se <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- cmb.cl1[18:27]
post2 <- cmb.cl2[18:27]
post3 <- cmb.cl3[18:27]
for (i in 1:length(cmb.cl.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  l1 <- nrow(pre)
  l2 <- nrow(post2)
  l3 <- nrow(post3)
  se1 <- sd1/sqrt(l1)
  se2 <- sd2/sqrt(l2)
  se3 <- sd3/sqrt(l3)
  cmb.cl.se[1,i] <- se1
  cmb.cl.se[2,i] <- se2
  cmb.cl.se[3,i] <- se3
}

cmb.cl.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
cmb.cl.se <- cmb.cl.se[,c(10,1:9)] 
cmb.cl.se <- gather(cmb.cl.se, "Behavior", "SE", 2:10)
SE <- c(cmb.cl.se$SE)
cmb.cl.mns <- gather(cmb.cl.mns, "Behavior", "Mean", 2:10)

cmb.cl.mns <- cbind(cmb.cl.mns, StdDev, SE)
N.order <-order(cmb.cl.mns$Behavior,decreasing=FALSE)
cmb.cl.mns <- cmb.cl.mns[N.order,]


cmb.cl.pre <- subset(cmb.cl.mns, Phase=="Pre")
cmb.cl.col <- subset(cmb.cl.mns, Phase=="Post (0-3)")
cmb.cl.post <- subset(cmb.cl.mns, Phase=="Post (4+)")

xax <- c(1,1.8,2.6,4,4.8,5.6,7,7.8,8.6,10,10.8,11.6,13,13.8,14.6,16,16.8,17.6,19,19.8,20.6,22,22.8,23.6,25,25.8,26.6)
xpre <- c(1,4,7,10,13,16,19,22,25)
xcol <- c(1.8,4.8,7.8,10.8,13.8,16.8,19.8,22.8,25.8)
xpost <- c(2.6,5.6,8.6,11.6,14.6,17.6,20.6,23.6,26.6)

avgpre3 <- c(cmb.cl.pre$Mean)  
sepre3 <- c(cmb.cl.pre$SE)
avgcol3 <- c(cmb.cl.col$Mean)  
secol3 <- c(cmb.cl.col$SE)
avgpost3 <- c(cmb.cl.post$Mean)  
sepost3 <- c(cmb.cl.post$SE)



# Collared Males......................................................................................................................
m.clr1 <- subset(bdata, Sex=="Male" & Treatment=="Vectronic" & AdjObTime==1)
m.clr2 <- subset(bdata, Sex=="Male" & Treatment=="Vectronic" & AdjObTime==2)
m.clr3 <- subset(bdata, Sex=="Male" & Treatment=="Vectronic" & AdjObTime==3)

m.clr.mns <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- m.clr1[18:27]
post2 <- m.clr2[18:27]
post3 <- m.clr3[18:27]
for (i in 1:length(m.clr.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  mn1 <- mean(behaviors.pre, na.rm=TRUE)
  mn2 <- mean(behaviors.post2, na.rm=TRUE)
  mn3 <- mean(behaviors.post3, na.rm=TRUE)
  m.clr.mns[1,i] <- mn1
  m.clr.mns[2,i] <- mn2
  m.clr.mns[3,i] <- mn3
}

m.clr.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.clr.mns <- m.clr.mns[,c(10,1:9)] 

m.clr.sd <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- m.clr1[18:27]
post2 <- m.clr2[18:27]
post3 <- m.clr3[18:27]
for (i in 1:length(m.clr.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  m.clr.sd[1,i] <- sd1
  m.clr.sd[2,i] <- sd2
  m.clr.sd[3,i] <- sd3
}

m.clr.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.clr.sd <- m.clr.sd[,c(10,1:9)] 
m.clr.sd <- gather(m.clr.sd, "Behavior", "StdDev", 2:10)
StdDev <- c(m.clr.sd$StdDev)

m.clr.se <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- m.clr1[18:27]
post2 <- m.clr2[18:27]
post3 <- m.clr3[18:27]
for (i in 1:length(m.clr.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  l1 <- nrow(pre)
  l2 <- nrow(post2)
  l3 <- nrow(post3)
  se1 <- sd1/sqrt(l1)
  se2 <- sd2/sqrt(l2)
  se3 <- sd3/sqrt(l3)
  m.clr.se[1,i] <- se1
  m.clr.se[2,i] <- se2
  m.clr.se[3,i] <- se3
}

m.clr.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
m.clr.se <- m.clr.se[,c(10,1:9)] 
m.clr.se <- gather(m.clr.se, "Behavior", "SE", 2:10)
SE <- c(m.clr.se$SE)

m.clr.mns <- gather(m.clr.mns, "Behavior", "Mean", 2:10)

m.clr.mns <- cbind(m.clr.mns, StdDev, SE)
N.order <-order(m.clr.mns$Behavior,decreasing=FALSE)
m.clr.mns <- m.clr.mns[N.order,]


m.clr.pre <- subset(m.clr.mns, Phase=="Pre")
m.clr.col <- subset(m.clr.mns, Phase=="Post (0-3)")
m.clr.post <- subset(m.clr.mns, Phase=="Post (4+)")

xax <- c(1,1.8,2.6,4,4.8,5.6,7,7.8,8.6,10,10.8,11.6,13,13.8,14.6,16,16.8,17.6,19,19.8,20.6,22,22.8,23.6,25,25.8,26.6)
xpre <- c(1,4,7,10,13,16,19,22,25)
xcol <- c(1.8,4.8,7.8,10.8,13.8,16.8,19.8,22.8,25.8)
xpost <- c(2.6,5.6,8.6,11.6,14.6,17.6,20.6,23.6,26.6)

avgpre2 <- c(m.clr.pre$Mean)  
sepre2 <- c(m.clr.pre$SE)
avgcol2 <- c(m.clr.col$Mean)  
secol2 <- c(m.clr.col$SE)
avgpost2 <- c(m.clr.post$Mean)  
sepost2 <- c(m.clr.post$SE)



# Control Females.......................................................................................................................
f.ctl1 <- subset(bdata, Sex=="Female" & Treatment=="control" & AdjObTime==1)
f.ctl2 <- subset(bdata, Sex=="Female" & Treatment=="control" & AdjObTime==2)
f.ctl3 <- subset(bdata, Sex=="Female" & Treatment=="control" & AdjObTime==3)

fct.mns <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- f.ctl1[18:27]
post2 <-f.ctl2[18:27]
post3 <-f.ctl3[18:27]
for (i in 1:length(fct.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  mn1 <- mean(behaviors.pre, na.rm=TRUE)
  mn2 <- mean(behaviors.post2, na.rm=TRUE)
  mn3 <- mean(behaviors.post3, na.rm=TRUE)
  fct.mns[1,i] <- mn1
  fct.mns[2,i] <- mn2
  fct.mns[3,i] <- mn3
}

fct.mns$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fct.mns <- fct.mns[,c(10,1:9)] 

fct.sd <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- f.ctl1[18:27]
post2 <-f.ctl2[18:27]
post3 <-f.ctl3[18:27]
for (i in 1:length(fct.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  fct.sd[1,i] <- sd1
  fct.sd[2,i] <- sd2
  fct.sd[3,i] <- sd3
}

fct.sd$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fct.sd <- fct.sd[,c(10,1:9)] 
fct.sd <- gather(fct.sd, "Behavior", "StdDev", 2:10)
StdDev <- c(fct.sd$StdDev)

fct.se <- data.frame(
  SHU=numeric(),
  SHD=numeric(),
  Laying=numeric(),
  Headshaking=numeric(),
  Walking=numeric(),
  FHU=numeric(),
  FHD=numeric(),
  Scratch=numeric(),
  Social=numeric()
)

# Add values to data frame
pre <- f.ctl1[18:27]
post2 <-f.ctl2[18:27]
post3 <-f.ctl3[18:27]
for (i in 1:length(fct.mns[,1:9])){
  behaviors.pre <- pre[,i]
  behaviors.post2 <- post2[,i]
  behaviors.post3 <- post3[,i]
  sd1 <- sd(behaviors.pre, na.rm=TRUE)
  sd2 <- sd(behaviors.post2, na.rm=TRUE)
  sd3 <- sd(behaviors.post3, na.rm=TRUE)
  l1 <- nrow(pre)
  l2 <- nrow(post2)
  l3 <- nrow(post3)
  se1 <- sd1/sqrt(l1)
  se2 <- sd2/sqrt(l2)
  se3 <- sd3/sqrt(l3)
  fct.se[1,i] <- se1
  fct.se[2,i] <- se2
  fct.se[3,i] <- se3
}

fct.se$Phase <- c("Pre", "Post (0-3)", "Post (4+)")
fct.se <- fct.se[,c(10,1:9)] 
fct.se <- gather(fct.se, "Behavior", "SE", 2:10)
SE <- c(fct.se$SE)
fct.mns <- gather(fct.mns, "Behavior", "Mean", 2:10)

fct.mns <- cbind(fct.mns, StdDev, SE)
N.order <-order(fct.mns$Behavior,decreasing=FALSE)
fct.mns <- fct.mns[N.order,]

fct.pre <- subset(fct.mns, Phase=="Pre")
fct.col <- subset(fct.mns, Phase=="Post (0-3)")
fct.post <- subset(fct.mns, Phase=="Post (4+)")

xax <- c(1,1.8,2.6,4,4.8,5.6,7,7.8,8.6,10,10.8,11.6,13,13.8,14.6,16,16.8,17.6,19,19.8,20.6,22,22.8,23.6,25,25.8,26.6)
xpre <- c(1,4,7,10,13,16,19,22,25)
xcol <- c(1.8,4.8,7.8,10.8,13.8,16.8,19.8,22.8,25.8)
xpost <- c(2.6,5.6,8.6,11.6,14.6,17.6,20.6,23.6,26.6)

avgpre1 <- c(fct.pre$Mean)  
sepre1 <- c(fct.pre$SE)
avgcol1 <- c(fct.col$Mean)  
secol1 <- c(fct.col$SE)
avgpost1 <- c(fct.post$Mean)  
sepost1 <- c(fct.post$SE)

#################################################################################################################################################################################################
#################################################################################################################################################################################################
# Put the plots together so that they are all in the same window.

par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(4,4.5,2,0))
# Collared Males
plot(xpre, m.clr.pre$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray60", axes=FALSE, xlab=" ", ylab="Proportion oF Observation Period", main="Collars", cex.main=1.5)
axis(1, at=seq(1.8,25.8,3), labels=c("FHD","FHU","HS","LAY","SCR","SHD","SHU","SOC","WLK"), cex.axis=1.3)
axis(2, at=seq(0,0.75,0.05), cex.axis=1.3)
arrows(xpre, avgpre2-sepre2, xpre, avgpre2+sepre2, length=0.05, angle=90, code=3, lwd=2)
par(new=TRUE)
plot(xcol, m.clr.col$Mean, type="p", ylim=c(0,0.75), xlim=c(1,27), pch=22, cex=2.5, col="black",cex.lab=1.4,  bg="gray80", axes=FALSE, ylab=" ", xlab=" ")
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












