### LOAD IMPORTANT LIBRARIES
library(matrixStats)

setwd("E:/Dropbox/NMDA/")
setwd("C:/Users/Kevin Caref/Dropbox/NMDA/")

### DEFINE ALL THE IMPORTANT FOLDERS

funcdirect <- paste(getwd(), "/R functions/", sep="")
datafolder <- paste(getwd(), "/EXP3_NAc FR acquisition/MedPC files/", sep="")
dataForRdir <- paste(getwd(), "/EXP3_NAc FR acquisition/Data for R/", sep="")
dataForRCumulative <- paste(getwd(), "/EXP3_NAc FR acquisition/Data for R cumulative/", sep="")
behGraphFolder <- paste(getwd(), "/EXP3_NAc FR acquisition/Graphs/Behavior/", sep="")
CPfuncFolder <- paste(funcdirect, 'Change_Point-master/', sep="")
CPGraphFolder <- paste(getwd(), "/EXP3_NAc FR acquisition/Graphs/Behavior/Change point/", sep="")
MixedGraphFolder <- paste(getwd(), "/EXP3_NAc FR acquisition/Graphs/Mixed/", sep="")
NEXfiles <- paste(getwd(), "/EXP3_NAc FR acquisition/NEX files/", sep="")

### Load necessary functions
load(file=paste(funcdirect, "MedPCextract.r", sep=""))
load(file=paste(funcdirect, "mpcextract_blockSingle.Rfunc", sep=""))
load(file=paste(funcdirect, "CPextract.R", sep=""))
load(file=paste(funcdirect, "cumulativeIndGraphs.R", sep=""))
load(file=paste(funcdirect, "PerformanceFromCP.R", sep=""))
load(file=paste(funcdirect, "PrePostCP_Perf.R", sep=""))
load(file=paste(funcdirect, "avgPerfByBin.R", sep=""))
load(file=paste(funcdirect, "CPextractMultipleCrit.R", sep=""))
load(file=paste(funcdirect, "neuralhist.r", sep=""))
load(file=paste(funcdirect, "FRbyNEURONbyBINcue.r", sep=""))
load(file=paste(funcdirect, "plotFRandCP.r", sep=""))
load(file=paste(funcdirect, "errBars.r", sep=""))
load(file=paste(funcdirect, "errCloud.r", sep=""))
load(file=paste(funcdirect, "masterDFsummary.r", sep=""))
#load(file=paste(funcdirect, "plotFRandCPhistogram.r", sep=""))


#### GENERATE IMPORTANT DATA OBJECTS
MedPCextract(oneITIwdw=T, funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative, cuelength=10, consumeRewWdw = 4)

# Load important behavior-related objects
files <- paste(dataForRdir, list.files(dataForRdir), sep="")
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep="")
for(i in 1:length(files)){load(files[[i]])}
for(i in 1:length(filesCum)){load(filesCum[[i]])}


#### BEHAVIORAL GRAPHS

# RASTERS
load(file=paste(funcdirect, "MAKERASTER.r", sep=""))
MAKERASTER(i=x, data=alldata, idxdata=csacqidx)
MAKERASTER(i=4, data=alldata, idxdata=csacqidx)


#### ALL OF THIS CODE SERVES THE PURPOSE OF FIGURING OUT WHAT SEGMENT OF THE ITI TO USE AS BASELINE BEHAVIOR FOR THE PERFORMANCE INDEX CALCULATION

# POST CUE-END ENTRIES (this is useful to decide how much of the ITI NOT to sample when calculating ITI 'latency')
PostDSentries <- lapply(seq(1, length(alldata)), function(x){
  targetEntries <- alldata[[x]]$PostCSplusEndEntries
  return(targetEntries)
})

PostDSexits <- lapply(seq(1, length(alldata)), function(x){
  targetExits <- alldata[[x]]$PostCSplusEndExits
  return(targetExits)
})

# Figure out in how many trials, the animal is inside the compartment when the 10s precue window starts
Impinge <- sapply(seq(1, length(alldata)), function(x){
  entries <- alldata[[x]]$receptacleentries
  exits <- alldata[[x]]$receptacleexits
  pairs <- lapply(seq(1, length(entries)), function(y){
    a <- c(entries[y], exits[y])
  })
  trialStarts <- alldata[[x]]$trialStart
  ITIprecue <- alldata[[x]]$allCues-10
  
  ITItoEntryIdx <- findInterval(ITIprecue, entries)
  toTest <- pairs[ITItoEntryIdx]
  
  #Select trials based on what cue came on before (modify accordingly). 1=DS, 2=NS.
  preTrialCue <- c(0, alldata[[x]]$orderCues[-length(alldata[[x]]$orderCues)])
  postDSonly=T
  postNSonly=F
  postBoth=F
  
  ImpingeTrials <- sum(sapply(seq(1, length(toTest)), function(z){
    if(postDSonly==T){
      
      test <- ITIprecue[z]>toTest[[z]][1] & ITIprecue[z]<toTest[[z]][2] & preTrialCue[z]==1
    }  
    
    if(postNSonly==T)
    {
      test <- ITIprecue[z]>toTest[[z]][1] & ITIprecue[z]<toTest[[z]][2] & preTrialCue[z]==2
    }
    
    if(postBoth==T){
      test <- ITIprecue[z]>toTest[[z]][1] & ITIprecue[z]<toTest[[z]][2] 
    }
    
    return(test)
  }))
  nTrials <- length(ITIprecue)
  
  prop <- ImpingeTrials/nTrials
  
})


# Figure out ITI latency in trials in which head entries impinge onto ITI period
Impinge <- sapply(seq(1, length(alldata)), function(x){
  entries <- alldata[[x]]$receptacleentries
  exits <- alldata[[x]]$receptacleexits
  pairs <- lapply(seq(1, length(entries)), function(y){
    a <- c(entries[y], exits[y])
  })
  trialStarts <- alldata[[x]]$trialStart
  ITIprecue <- alldata[[x]]$allCues-10
  
  #This tells me the pre-cue periods that start in between two entries. Select the entry-exit pairs that corresponds to those ITI periods
  ITItoEntryIdx <- findInterval(ITIprecue, entries)
  toTest <- pairs[ITItoEntryIdx]
  
  #Select trials based on what cue came on before (modify accordingly). 1=DS, 2=NS.
  preTrialCue <- c(0, alldata[[x]]$orderCues[-length(alldata[[x]]$orderCues)])
  postDSonly=T
  postNSonly=F
  postBoth=F
  
  ImpingeTrialsIdx <- sapply(seq(1, length(toTest)), function(z){
    if(postDSonly==T){
      test <- ITIprecue[z]>toTest[[z]][1] & ITIprecue[z]<toTest[[z]][2] & preTrialCue[z]==1
    }  
    
    if(postNSonly==T)
    {
      test <- ITIprecue[z]>toTest[[z]][1] & ITIprecue[z]<toTest[[z]][2] & preTrialCue[z]==2
    }
    
    if(postBoth==T){
      test <- ITIprecue[z]>toTest[[z]][1] & ITIprecue[z]<toTest[[z]][2] 
    }
    
    return(test)
  })
  
  LatImpingedTrials <- alldata[[x]]$ITIlatency[ImpingeTrialsIdx]
  LatNonImpingedTrials <- alldata[[x]]$ITIlatency[-ImpingeTrialsIdx]
  
  return(data.frame(ITIlatImpTrials=mean(LatImpingedTrials, na.rm=T), ITIlatNONImpTrials=mean(LatNonImpingedTrials, na.rm=T)))
})

toplot <- list(unlist(Impinge[1,]), unlist(Impinge[2,]))

boxplot(x=toplot, whisklty = 0, staplelty = 0, names=c("Mean ITIlat Impinged trials", "ITIlat non-impinged trials"), outline=F)

points(x=rep(1, length(toplot[[1]])), y=toplot[[1]], pch=16)
points(x=rep(2, length(toplot[[2]])), y=toplot[[2]], pch=16)

sapply(seq(1, length(toplot[[1]])), function(y){
  lines(x=c(1, 2), y=c(toplot[[1]][y], toplot[[2]][y]))
})

#PostDS
t.test(x=toplot[[1]], y=toplot[[2]], paired=T) #t = 1.8889, df = 30, p-value = 0.06861
#PostNS
t.test(x=toplot[[1]], y=toplot[[2]], paired=T) #t = -0.23975, df = 26, p-value = 0.8124
#All
t.test(x=toplot[[1]], y=toplot[[2]], paired=T) #t = 1.3213, df = 33, p-value = 0.1955


title("Mean ITI latency in post-DS impinged vs. non-impinged trials")

#I calculated this separately
ImpingePostDS <- Impinge
ImpingePostNS <- Impinge



PostNSentries <- lapply(seq(1, length(alldata)), function(x){
  targetEntries <- alldata[[x]]$PostCSMinusEndEntries
  return(targetEntries)
})

PostNSexits <- lapply(seq(1, length(alldata)), function(x){
  targetExits <- alldata[[x]]$PostCSMinusEndExits
  return(targetExits)
})

PostREWDSentries <- lapply(seq(1, length(alldata)), function(x){
  RewEntries <- !is.na(alldata[[x]]$CSplusresponse)
  return(RewEntries)
})

#### Post DS end entries
plot.new()
plot.window(ylim=c(1,length(alldata)*2), xlim=c(0, 45))

#Plot dots with the entry time
lapply(seq(1, length(PostDSentries)), function(x){
  sapply(seq(1, length(PostDSentries[[x]])), function(y){
    if(PostREWDSentries[[x]][y]==T){col <- "blue"; ypos <- +0.3} else {col <- "black"; ypos <- -0.3}
    
    points(x=PostDSentries[[x]][[y]], y=(rep(x, length(PostDSentries[[x]][[y]])))*2+ypos, pch=16, cex=0.4, col=col)
  })
})

axis(side=1, cex=1.4, at=seq(0, 45, 1))
axis(side=2, cex=1, at=seq(2, length(alldata)*2, by=2), labels=seq(1:length(alldata)), las=2)

mtext(side=1, line=2.5, text="Time from previous cue end (s)", cex=1.5, font=2)
mtext(side=2, line=2.5, text="Session index", cex=1.5, font=2)

legend("topright", col=c("blue", "black"), pch=16, legend=c("Entry after REWARD", "Entry after miss"))


#NS: Plot duration of entries as a function of time from previous cue and whether the previous trial was rewarded or not
plot.new()
plot.window(xlim=c(0, 15), ylim=c(0, 25))

PostDSdurEnt <- lapply(seq(1, length(PostDSentries)), function(x){
  sapply(seq(1, length(PostDSentries[[x]])), function(y){
    entries <- PostDSentries[[x]][[y]]
    exits <- PostDSexits[[x]][[y]]
    
    if(sum(c(!is.na(entries), !is.na(exits)))>1){
      if(exits[1]<entries[1]){exits <- exits[-1]}
      if(last(entries)>last(exits)){entries <- entries[-length(entries)]}
      
      entries <- entries[entries <= 15]
      
      sapply(seq(1, length(entries)), function(z){
        d <- exits[z]-entries[z]
        if(PostREWDSentries[[x]][[y]]==T){color="blue"} else {color="black"}
        points(x=entries[z], y=d, pch=16, col=color, cex=0.5)
        return(d)
      })
    }
  })
})

axis(side=1, at=seq(0, 15, 5), cex=1.4)
axis(side=2, at=seq(0, 26, 2), las=2)

mtext(side=1, line=2.5, text="Time from previous S+ end (s)", cex=1.5, font=2)
mtext(side=2, line=2.5, text="Duration of entry (s)", cex=1.5, font=2)

legend("topright", col=c("blue", "black"), pch=16, legend=c("Entry after REWARD", "Entry after miss"), cex=1.5)


#NS: Plot duration of entries as a function of time from previous cue
plot.new()
plot.window(xlim=c(0, 15), ylim=c(0, 25))

PostNSdurEnt <- lapply(seq(1, length(PostNSentries)), function(x){
  sapply(seq(1, length(PostNSentries[[x]])), function(y){
    entries <- PostNSentries[[x]][[y]]
    exits <- PostNSexits[[x]][[y]]
    
    if(sum(c(!is.na(entries), !is.na(exits)))>1){
      if(exits[1]<entries[1]){exits <- exits[-1]}
      if(last(entries)>last(exits)){entries <- entries[-length(entries)]}
      
      entries <- entries[entries <= 15]
      
      sapply(seq(1, length(entries)), function(z){
        d <- exits[z]-entries[z]
        points(x=entries[z], y=d, pch=16, cex=0.5)
      })
    }
  })
})

axis(side=1, at=seq(0, 15, 5), cex=1.4)
axis(side=2, at=seq(0, 26, 2), las=2)

mtext(side=1, line=2.5, text="Time from previous S- end (s)", cex=1.5, font=2)
mtext(side=2, line=2.5, text="Duration of entry (s)", cex=1.5, font=2)


### HISTOGRAM: PROBABILITY OF POKING HEAD INSIDE RECEPTACLE AFTER s+ (REWARDED)
limXaxis <- 15
binCuts <- seq(0, limXaxis, by=1)
binIdx <- 0:length(binCuts)
POSTREW=T #Whether I want the entriews after rewarded DS (POSTREW=T) or after all DS cues (POSTREW=F)

ProbENTperBin <- sapply(seq(1, length(PostDSentries)), function(x){
  
  
  
  sortedEntriesPerSess <- sort(unlist(PostDSentries[[x]]))
  sortedEntriesPerSess <- sortedEntriesPerSess[sortedEntriesPerSess<limXaxis]
  entPerBinIdx <- findInterval(sortedEntriesPerSess, binCuts)
  
})

allEntriesPerBin <- unlist(ProbENTperBin)
#allEntriesPerBin <- allEntriesPerBin[allEntriesPerBin!=1]

hist(x=allEntriesPerBin, breaks=binIdx, xlab="Time from reward (s)", labels=binCuts, main="Frequency of entries at different times after reward")

legend("topright", legend="n sess=43 (n rats=7)")


### HISTOGRAM: PROBABILITY OF POKING HEAD INSIDE RECEPTACLE AFTER S-
limXaxis <- 15
binCuts <- seq(0, limXaxis, by=1)
binIdx <- 0:length(binCuts)
POSTREW=T #Whether I want the entriews after rewarded DS (POSTREW=T) or after all DS cues (POSTREW=F)

ProbENTperBin <- sapply(seq(1, length(PostNSentries)), function(x){
  
  
  
  sortedEntriesPerSess <- sort(unlist(PostNSentries[[x]]))
  sortedEntriesPerSess <- sortedEntriesPerSess[sortedEntriesPerSess<limXaxis]
  entPerBinIdx <- findInterval(sortedEntriesPerSess, binCuts)
  
})

allEntriesPerBin <- unlist(ProbENTperBin)
#allEntriesPerBin <- allEntriesPerBin[allEntriesPerBin!=1]

hist(x=allEntriesPerBin, breaks=binIdx, xlab="Time from reward (s)", labels=binCuts, main="Frequency of entries at different times after S-", ylim=c(0, 1500))

legend("topright", legend="n sess=43 (n rats=7)")


# HISTOGRAM OF TIME SPENT INSIDE RECEPTACLE 12S AFTER REWARD

POSTREW=T #Whether I want the entries after rewarded DS (POSTREW=T) or after all DS cues (POSTREW=F)

meanDur <- sapply(seq(1, length(PostDSdurEnt)), function(x){
  
  if(POSTREW==T){trialSel <- PostREWDSentries[[x]]} else {trialSel <- 1:length(PostDSdurEnt[[x]])}
  PostDSdurEntSel <- PostDSdurEnt[[x]][trialSel]
  
  nTrials <- length(PostDSdurEntSel)
  sum(unlist(PostDSdurEntSel), na.rm=T)/nTrials
  
})

hist(meanDur)       



#Percentage of trials in which the animal was still in at the time the ITI window started
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 0.5))

ImpingeAll <- list(ImpingePostDS, ImpingePostNS)

boxplot(x=ImpingeAll, whisklty = 0, staplelty = 0, names=c("Post S+ trials", "Post S- trials"), outline=F)
title("Percentage of trials in which rat was INSIDE comp. at the time of 10s pre-cue period")

lapply(seq(1, length(ImpingeAll)), function(x){
  
  points(x=jitter(rep(x, length(ImpingeAll[[x]]))), y=ImpingeAll[[x]], pch=16, cex=0.5)
})











###########################################################
###############################
### I'm going to use the 10s pre-cue period as a baseline. 
### In order to use more than one trial for baseline calculation purposes, I'm going to use +/- n trials around the trial of interest.

cuesOrderAll <- lapply(seq(1, length(rats)), function(k){
  do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
    sessIdx <- idx[[k]][m]
    cuesOrder <- alldata[[sessIdx]]$orderCues
  }))
})

ITIlatMovingAvg <- sapply(seq(1, length(rats)), function(k){
  nTrialWdw <- 5 #Number of trials around the trial of interest that I want to use for the moving window
  nTrials <- length(ITIlatency[[k]])
  sapply(seq(1, nTrials), function(l){
    if(l<(1+nTrialWdw)){
      trlIdx <- 1:((nTrialWdw*2)+1)
    } 
    if(l>(nTrials-nTrialWdw)){
      trlIdx <- (nTrials-(nTrialWdw*2)):nTrials
    }
    if(l>nTrialWdw & l<=(nTrials-nTrialWdw)){
      trlIdx <- (l-nTrialWdw):(l+nTrialWdw)
    }
    
    avg <- mean(ITIlatency[[k]][trlIdx], na.rm=T)
  })
})


ITIlatMovingSD <- sapply(seq(1, length(rats)), function(k){
  nTrialWdw <- 5 #Number of trials around the trial of interest that I want to use for the moving window
  nTrials <- length(ITIlatency[[k]])
  sapply(seq(1, nTrials), function(l){
    if(l<(1+nTrialWdw)){
      trlIdx <- 1:((nTrialWdw*2)+1)
    } 
    if(l>(nTrials-nTrialWdw)){
      trlIdx <- (nTrials-(nTrialWdw*2)):nTrials
    }
    if(l>nTrialWdw & l<=(nTrials-nTrialWdw)){
      trlIdx <- (l-nTrialWdw):(l+nTrialWdw)
    }
    
    sd <- sd(ITIlatency[[k]][trlIdx], na.rm=T)
  })
})




#Compare raw and moving average ITI latency
plot.new()
plot.window(xlim=c(1, length(ITIlatency[[k]])), ylim=c(0, 90))

lapply(seq(1, length(rats)), function(k){
  ylevel <- (k-1)*15
  lines(x=1:length(ITIlatency[[k]]), y=ylevel+ITIlatency[[k]], col="black")
  lines(x=1:length(ITIlatCorrected[[k]]), y=ylevel+ITIlatCorrected[[k]], col="red", lwd=2)
  axis(side=2, at=seq(ylevel, ylevel+10, by=5), labels=seq(0, 10, 5), las=2)
  })

axis(side=1, at=seq(0, 560, by=80))
mtext(side=2, line=2.5, text="ITI latency (s)", cex=1.5)
mtext(side=1, line=2.5, text="Trials", cex=1.5)

legend("topright", lty=1, lwd=2, col="red", legend="Moving average (+/- 5 trials)")
title("ITI latency per animal (raw and moving average)")

####
# DS and NS performance index using moving average for ITI latency

DStaskAccMA <- lapply(seq(1, length(ITIlatMovingAvg)), function(k){
  DS_ITIlatMA <- ITIlatMovingAvg[[k]][cuesOrderAll[[k]]==1]
  DStaskAccMA <- DS_ITIlatMA-DSlatency[[k]]
  })

NStaskAccMA <- lapply(seq(1, length(ITIlatMovingAvg)), function(k){
  NS_ITIlatMA <- ITIlatMovingAvg[[k]][cuesOrderAll[[k]]==2]
  NStaskAccMA <- NS_ITIlatMA-NSlatency[[k]]
})
  
#Compare raw and moving average DS task accuracy (performance index or ITI latency - DS latency)
plot.new()
plot.window(xlim=c(1, length(DStaskAccMA[[k]])), ylim=c(0, 90))

lapply(seq(1, length(rats)), function(k){
  ylevel <- ((k-1)*15)+5
  lines(x=1:length(DStaskAcc[[k]]), y=ylevel+DStaskAcc[[k]], col="black")
  lines(x=1:length(DStaskAccMA[[k]]), y=ylevel+DStaskAccMA[[k]], col="red", lwd=1.5)
  axis(side=2, at=seq(ylevel-5, ylevel+5, by=5), labels=seq(-5, 5, 5), las=2)
})

axis(side=1, at=seq(0, 280, by=40))
mtext(side=2, line=2.5, text="Performance index (s)", cex=1.5)
mtext(side=1, line=2.5, text="Trials", cex=1.5)

legend("topright", lty=1, lwd=2, col="red", legend="Perf.Idx. calculated using moving average of ITI latency (+/- 5 trials)")
title("Performance index per animal (Using per trial or moving average ITI latency")

#Compare raw and corrected DS task accuracy (performance index or ITI latency - DS latency)
plot.new()
plot.window(xlim=c(1, length(DStaskAccCorrected[[k]])), ylim=c(0, 90))

lapply(seq(1, length(rats)), function(k){
  ylevel <- ((k-1)*15)+5
  lines(x=1:length(DStaskAcc[[k]]), y=ylevel+DStaskAcc[[k]], col="black")
  lines(x=1:length(DStaskAccCorrected[[k]]), y=ylevel+DStaskAccCorrected[[k]], col="red", lwd=1.5)
  axis(side=2, at=seq(ylevel-5, ylevel+5, by=5), labels=seq(-5, 5, 5), las=2)
})

axis(side=1, at=seq(0, 280, by=40))
mtext(side=2, line=2.5, text="Performance index (s)", cex=1.5)
mtext(side=1, line=2.5, text="Trials", cex=1.5)



#Compare raw and moving average NS task accuracy (performance index or ITI latency - NS latency)
plot.new()
plot.window(xlim=c(1, length(NStaskAccMA[[k]])), ylim=c(0, 90))

lapply(seq(1, length(rats)), function(k){
  ylevel <- ((k-1)*15)+5
  lines(x=1:length(NStaskAcc[[k]]), y=ylevel+NStaskAcc[[k]], col="black")
  lines(x=1:length(NStaskAccMA[[k]]), y=ylevel+NStaskAccMA[[k]], col="red", lwd=1.5)
  axis(side=2, at=seq(ylevel-5, ylevel+5, by=5), labels=seq(-5, 5, 5), las=2)
})

axis(side=1, at=seq(0, 280, by=40))
mtext(side=2, line=2.5, text="Performance index (s)", cex=1.5)
mtext(side=1, line=2.5, text="Trials", cex=1.5)

legend("topright", lty=1, lwd=2, col="red", legend="Perf.Idx. calculated using moving average of ITI latency (+/- 5 trials)")
title("Performance index per animal (Using per trial or moving average ITI latency")

#Compare raw and corrected NS task accuracy (performance index or ITI latency - NS latency)
plot.new()
plot.window(xlim=c(1, length(NStaskAccCorrected[[k]])), ylim=c(0, 90))

lapply(seq(1, length(rats)), function(k){
  ylevel <- ((k-1)*15)+5
  lines(x=1:length(NStaskAcc[[k]]), y=ylevel+NStaskAcc[[k]], col="black")
  lines(x=1:length(NStaskAccCorrected[[k]]), y=ylevel+NStaskAccCorrected[[k]], col="red", lwd=1.5)
  axis(side=2, at=seq(ylevel-5, ylevel+5, by=5), labels=seq(-5, 5, 5), las=2)
})

axis(side=1, at=seq(0, 280, by=40))
mtext(side=2, line=2.5, text="Performance index (s)", cex=1.5)
mtext(side=1, line=2.5, text="Trials", cex=1.5)


### COMPARE ITI LATENCY RAW AND USING MOVING AVERAGE ON ALL TRIALS
plot.new()
toplot <- list(unlist(ITIlatency), unlist(ITIlatMovingAvg))
boxplot(x=toplot, whisklty = 0, staplelty = 0, names=c("ITI lat raw", "ITI lat moving average"), outline=F)

jit <- jitter(rep(1, length(toplot[[1]])), factor=10)
points(x=jit, y=toplot[[1]], pch=16, cex=0.3, col="red")
jit <- jitter(rep(2, length(toplot[[1]])), factor=5)
points(x=jit, y=toplot[[2]], pch=16, cex=0.3, col="red")



### COMPARE ITI LATENCY RAW AND USING MOVING AVERAGE ON impinged trials only
plot.new()
toplot <- list(unlist(ITIlatency), unlist(ITIlatCorrected))
boxplot(x=toplot, whisklty = 0, staplelty = 0, names=c("ITI lat raw", "ITI lat corrected"), outline=F)

jit <- jitter(rep(1, length(toplot[[1]])), factor=10)
points(x=jit, y=toplot[[1]], pch=16, cex=0.3, col="red")
jit <- jitter(rep(2, length(toplot[[1]])), factor=5)
points(x=jit, y=toplot[[2]], pch=16, cex=0.3, col="red")

mtext(side=2, line=2.5, text="ITI latency (s)", cex=1.5)


# How many trials had to be corrected?
allCorrected <- unlist(ImpingeIdx)
sum(allCorrected)/length(allCorrected) #0.09391026

# ITI latency in trials before vs after correction
colorIndex <- brewer.pal(length(rats), "Accent")

plot.new()
plot.window(xlim=c(1, 2), ylim=c(0, 10))
lapply(seq(1, length(rats)), function(k){
  preCorrection <- ITIlatency[[k]][ImpingeIdx[[k]]]
  postCorrection <- ITIlatCorrected[[k]][ImpingeIdx[[k]]]
  colpick <- colorIndex[k]
  sapply(seq(1, length(preCorrection)), function(l){
    lines(x=c(1, 2), y=c(preCorrection[l], postCorrection[l]), col=colpick)
  })

})

axis(side=1, at=c(1, 2), labels=c("Pre-correction", "Post-correction"), cex=1.4, font=2)
axis(side=2, las=2)
mtext(side=2, line=2.5, text="ITI latency(s)", cex=1.4, font=2)
