

#############################################################
### EXPERIMENT 10: LEARNING  (Down and not down)          ###
#############################################################

 Exp10folder <- "E:/Dropbox/DISSERTATION/Learning Exp 10/"

### LOAD FUNCTIONS
funcdirect <- "E:/Dropbox/DISSERTATION/R functions/"
CPfuncFolder <- paste(funcdirect, 'Change_Point-master/', sep="")

#Load functions
load(file=paste(funcdirect, "MedPCextract.r", sep=""))
load(file=paste(funcdirect, "mpcextract_blockSingle.Rfunc", sep=""))
load(file=paste(funcdirect, "CPextract.r", sep=""))
load(file=paste(funcdirect, "neuralhist.r", sep=""))
load(file=paste(funcdirect, "FRbyNEURONbyBINcue.r", sep=""))
load(file=paste(funcdirect, "plotFRandCP.r", sep=""))
load(file=paste(funcdirect, "errBars.r", sep=""))
load(file=paste(funcdirect, "errCloud.r", sep=""))
load(file=paste(funcdirect, "psthInf.r", sep=""))
#load(file=paste(funcdirect, "plotFRandCPhistogram.r", sep=""))



### Late #################################################################
#Define folders
subTestFolder <- paste(Exp10folder, "Late/", sep="")
datafolder <- paste(subTestFolder, "MedPC/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
IndcumPerfGraphFolder <- paste(behGraphFolder, "Cumulative/", sep="")
PerfRelToCPFolder <- paste(behGraphFolder, "Beh rel to CP/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")

# ### Early VEH ##############################################################
# #Define folders
# subTestFolder <- paste(Exp10folder, "Early VEH/", sep="")
# datafolder <- paste(subTestFolder, "MedPC/", sep="")
# dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
# dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
# dataForRCumulativeEarlyVEH <- dataForRCumulative
# behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
# IndcumPerfGraphFolder <- paste(behGraphFolder, "Cumulative/", sep="")
# PerfRelToCPFolder <- paste(behGraphFolder, "Beh rel to CP/", sep="")
# MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
# CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
# NEXfiles <- paste(subTestFolder, "NEX files/", sep="")


### ALL ####################################################################
# subTestFolder <- paste(Exp10folder, "All/", sep="")
# datafolder <- paste(subTestFolder, "MedPC/", sep="")
# dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
# dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
# behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
# IndcumPerfGraphFolder <- paste(behGraphFolder, "Cumulative/", sep="")
# PerfRelToCPFolder <- paste(behGraphFolder, "Beh rel to CP/", sep="")
# MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
# CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
# NEXfiles <- paste(subTestFolder, "NEX files/", sep="")

#Create objects
MedPCextract(funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative, cuelength=10)

# Load important behavior-related objects
files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); for(i in 1:length(filesCum)){load(filesCum[[i]])}

CPdata <- CPextract(GallCrit=1.3, minSlope=0, adjDrop=-0.5, CPfuncFolder=CPfuncFolder, idx=idx, CPGraphFolder=CPGraphFolder, dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative)


#allNeuronsDSDown <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=T, allResults=T, side="both")
#allNeuronsDSNotDown <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="both")
#allNeuronsDSAll <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="both")

#masterDF <- FRbyNEURONbyBINcue(neudata=allNeuronsDSDown, funcdirect=funcdirect, dataForRdir=dataForRdir, format="freq", BLduration=psthmin)

### GIVE THESE OBJECTS A UNIQUE NAME

## LATE
csacqidxLate <- csacqidx
alldataLate <- alldata
ratsLate <- rats
idxLate <- idx
cumDataLate <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)
CPdataLate <- CPdata
masterDFLate <- masterDF
day7DataLate <- filter(csacqidx, session==7)        
       


# ## ALL
# csacqidxAll <- csacqidx
# alldataAll <- alldata
# ratsAll <- rats
# idxAll <- idx
# cumDataAll <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)
# CPdataAll <- CPdata
# masterDFAll <- masterDF



### PLOTTING IMPORTANT GRAPHS

### DAY 7 PERFORMANCE
plot.new()
plot.window(xlim=c(0, 1), ylim=c(-1, 5))

DStaskAccmeanAllD7 <- mean(day7DataLate$CSplusTaskAcc)
DStaskAccSEMAllD7 <- sd(day7DataLate$CSplusTaskAcc)/sqrt(4)

NStaskAccmeanAllD7 <- mean(day7DataLate$CSminusTaskAcc)
NStaskAccSEMAllD7 <- sd(day7DataLate$CSminusTaskAcc)/sqrt(4)

points(x=0.1, y=DStaskAccmeanAllD7, pch=19, cex=2)
errBars(x=0.1, y=DStaskAccmeanAllD7, err=DStaskAccSEMAllD7)

points(x=0.2, y=NStaskAccmeanAllD7, pch=19, cex=2, col="darkred")
errBars(x=0.2, y=NStaskAccmeanAllD7, err=NStaskAccSEMAllD7)

axis(side=1, at=c(0.1, 0.2), labels=c("S+", "S-"), cex.axis=1.7, font=2)
axis(side=2, las=2, cex.axis=1.8)
mtext(side=2, text="Performance index (s)", line=2.5, cex=1.5, font=2)
abline(h=0, lty=3)

# 
# points(x=0.3, y=mean(day7DataEarlyAP5$CSplusTaskAcc, na.rm=T), pch=19, cex=2)
# errBars(x=0.3, y=mean(day7DataEarlyAP5$CSplusTaskAcc, na.rm=T), err=sd(day7DataEarlyAP5$CSplusTaskAcc, na.rm=T)/sqrt(nrow(day7DataEarlyAP5)))
# points(x=0.3, y=mean(day7DataEarlyVEH$CSplusTaskAcc, na.rm=T), pch=19, cex=2)
# errBars(x=0.3, y=mean(day7DataEarlyVEH$CSplusTaskAcc, na.rm=T), err=sd(day7DataEarlyVEH$CSplusTaskAcc, na.rm=T)/sqrt(nrow(day7DataEarlyVEH)))


### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point DOWN
plotFRandCP(experiment="TrainingTask2Down", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="raw", correctOnly=FALSE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSDown)
plotFRandCP(experiment="TrainingTask2Down", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="Zscores", correctOnly=FALSE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSDown)
plotFRandCP(experiment="TrainingTask2Down", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="raw", correctOnly=TRUE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSDown)
plotFRandCP(experiment="TrainingTask2Down", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="Zscores", correctOnly=TRUE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSDown)

### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point NOT DOWN
plotFRandCP(experiment="TrainingTask2NotDown", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="raw", correctOnly=FALSE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSNotDown)
plotFRandCP(experiment="TrainingTask2NotDown", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="Zscores", correctOnly=FALSE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSNotDown)
plotFRandCP(experiment="TrainingTask2NotDown", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="raw", correctOnly=TRUE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSNotDown)
plotFRandCP(experiment="TrainingTask2NotDown", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="Zscores", correctOnly=TRUE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSNotDown)

### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point ALL
plotFRandCP(experiment="TrainingTask2All", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="raw", correctOnly=FALSE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSAll)
plotFRandCP(experiment="TrainingTask2All", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="Zscores", correctOnly=FALSE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSAll)
plotFRandCP(experiment="TrainingTask2All", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="raw", correctOnly=TRUE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSAll)
plotFRandCP(experiment="TrainingTask2All", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="Zscores", correctOnly=TRUE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 7, neudata=allNeuronsDSAll)


### PLOT FR AROUND THE CUE (PTSH) as a function of distance to change point
plotFRandCPhistogram(experiment="EXTINCTION", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=15, dataProcess="Zscores", correctOnly=TRUE, color="black", capped=T, capValue = c(-30, 30), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, imgFormat="pdf")


#### BEHAVIOR

#CUMULATIVE INDIVIDUAL PERFORMANCE

cumulativeIndGraphs(numSess=6, sessLines=T, dataForRCumulative, dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf")

PerformanceFromCP(relTrialMin=-100, relTrialMax=100, idx=idx, trialBinSize=5, typegraph="both", imgFormat="pdf", dataForRCumulative, dataForRdir, graphFolder=PerfRelToCPFolder, CPdata, rats=rats, alldata=alldata)


#COMPARE WITH OTHER GROUPS
compareCPs(data=list(CPdataExp4, CPdataExp10tr, CPdataHYBunilAP5, CPdataHYBbilAP5) , imgFormat="png", expNames=c("Task1", "Task2", "UnilAP5", "EXTINCTION"), colindx=c("gray30", "gray45", "gray60", "gray75"), graphFolder=behGraphFolder, minSess=5, graphs=c("Change point", "Non learners", "Dynamic Interval", "Asymptote"))


















### INFUSION DATA
#plotFRandINF(experiment="Exp10 Down", masterDF=masterDF, infTime=1800, infDur=12*60, graphFolder=MixedGraphFolder, timeBinSize=300, postCueWdw=400, dataProcess="raw", correctOnly=FALSE, color="black", capped=F, capValue = c(1, 7000), yAxMinZ = -1, yAxMaxZ = 10, yAxMaxRaw = 20, neudata=allNeuronsDS)


#Load cumulative files of both tests if necessary
#PerformanceFromInf(timeBinSize=600, maxSess=9000, comp=c("VEH", "AP5"), colindx=c("blue", "red"), moreInfo = "Late", 
                   #imgFormat="pdf", graphFolder=behGraphFolder, alldata=list(alldataLateVEH, alldataLateAP5))
        

#PSTH
#AP5
psthInf(formatDat="Zscores", group="AP5", expName="Late", errShade = F, comp=c("Pre-infusion", "Post-infusion"), graphFolder=MixedGraphFolder, colindx=c("black", "red"), infTime=1800, infDur=12*60, psthmin=1.5, psthmax=5, binw=50, neudata=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf))

#VEH
#psthInf(formatDat="Zscores", group="VEH", expName="Late", comp=c("Pre-infusion", "Post-infusion"), graphFolder=MixedGraphFolder, colindx=c("black", "blue"), infTime=1800, infDur=12*60, psthmin=1.5, psthmax=5, binw=50, neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf))

        