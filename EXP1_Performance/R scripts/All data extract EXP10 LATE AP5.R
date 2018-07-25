

#############################################################
### EXPERIMENT 10: LATE AP5 VS VEH TEST                  ###
#############################################################

Exp10folder <- "E:/Dropbox/DISSERTATION/Experiment 10/"

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



### LATE AP5 #################################################################
#Define folders
subTestFolder <- paste(Exp10folder, "Late AP5 test/", sep="")
datafolder <- paste(subTestFolder, "MedPC/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
dataForRCumulativeEarlyAP5 <- dataForRCumulative
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")

### LATE VEH #################################################################
#Define folders
subTestFolder <- paste(Exp10folder, "Late VEH test/", sep="")
datafolder <- paste(subTestFolder, "MedPC/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
dataForRCumulativeEarlyVEH <- dataForRCumulative
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")

#Create objects
MedPCextract(funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative, cuelength=10)

# Load important behavior-related objects
files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); for(i in 1:length(filesCum)){load(filesCum[[i]])}

csacqidx$session <- rep(1, nrow(csacqidx))

CPdata <- CPextract(GallCrit=1.3, minSlope=0, adjDrop=-0.5, CPfuncFolder=CPfuncFolder, idx=idx, CPGraphFolder=CPGraphFolder, dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative)
CPdata$CP <- c(rep(0, nrow(CPdata)))

allNeuronsDS <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=F, side="both")
allNeuronsDSLateAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=T, allResults=T, side="both")
allNeuronsDSLateAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=1800, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=T, allResults=T, side="both")
allNeuronsNSLateAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSLateAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=1800, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")

allNeuronsDSLateVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=T, allResults=T, side="both")
allNeuronsDSLateVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=1800, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=T, allResults=T, side="both")
allNeuronsNSLateVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSLateVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=1800, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")

masterDF <- FRbyNEURONbyBINcue(neudata=allNeuronsDSLateAP5PreInf, funcdirect=funcdirect, dataForRdir=dataForRdir, format="freq", BLduration=psthmin, cueExcOnly=F, cueKind=1)

### GIVE THESE OBJECTS A UNIQUE NAME

## AP5 SIDE
csacqidxLateAP5 <- csacqidx
alldataLateAP5 <- alldata
ratsLateAP5 <- rats
idxLateAP5 <- idx
cumDataLateAP5 <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)
CPdataLateAP5 <- CPdata
masterDFLateAP5 <- masterDF

## VEH SIDE
csacqidxLateVEH <- csacqidx
alldataLateVEH <- alldata
ratsLateVEH <- rats
idxLateVEH <- idx
cumDataLateVEH <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)
CPdataLateVEH <- CPdata
masterDFLateVEH <- masterDF


#Performance throughout the session per bin
dfAllPerBinLateAP5 <- PerformanceFromInf(timeBinSize=600, maxSess=9000, makeGraph=F, comp=c("VEH", "AP5"), colindx=c("blue", "red"), moreInfo = "Late", 
                                          imgFormat="pdf", graphFolder=behGraphFolder, alldata=list(alldataLateAP5))
dfAllPerBinLateVEH <- PerformanceFromInf(timeBinSize=600, maxSess=9000, makeGraph=F, comp=c("VEH", "AP5"), colindx=c("blue", "red"), moreInfo = "Late", 
                                          imgFormat="pdf", graphFolder=behGraphFolder, alldata=list(alldataLateVEH))


#Behavior per bin
plotBehPerBin(data=list(dfAllPerBinLateVEH, dfAllPerBinLateAP5), percOfBL=T, binDur=600, infStart=1800, infDur=12*60, imgFormat="pdf", graphFolder=behGraphFolder, toplot=c("S- specificity"), colGroups=c("blue", "red"), lty=c(1, 2), groups = c("VEH", "AP5"), exp="Late")
        
### PLOTTING IMPORTANT GRAPHS
plotFRandINF(experiment="Exp10 Late AP5", masterDF=masterDF, infTime=1800, infDur=12*60, graphFolder=MixedGraphFolder, timeBinSize=300, postCueWdw=400, dataProcess="raw", correctOnly=FALSE, color="black", capped=F, capValue = c(1, 7000), yAxMinZ = -1, yAxMaxZ = 10, yAxMaxRaw = 20, neudata=allNeuronsDS)

#Load cumulative files of both tests if necessary
PerformanceFromInf(timeBinSize=600, maxSess=9000, comp=c("VEH", "AP5"), colindx=c("blue", "red"), moreInfo = "Late", 
                   imgFormat="pdf", graphFolder=behGraphFolder, alldata=list(alldataLateVEH, alldataLateAP5))
        
     

#PSTH
#AP5
psthInf(formatDat="Zscores", group="AP5", expName="Late", errShade=F, ymax=14, comp=c("Pre-infusion", "Post-infusion"), graphFolder=MixedGraphFolder, colindx=c("black", "red"), infTime=1800, infDur=12*60, psthmin=1.5, psthmax=1.5, binw=50, neudata=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf))
psthInf(formatDat="Zscores", group="AP5", expName="Late", errShade=F, ymax=14, comp=c("Pre-infusion", "Post-infusion"), graphFolder=MixedGraphFolder, colindx=c("black", "red"), infTime=1800, infDur=12*60, psthmin=1.5, psthmax=1.5, binw=50, neudata=list(allNeuronsNSLateAP5PreInf, allNeuronsNSLateAP5PostInf))

#VEH
psthInf(formatDat="Zscores", group="VEH", expName="Late", errShade=F, ymax=14, comp=c("Pre-infusion", "Post-infusion"), graphFolder=MixedGraphFolder, colindx=c("black", "blue"), infTime=1800, infDur=12*60, psthmin=1.5, psthmax=1.5, binw=50, neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf))
psthInf(formatDat="Zscores", group="VEH", expName="Late", errShade=F, ymax=14, comp=c("Pre-infusion", "Post-infusion"), graphFolder=MixedGraphFolder, colindx=c("black", "blue"), infTime=1800, infDur=12*60, psthmin=1.5, psthmax=1.5, binw=50, neudata=list(allNeuronsNSLateVEHPreInf, allNeuronsNSLateVEHPostInf))
