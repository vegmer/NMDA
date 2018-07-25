

#############################################################
### EXPERIMENT 1A: EARLY AP5 VS VEH TEST                  ###
#############################################################

Exp1folder <- "E:/Dropbox/NMDA/EXP1_Performance/"

##########################
##########################
### LOAD FUNCTIONS     ###
##########################
##########################

funcdirect <- "E:/Dropbox/NMDA/R functions/"
CPfuncFolder <- paste(funcdirect, 'Change_Point-master/', sep="")

#Load functions
load(file=paste(funcdirect, "MedPCextract.r", sep=""))
load(file=paste(funcdirect, "mpcextract_blockSingle.Rfunc", sep=""))
load(file=paste(funcdirect, "CPextract.r", sep=""))
load(file=paste(funcdirect, "neuralhist.r", sep=""))
load(file=paste(funcdirect, "FRbyNEURONbyBINcue.r", sep=""))
load(file=paste(funcdirect, "errBars.r", sep=""))
load(file=paste(funcdirect, "errCloud.r", sep=""))
load(file=paste(funcdirect, "psthInf.r", sep=""))


###########################
###########################
### DEFINE FOLDERS      ###
###########################
###########################


# Define folders for one group *OR* the other before running the rest of the code (not both because then you'll just rewrite the folders you defined for the first group)

### EARLY AP5 #################################################################
subTestFolder <- paste(Exp1folder, "Early AP5/", sep="")
datafolder <- paste(subTestFolder, "MedPC files/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
dataForRCumulativeEarlyAP5 <- dataForRCumulative
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")

### EARLY VEH #################################################################
subTestFolder <- paste(Exp1folder, "Early VEH/", sep="")
datafolder <- paste(subTestFolder, "MedPC files/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
dataForRCumulativeEarlyVEH <- dataForRCumulative
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")


################################
################################
### CREATE IMPORTANT OBJECTS ###
################################
################################

# Extract behavioral data from MedPC files. This function saves the generated objects in the "dataForRdir". You have to load them (see next line of code) to bring them to your environment.
# This will give you a few error messages if, in any file, the first cue comes on after 5s of session onset. Ignore it, it just assigns NA to that trial, which is what you want.
# The parameter 'consumeRewWdw' is just the segment of the ITI that we discard (for ITI latency calculations) bc we assume that, if the animal got a reward on the previous trial, he might still be consuming the reward.
MedPCextract(cuelength=10, consumeRewWdw=4.9, funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative)

# Load the behavior-related objects that you generated with the previous function
files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); for(i in 1:length(filesCum)){load(filesCum[[i]])}

#The main objects that we loaded are 'alldata' (detailed data by session by animal) and 'csacqidx' (an index of all the files). Name all sessions on csacqidx the same (i.e. '1')
csacqidx$session <- rep(1, nrow(csacqidx))

# Extract neuronal data from NEX files. 

#AP5 test: data aligned to DS and NS onset BEFORE and AFTER the infusion.
allNeuronsDSEarlyAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsDSEarlyAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=2520, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=2520, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")


#VEH test: data aligned to DS and NS onset BEFORE and AFTER the infusion.
allNeuronsDSEarlyVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsDSEarlyVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=1800, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=2520, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")


### GIVE THESE OBJECTS A UNIQUE NAME

## AP5 SIDE
csacqidxEarlyAP5 <- csacqidx
alldataEarlyAP5 <- alldata
ratsEarlyAP5 <- rats
idxEarlyAP5 <- idx
cumDataEarlyAP5 <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)

## VEH SIDE
csacqidxEarlyVEH <- csacqidx
alldataEarlyVEH <- alldata
ratsEarlyVEH <- rats
idxEarlyVEH <- idx
cumDataEarlyVEH <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)


######################################################
######################################################
### PLOT GRAPHS                                    ###
######################################################
######################################################


###################
### 1. BEHAVIOR ###
###################

### 1.1. LEARNING CURVE














plotFRandINF(experiment="Exp10 Late AP5", masterDF=masterDF, infTime=1800, infDur=12*60, graphFolder=MixedGraphFolder, timeBinSize=300, postCueWdw=400, dataProcess="raw", correctOnly=FALSE, color="black", capped=F, capValue = c(1, 7000), yAxMinZ = -1, yAxMaxZ = 10, yAxMaxRaw = 20, neudata=allNeuronsDS)


#Performance throughout the session
dfAllPerBinEarlyAP5 <- PerformanceFromInf(timeBinSize=600, maxSess=9000, makeGraph=F, comp=c("VEH", "AP5"), colindx=c("blue", "red"), moreInfo = "Early", 
                   imgFormat="pdf", graphFolder=behGraphFolder, alldata=list(alldataEarlyAP5))
dfAllPerBinEarlyVEH <- PerformanceFromInf(timeBinSize=600, maxSess=9000, makeGraph=F, comp=c("VEH", "AP5"), colindx=c("blue", "red"), moreInfo = "Early", 
                                  imgFormat="pdf", graphFolder=behGraphFolder, alldata=list(alldataEarlyVEH))


#Behavior per bin
plotBehPerBin(data=list(dfAllPerBinEarlyVEH, dfAllPerBinEarlyAP5), percOfBL=T, binDur=600, infStart=1800, infDur=12*60, imgFormat="pdf", graphFolder=behGraphFolder, toplot=c("S+ specificity"), colGroups=c("blue", "red"), lty=c(1, 2), groups = c("VEH", "AP5"), exp="Early")

# PLOT S+ AND S- PER BIN
#PSTH pre and post infusion
#AP5
psthInf(formatDat="Zscores", group="AP5", expName = "Early", errShade=F, ymax=14, graphFolder=MixedGraphFolder, colindx=c("black", "red"), infTime=1800, infDur=12*60, psthmin=1.5, psthmax=1.5, binw=50, neudata=list(allNeuronsNSEarlyAP5PreInf, allNeuronsNSEarlyAP5PostInf))       

#VEH
psthInf(formatDat="Zscores", group="VEH", expName = "Early", errShade=F, ymax=14, graphFolder=MixedGraphFolder, colindx=c("black", "blue"), infTime=1800, infDur=12*60, psthmin=1.5, psthmax=1.5, binw=50, neudata=list(allNeuronsNSEarlyVEHPreInf, allNeuronsNSEarlyVEHPostInf))
