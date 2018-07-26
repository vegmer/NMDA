

#############################################################
### EXPERIMENT 1A: EARLY AP5 VS VEH TEST                  ###
#############################################################

library(matrixStats)

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
load(file=paste(funcdirect, "BinIndexCalculator.R", sep=""))


###########################
###########################
### DEFINE FOLDERS      ###
###########################
###########################


# Define folders for one group *OR* the other before running the rest of the code (not both because then you'll just rewrite the folders you defined for the first group)

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


################################
################################
### CREATE IMPORTANT OBJECTS ###
################################
################################

# Extract behavioral data from MedPC files. This function saves the generated objects in the "dataForRdir". You have to load them (see next line of code) to bring them to your environment.
# This will give you a few error messages if, in any file, the first cue comes on after 5s of session onset. Ignore it, it just assigns NA to that trial, which is what you want.
# The parameter 'consumeRewWdw' is just the segment of the ITI that we discard (for ITI latency calculations) bc we assume that, if the animal got a reward on the previous trial, he might still be consuming the reward.
MedPCextract(cuelength=10, consumeRewWdw=4.9, funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative)

# Load the behavior-related objects that you generated with the previous function. The main objects that we loaded are 'alldata' (detailed data by session by animal) and 'csacqidx' (an index of all the files). Name all sessions on csacqidx the same (i.e. '1')
files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); for(i in 1:length(filesCum)){load(filesCum[[i]])}

#This function will create the following objects: DSbinIdx, NSbinIdx and AllCueBinIdx. These are indexes indicating, for each rat and each kind of event, to what bin the events belong. The were generated and saved in the dataForRCumulative folder, so I need to load them
binsize <- 600
BinIndexCalculator(data=alldata, binsize=binsize, sessLength = 9000); filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); for(i in 1:length(filesCum)){load(filesCum[[i]])}

# Create an object with data per bin for each one of our behavioral parameters
# Response ratio:
minBinNo <- min(sapply(DSbinIdx, max))

DSrespRatioByBin <- lapply(seq(1, length(DSrespAll)), function(x){
        sapply(seq(1, minBinNo), function(y){
                DSinBin <- DSrespAll[[x]][DSbinIdx[[x]]==y]
                DSrespRatio <- sum(DSinBin)/length(DSinBin)
        })
})

NSrespRatioByBin <- lapply(seq(1, length(NSrespAll)), function(x){
        sapply(seq(1, minBinNo), function(y){
                NSinBin <- NSrespAll[[x]][NSbinIdx[[x]]==y]
                NSrespRatio <- sum(NSinBin)/length(NSinBin)
        })
})


# Latency:
DSlatencyByBin <- lapply(seq(1, length(DSlatency)), function(x){
        sapply(seq(1, minBinNo), function(y){
                DSinBin <- DSlatency[[x]][DSbinIdx[[x]]==y]
                DSlatencyByBin <- mean(DSinBin, na.rm=T)
        })
})

NSlatencyByBin <- lapply(seq(1, length(NSlatency)), function(x){
        sapply(seq(1, minBinNo), function(y){
                NSinBin <- NSlatency[[x]][NSbinIdx[[x]]==y]
                NSlatencyByBin <- mean(NSinBin, na.rm=T)
        })
})

# Task Accuracy
DStaskAccByBin <- lapply(seq(1, length(DStaskAcc)), function(x){
        sapply(seq(1, minBinNo), function(y){
                DSinBin <- DStaskAcc[[x]][DSbinIdx[[x]]==y]
                DStaskAccByBin <- mean(DSinBin, na.rm=T)
        })
})

NStaskAccByBin <- lapply(seq(1, length(NStaskAcc)), function(x){
        sapply(seq(1, minBinNo), function(y){
                NSinBin <- NStaskAcc[[x]][NSbinIdx[[x]]==y]
                NStaskAccByBin <- mean(NSinBin, na.rm=T)
        })
})

# ITI latency
ITIlatByBin <- lapply(seq(1, length(ITIlatency)), function(x){
        sapply(seq(1, minBinNo), function(y){
                ITIlatInBin <- ITIlatency[[x]][AllCueBinIdx[[x]]==y]
                ITIlatByBin <- mean(ITIlatInBin, na.rm=T)
        })
})



# Extract neuronal data from NEX files. 

#VEH test: data aligned to DS and NS onset BEFORE and AFTER the infusion.
allNeuronsDSEarlyVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsDSEarlyVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=1800, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=2520, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")

#AP5 test: data aligned to DS and NS onset BEFORE and AFTER the infusion.
allNeuronsDSEarlyAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsDSEarlyAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=2520, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=2520, endt=3600, binw=50, psthmin=1.5, psthmax=10, cueexonly=F, allResults=T, side="both")


### GIVE THESE OBJECTS A UNIQUE NAME

## VEH SIDE
csacqidxEarlyVEH <- csacqidx
alldataEarlyVEH <- alldata
ratsEarlyVEH <- rats
idxEarlyVEH <- idx
cumDataEarlyVEH <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)
byBinDataEarlyVEH <- list(DSrespRatioByBin, NSrespRatioByBin, DSlatencyByBin, NSlatencyByBin, DStaskAccByBin, NStaskAccByBin, ITIlatByBin)

## AP5 SIDE
csacqidxEarlyAP5 <- csacqidx
alldataEarlyAP5 <- alldata
ratsEarlyAP5 <- rats
idxEarlyAP5 <- idx
cumDataEarlyAP5 <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)
byBinDataEarlyAP5 <- list(DSrespRatioByBin, NSrespRatioByBin, DSlatencyByBin, NSlatencyByBin, DStaskAccByBin, NStaskAccByBin, ITIlatByBin)



######################################################
######################################################
### PLOT GRAPHS                                    ###
######################################################
######################################################


###################
### 1. BEHAVIOR ###
###################

# Let's create objects to help us select the bins of interest for the pre and the post
# The infusion took place after 30min and it lasted 12min. I'm going to use the 30min before the infusion as baseline and the 30min after the infusion as the post.
PreInfLength <- 30*60            #In sec
PostInfStart <- (30*60)+12*60    #In sec
PostInfEnd <- PostInfStart+30*60 #In sec

BLbinIndex <- (1:minBinNo)[1:(PreInfLength/binsize)]
PostInfBinIndex <- (1:minBinNo)[ceiling(PostInfStart/binsize):(PostInfEnd/binsize)]

#Function for plotting lines more easily. I just need to adjust the data I feed the function, the color and the points
plotPrePostLines <- function(data, color, pch, scores, jitter=0){
        mat <- do.call("rbind", data) #Create matrix in which rows are different rats and columns are bins
        
        if(scores=="absolute"){
                BLmean <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
                PostMean <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
        }
        
        if(scores=="percentBL"){
                BLmeanAll <- mean(rowMeans(mat[,BLbinIndex], na.rm=T), na.rm=T) #Mean all subjects PRE infusion
                PostMeanAll <- mean(rowMeans(mat[,PostInfBinIndex], na.rm=T), na.rm=T) #Mean all subjects POST infusion
                
                BLmeanEach <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
                PostMeanEach <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
                
                BLmean <- (BLmeanEach/BLmeanEach)*100 #Mean by subject PRE infusion in terms of percentage of BL performance of that same subject (it has to be 100%)
                PostMean <- (PostMeanEach/BLmeanEach)*100 #Mean by subject POST infusion in terms of percentage of BL performance 
        }
        
        lines(x=c(0, 1), y=c(mean(BLmean), mean(PostMean)), col=color, cex=2)
        errBars(x=c(0, 1), y=c(mean(BLmean), mean(PostMean)), err=c(sd(BLmean)/sqrt(length(BLmean)), sd(PostMean)/sqrt(length(PostMean))), color=color, jitter=jitter)
        points(x=c(0, 1), y=c(mean(BLmean), mean(PostMean)), pch=pch, col=color, cex=2)
        if(pch==22){points(x=c(0, 1), y=c(mean(BLmean), mean(PostMean)), pch=pch, col=color, cex=2, bg="white")}
        
}

colindx <- c("blue", "red")


#### 1.1. RESPONSE RATIO

### 1.1.1. Response ratio: S+ and S- responding pre vs. post infusion in AP5 vs. VEH

## 1.1.1.1. Absolute scores
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 1))

plotPrePostLines(data=byBinDataEarlyVEH[[1]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[1]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

plotPrePostLines(data=byBinDataEarlyVEH[[2]], color=colindx[1], pch=22, scores="absolute", jitter=0.015) #VEH group, S-
plotPrePostLines(data=byBinDataEarlyAP5[[2]], color=colindx[2], pch=22, scores="absolute", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1, font=2)
axis(side=2, at=seq(0, 1, by=0.2, labels=seq(0, 1, 0.2)), font=2, las=2, pos=-0.1)

legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")

## 1.1.1.2. Percentage of BL 
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 120))

plotPrePostLines(data=byBinDataEarlyVEH[[1]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[1]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

#plotPrePostLines(data=byBinDataEarlyVEH[[2]], color=colindx[1], pch=22, scores="percentBL", jitter=0.015) #VEH group, S-. It's confusing so I'm not plotting it
#plotPrePostLines(data=byBinDataEarlyAP5[[2]], color=colindx[2], pch=22, scores="percentBL", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1, font=2)
axis(side=2, at=seq(0, 120, by=20), labels=seq(0, 120, 20), font=2, las=2, pos=-0.1)
mtext(side=2, line=3, text="% of BL response ratio", font=2, cex.axis=1.5)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")


### 1.1.2. Response ratio: S+ and S- responding by bin on test day in AP5 vs. VEH

## 1.1.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 1))

#lapply(seq(1, length(ratsEarlyVEH)), function(x) {lines(byBinDataEarlyVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsEarlyAP5)), function(x) {lines(byBinDataEarlyAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataEarlyVEH[[1]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[1]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)


## 1.1.2.2. Percentage of BL



#### 1.2. CUED LATENCY

### 1.2.1. Cued latency: S+ and S- latency pre vs. post infusion in AP5 vs. VEH

## 1.2.1.1. Absolute scores

## 1.2.1.2. Percentage of BL 

### 1.2.2. Cued latency: S+ and S- responding by bin on test day in AP5 vs. VEH

## 1.2.2.1. Absolute scores

## 1.2.2.2. Percentage of BL




#### 1.3. ITI latency

### 1.3.1. ITI latency: ITI latency pre vs. post infusion in AP5 vs. VEH

## 1.3.1.1. Absolute scores

## 1.3.1.2. Percentage of BL 

### 1.3.2. ITI latency: ITI latency by bin on test day in AP5 vs. VEH

## 1.3.2.1. Absolute scores

## 1.3.2.2. Percentage of BL



#### 1.4. CUED SPECIFICITY

### 1.4.1. Cue specificity: S+ and S- specificity pre vs. post infusion in AP5 vs. VEH

## 1.4.1.1. Absolute scores

## 1.4.1.2. Percentage of BL 

### 1.4.2.  Cue specificity: S+ and S- specificity by bin on test day in AP5 vs. VEH

## 1.4.2.1. Absolute scores

## 1.4.2.2. Percentage of BL














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
