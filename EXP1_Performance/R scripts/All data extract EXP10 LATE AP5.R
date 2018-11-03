

#############################################################
### EXPERIMENT 1A: Late AP5 VS VEH TEST                  ###
#############################################################

library(matrixStats)
library(ez)

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
load(file=paste(funcdirect, "dotplot.r", sep=""))
load(file=paste(funcdirect, "KC.sigbins.R", sep=""))
load(file=paste(funcdirect, "KC.inhib.sigbins.R", sep=""))
load(file=paste(funcdirect, "prePostInf_FR.r", sep=""))

###########################
###########################
### DEFINE FOLDERS      ###
###########################
###########################


# Define folders for one group *OR* the other before running the rest of the code (not both because then you'll just rewrite the folders you defined for the first group)

### LATE VEH #################################################################
subTestFolder <- paste(Exp1folder, "Late VEH/", sep="")
datafolder <- paste(subTestFolder, "MedPC files/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
dataForRCumulativeLateVEH <- dataForRCumulative
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")
preVsPostFRFolder <- "E:/Dropbox/NMDA/EXP1_Performance/Graphs/Neuronal/FR pre vs post scatterplot/"


### LATE AP5 #################################################################
subTestFolder <- paste(Exp1folder, "LATE AP5/", sep="")
datafolder <- paste(subTestFolder, "MedPC files/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
dataForRCumulativeLateAP5 <- dataForRCumulative
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")
preVsPostFRFolder <- "E:/Dropbox/NMDA/EXP1_Performance/Graphs/Neuronal/FR pre vs post scatterplot/"


################################
################################
### CREATE IMPORTANT OBJECTS ###
################################
################################

# Extract behavioral data from MedPC files. This function saves the generated objects in the "dataForRdir". You have to load them (see next line of code) to bring them to your environment.
# This will give you a few error messages if, in any file, the first cue comes on after 5s of session onset. Ignore it, it just assigns NA to that trial, which is what you want.
# The parameter 'consumeRewWdw' is just the segment of the ITI that we discard (for ITI latency calculations) bc we assume that, if the animal got a reward on the previous trial, he might still be consuming the reward.
MedPCextract(cuelength=10, MovAvg = "Impinged only", funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative)

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

# Make an object with the DStaskAccByBin for later statistical analyses
# Change "drug="VEH"" when I run this for AP5
DStaskAccByBin_LongFormat <- do.call("rbind", lapply(seq(1, length(DStaskAccByBin)), function(k){
        rat=as.character(rats[[k]])
        a <- 1:length(DStaskAccByBin[[k]])
        data.frame(rat=rat, bin=a, drug="AP5", perf=DStaskAccByBin[[k]])
}))


### Make a long-format object with all these data for statistical analyses

#Run all of the above lines FIRST for VEH rats and then this line:
byBinDataLateVEH <- list(DSrespRatioByBin, NSrespRatioByBin, DSlatencyByBin, NSlatencyByBin, DStaskAccByBin, NStaskAccByBin, ITIlatByBin)
IndexLabel <- c("S+.RR", "S-.RR", "S+.Latency", "S-.Latency", "S+.Spec.", "S-.Spec.", "ITI.Latency.")

LateVEH_LongFormat <- do.call("rbind", lapply(seq(1, length(byBinDataLateVEH)), function(x){ #For each index
        mat <- do.call("rbind", byBinDataLateVEH[[x]])
        BLmean <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
        PostMean <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
        ratnames <- 1:nrow(mat)
        return(data.frame(Drug="VEH", Rat=ratnames, Index=IndexLabel[x], Infusion=c(rep("Pre", nrow(mat)), rep("Post", nrow(mat))), Performance=c(BLmean, PostMean)))
}))

LateVEH_DStaskAccByBin_LongFormat <- DStaskAccByBin_LongFormat
save(LateVEH_LongFormat, file=paste(dataForRdir, "LateVEH_LongFormat.rdat", sep=""))
save(LateVEH_DStaskAccByBin_LongFormat, file=paste(dataForRdir, "LateVEH_DStaskAccByBin_LongFormat.rdat", sep=""))


#Then repeat for AP5 rats and run these lines
byBinDataLateAP5 <- list(DSrespRatioByBin, NSrespRatioByBin, DSlatencyByBin, NSlatencyByBin, DStaskAccByBin, NStaskAccByBin, ITIlatByBin)
IndexLabel <- c("S+.RR", "S-.RR", "S+.Latency", "S-.Latency", "S+.Spec.", "S-.Spec.", "ITI.Latency.")

LateAP5_LongFormat <- do.call("rbind", lapply(seq(1, length(byBinDataLateAP5)), function(x){ #For each index
        mat <- do.call("rbind", byBinDataLateAP5[[x]])
        BLmean <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
        PostMean <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
        ratnames <- 1:nrow(mat)
        return(data.frame(Drug="AP5", Rat=ratnames, Index=IndexLabel[x], Infusion=c(rep("Pre", nrow(mat)), rep("Post", nrow(mat))), Performance=c(BLmean, PostMean)))
}))

LateAP5_DStaskAccByBin_LongFormat <- DStaskAccByBin_LongFormat

save(LateAP5_LongFormat, file=paste(dataForRdir, "LateAP5_LongFormat.rdat", sep=""))
save(LateAP5_DStaskAccByBin_LongFormat, file=paste(dataForRdir, "LateAP5_DStaskAccByBin_LongFormat.rdat", sep=""))


###
Late_LongFormat <- rbind(LateVEH_LongFormat, LateAP5_LongFormat)
Late_LongFormatByBin <- rbind(LateVEH_DStaskAccByBin_LongFormat, LateAP5_DStaskAccByBin_LongFormat)

save(Late_LongFormat, file=paste(dataForRdir, "Late_LongFormat.rdat", sep=""))
save(Late_LongFormatByBin, file=paste(dataForRdir, "Late_LongFormatByBin.rdat", sep=""))



# Extract neuronal data from NEX files. 

#VEH test: data aligned to DS and NS onset BEFORE and AFTER the infusion.
postInfTargetWdw <- 1800+12*60+30*60 #For the post infusion window, I'll choose the period between the end of the infusion +30'.
allNeuronsDSLateVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsDSLateVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=1800, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSLateVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSLateVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryDSLateVEHPreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryDSLateVEHPostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryNSLateVEHPreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryNSLateVEHPostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryITILateVEHPreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryITILateVEHPostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")





#AP5 test: data aligned to DS and NS onset BEFORE and AFTER the infusion.
postInfTargetWdw <- 1800+12*60+30*60 #For the post infusion window, I'll choose the period between the end of the infusion +30'.
allNeuronsDSLateAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsDSLateAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSLateAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSLateAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryDSLateAP5PreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryDSLateAP5PostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryNSLateAP5PreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryNSLateAP5PostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryITILateAP5PreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryITILateAP5PostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")



### GIVE THESE OBJECTS A UNIQUE NAME

# ## VEH SIDE
# csacqidxLateVEH <- csacqidx
# alldataLateVEH <- alldata
# ratsLateVEH <- rats
# idxLateVEH <- idx
# cumDataLateVEH <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)
# 
# ## AP5 SIDE
# csacqidxLateAP5 <- csacqidx
# alldataLateAP5 <- alldata
# ratsLateAP5 <- rats
# idxLateAP5 <- idx
# cumDataLateAP5 <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)



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


#Function for plotting bars more easily. I just need to adjust the data I feed the function, the color and the points
plotPrePostBars <- function(data, color, xmiddle, barwidth, labelY, colLabel){
        mat <- do.call("rbind", data) #Create matrix in which rows are different rats and columns are bins
        BLmean <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
        PostMean <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
        
        #Pre
        rect(xleft=xmiddle-barwidth, xright = xmiddle, ybottom=0, ytop=mean(BLmean), col=color, border="white")
        text(x=xmiddle-barwidth/2, y=labelY, labels = "Pre", col=colLabel, font=2)
        #Post
        rect(xleft=xmiddle, xright=xmiddle+barwidth, ybottom=0, ytop=mean(PostMean), col=color, border="white")
        text(x=xmiddle+barwidth/2, y=labelY, labels = "Post", col=colLabel, font=2)
        #Individual lines
        for(i in 1:length(data)){lines(x=c(xmiddle-barwidth/2, xmiddle+barwidth/2), y=c(BLmean[i], PostMean[i]))}
}

colindx <- c("#2171b5", "#cb181d") #Strong blue and red
colindxB <- c("#bdd7e7", "#fcae91") #Less strong blue and red
colindxC <- c("#eff3ff", "#fb6a4a") #Even less strong blue and red
colindxD <- c("#6baed6", "#fee5d9") #Lightest blue and red



#### 1.1. RESPONSE RATIO

### 1.1.1. Response ratio: S+ and S- responding pre vs. post infusion in AP5 vs. VEH
# In the objects 'byBinDataLateVEH' and 'byBinDataLateAP5', the first and second items are DSrespratio and NSrespratio by subject by bin

## 1.1.1.1. Absolute scores
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 1))

plotPrePostLines(data=byBinDataLateVEH[[1]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataLateAP5[[1]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

plotPrePostLines(data=byBinDataLateVEH[[2]], color=colindx[1], pch=22, scores="absolute", jitter=0.015) #VEH group, S-
plotPrePostLines(data=byBinDataLateAP5[[2]], color=colindx[2], pch=22, scores="absolute", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 1, by=0.2, labels=seq(0, 1, 0.2)), font=2, las=2, pos=-0.1)
mtext(side=2, line=4, text="Proportion", cex=1.4, font=2)

legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")

## 1.1.1.2. Percentage of BL 
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 120))

plotPrePostLines(data=byBinDataLateVEH[[1]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataLateAP5[[1]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

#plotPrePostLines(data=byBinDataLateVEH[[2]], color=colindx[1], pch=22, scores="percentBL", jitter=0.015) #VEH group, S-. It's confusing so I'm not plotting it
#plotPrePostLines(data=byBinDataLateAP5[[2]], color=colindx[2], pch=22, scores="percentBL", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 120, by=20), labels=seq(0, 120, 20), font=2, las=2, pos=-0.1)
mtext(side=2, line=4, text="% of BL response ratio", font=2, cex=1.4)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")


### 1.1.2. Response ratio: S+ and S- responding by bin on test day in AP5 vs. VEH

## 1.1.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 1))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=1.2, col="gray95", border="white")

#lapply(seq(1, length(ratsLateVEH)), function(x) {lines(byBinDataLateVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsLateAP5)), function(x) {lines(byBinDataLateAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataLateVEH[[1]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[1]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 1, by=0.2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.4)
mtext(side=2, line=1, text="Proportion", font=2, cex=1.4)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)



## 1.1.2.2. Percentage of BL
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 120))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=120, col="gray95", border="white")

#Get data ready
matVEH <- do.call("rbind", byBinDataLateVEH[[1]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[1]])

VEHbl <- colMeans(matVEH[,BLbinIndex], na.rm=T); AP5bl <- colMeans(matAP5[,BLbinIndex], na.rm=T)

matVEHperc <- (matVEH/VEHbl)*100
matAP5perc <- (matAP5/AP5bl)*100

#Plot
lines(colMeans(matVEHperc), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEHperc), err=colSds(matVEHperc)/sqrt(nrow(matVEHperc)), color=colindx[1])
points(colMeans(matVEHperc), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5perc), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5perc), err=colSds(matAP5perc)/sqrt(nrow(matAP5perc)), color=colindx[2])
points(colMeans(matAP5perc), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 120, by=20), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.4)
mtext(side=2, line=1, text="% of baseline", font=2, cex=1.4)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)


### 1.1.3. Response ratio: barplots of pre and post infusion, S+ vs S- and VEH vs AP5
plot.new()

plot.window(xlim=c(0, 6), ylim=c(0, 1))

#S+ both groups pre and post infusion
plotPrePostBars(data=byBinDataLateVEH[[1]], color=colindx[1], xmiddle=1, barwidth=0.5, colLabel = "white", labelY = 0.05)
plotPrePostBars(data=byBinDataLateAP5[[1]], color=colindx[2], xmiddle=2, barwidth=0.5, colLabel = "white", labelY = 0.05)

#S- both groups pre and post infusion
plotPrePostBars(data=byBinDataLateVEH[[2]], color=colindxB[1], xmiddle=3.25, barwidth=0.5, colLabel = "black", labelY = 0.05)
plotPrePostBars(data=byBinDataLateAP5[[2]], color=colindxB[2], xmiddle=4.25, barwidth=0.5, colLabel = "black", labelY = 0.05)

#Axis and labels
axis(side=1, tick = F, at=c(1.5, 3.75), labels=c("S+", "S-"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 1, 0.2), cex.axis=1, font=2, las=2)
mtext(side=2, line=2.5, text="Response ratio", cex=1.4, font=2)
rect(xleft=3, xright=3.5, ybottom=1.15, ytop=1.2, col=colindx[1], border="white")
rect(xleft=3.5, xright=4, ybottom=1.15, ytop=1.2, col=colindxB[1], border="white")
rect(xleft=3, xright=3.5, ybottom=1.05, ytop=1.1, col=colindx[2], border="white")
rect(xleft=3.5, xright=4, ybottom=1.05, ytop=1.1, col=colindxB[2], border="white")
text(x=4.5, y=1.18, labels="VEH", cex=1.5)
text(x=4.5, y=1.08, labels="AP5", cex=1.5)



#### 1.2. CUED LATENCY
# In the objects 'byBinDataLateVEH' and 'byBinDataLateAP5', the 3rd and 4th items are DSlatency and NSlatency by subject by bin

### 1.2.1. Cued latency: S+ and S- latency pre vs. post infusion in AP5 vs. VEH

## 1.2.1.1. Absolute scores
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 10))

plotPrePostLines(data=byBinDataLateVEH[[3]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataLateAP5[[3]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

plotPrePostLines(data=byBinDataLateVEH[[4]], color=colindx[1], pch=22, scores="absolute", jitter=0.015) #VEH group, S-
plotPrePostLines(data=byBinDataLateAP5[[4]], color=colindx[2], pch=22, scores="absolute", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 10, by=2, labels=seq(0, 10, 2)), font=2, las=2, pos=-0.05)

mtext(side=2, line=2.5, text="Latency (s)", cex=1.4, font=2)

legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")


## 1.2.1.2. Percentage of BL 
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 400))

plotPrePostLines(data=byBinDataLateVEH[[3]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataLateAP5[[3]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

#plotPrePostLines(data=byBinDataLateVEH[[2]], color=colindx[1], pch=22, scores="percentBL", jitter=0.015) #VEH group, S-. It's confusing so I'm not plotting it
#plotPrePostLines(data=byBinDataLateAP5[[2]], color=colindx[2], pch=22, scores="percentBL", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 400, by=50), labels=seq(0, 400, 50), font=2, las=2, pos=-0.05)
mtext(side=2, line=3, text="% of BL latency", font=2, cex=1.4)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n", cex=1.5)

### 1.2.2. Cued latency: S+ and S- responding by bin on test day in AP5 vs. VEH

## 1.2.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 10))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=10, col="gray95", border="white")

#lapply(seq(1, length(ratsLateVEH)), function(x) {lines(byBinDataLateVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsLateAP5)), function(x) {lines(byBinDataLateAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataLateVEH[[3]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[3]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 10, by=2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.4)
mtext(side=2, line=1, text="Latency (s)", font=2, cex=1.4)

legend("topright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)



## 1.2.2.2. Percentage of BL
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 300))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=300, col="gray95", border="white")

#Get data ready
matVEH <- do.call("rbind", byBinDataLateVEH[[3]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[3]])

VEHbl <- colMeans(matVEH[,BLbinIndex], na.rm=T); AP5bl <- colMeans(matAP5[,BLbinIndex], na.rm=T)

matVEHperc <- (matVEH/VEHbl)*100
matAP5perc <- (matAP5/AP5bl)*100

#Plot
lines(colMeans(matVEHperc), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEHperc), err=colSds(matVEHperc)/sqrt(nrow(matVEHperc)), color=colindx[1])
points(colMeans(matVEHperc), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5perc), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5perc), err=colSds(matAP5perc)/sqrt(nrow(matAP5perc)), color=colindx[2])
points(colMeans(matAP5perc), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 300, by=50), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="% of baseline latency", font=2, cex=1.2)

legend("topright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)

### 1.2.3. Cued latency: barplots of pre and post infusion, S+ vs S- and VEH vs AP5
plot.new()
plot.window(xlim=c(0, 6), ylim=c(0, 10))

#S+ both groups pre and post infusion
plotPrePostBars(data=byBinDataLateVEH[[3]], color=colindx[1], xmiddle=1, barwidth=0.5, colLabel = "white", labelY = 0.5)
plotPrePostBars(data=byBinDataLateAP5[[3]], color=colindx[2], xmiddle=2, barwidth=0.5, colLabel = "white", labelY = 0.5)

#S- both groups pre and post infusion
plotPrePostBars(data=byBinDataLateVEH[[4]], color=colindxB[1], xmiddle=3.25, barwidth=0.5, colLabel = "black", labelY = 0.5)
plotPrePostBars(data=byBinDataLateAP5[[4]], color=colindxB[2], xmiddle=4.25, barwidth=0.5, colLabel = "black", labelY = 0.5)

#Axis and labels
axis(side=1, tick = F, at=c(1.5, 3.75), labels=c("S+", "S-"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 10, 2), cex.axis=1.4, font=2, las=2)
mtext(side=2, line=2.5, text="Latency (s)", cex=1.4, font=2)
rect(xleft=0, xright=0.5, ybottom=9.5, ytop=10, col=colindx[1], border="white")
rect(xleft=0.5, xright=1, ybottom=9.5, ytop=10, col=colindxB[1], border="white")
rect(xleft=0, xright=0.5, ybottom=8.5, ytop=9, col=colindx[2], border="white")
rect(xleft=0.5, xright=1, ybottom=8.5, ytop=9, col=colindxB[2], border="white")
text(x=1.5, y=9.8, labels="VEH", cex=1.5)
text(x=1.5, y=8.8, labels="AP5", cex=1.5)



#### 1.3. ITI latency

### 1.3.1. ITI latency: ITI latency pre vs. post infusion in AP5 vs. VEH

## 1.3.1.1. Absolute scores
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 10))

plotPrePostLines(data=byBinDataLateVEH[[7]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataLateAP5[[7]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 10, by=2), font=2, las=2, pos=-0.05)

mtext(side=2, line=2.5, text="ITI latency (s)", cex=1.4, font=2)

legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")

## 1.3.1.2. Percentage of BL 
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 150))

plotPrePostLines(data=byBinDataLateVEH[[7]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataLateAP5[[7]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 200, by=50), font=2, las=2, pos=-0.05)
mtext(side=2, line=3, text="% of BL ITI latency", font=2, cex=1.4)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomright", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n", cex=1.5)


### 1.3.2. ITI latency: ITI latency by bin on test day in AP5 vs. VEH

## 1.3.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 10))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=10, col="gray95", border="white")

#lapply(seq(1, length(ratsLateVEH)), function(x) {lines(byBinDataLateVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsLateAP5)), function(x) {lines(byBinDataLateAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataLateVEH[[7]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[7]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 10, by=2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.4)
mtext(side=2, line=1, text="ITI latency (s)", font=2, cex=1.4)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)


## 1.3.2.2. Percentage of BL
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 150))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=150, col="gray95", border="white")

#Get data ready
matVEH <- do.call("rbind", byBinDataLateVEH[[7]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[7]])

VEHbl <- colMeans(matVEH[,BLbinIndex], na.rm=T); AP5bl <- colMeans(matAP5[,BLbinIndex], na.rm=T)

matVEHperc <- (matVEH/VEHbl)*100
matAP5perc <- (matAP5/AP5bl)*100

#Plot
lines(colMeans(matVEHperc), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEHperc), err=colSds(matVEHperc)/sqrt(nrow(matVEHperc)), color=colindx[1])
points(colMeans(matVEHperc), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5perc), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5perc), err=colSds(matAP5perc)/sqrt(nrow(matAP5perc)), color=colindx[2])
points(colMeans(matAP5perc), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 150, by=50), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="% of baseline ITI latency", font=2, cex=1.2)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)


### 1.3.3. ITI latency: barplots of pre and post infusion, S+ vs S- and VEH vs AP5
plot.new()
plot.window(xlim=c(0, 6), ylim=c(0, 10))

#S+ both groups pre and post infusion
plotPrePostBars(data=byBinDataLateVEH[[7]], color=colindx[1], xmiddle=1, barwidth=0.5, colLabel = "white", labelY = 0.5)
plotPrePostBars(data=byBinDataLateAP5[[7]], color=colindx[2], xmiddle=2, barwidth=0.5, colLabel = "white", labelY = 0.5)

#Axis and labels
axis(side=2, at=seq(0, 10, 2), cex.axis=1, font=2, las=2)
mtext(side=2, line=2.5, text="ITI latency (s)", cex=1.4, font=2)
rect(xleft=0, xright=0.5, ybottom=9.5, ytop=10, col=colindx[1], border="white")
rect(xleft=0, xright=0.5, ybottom=8.5, ytop=9, col=colindx[2], border="white")
text(x=1, y=9.8, labels="VEH", cex=1.5)
text(x=1, y=8.8, labels="AP5", cex=1.5)







#### 1.4. CUED SPECIFICITY

### 1.4.1. Cue specificity: S+ and S- specificity pre vs. post infusion in AP5 vs. VEH
# In the objects 'byBinDataLateVEH' and 'byBinDataLateAP5', the 5th and 6th items are DStaskAccuracy and NStaskAccuracy by subject by bin

## 1.4.1.1. Absolute scores
plot.new()
plot.window(xlim=c(0, 1), ylim=c(-2, 6))

abline(h=0, lty=3)

plotPrePostLines(data=byBinDataLateVEH[[5]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataLateAP5[[5]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

plotPrePostLines(data=byBinDataLateVEH[[6]], color=colindx[1], pch=22, scores="absolute", jitter=0.015) #VEH group, S-
plotPrePostLines(data=byBinDataLateAP5[[6]], color=colindx[2], pch=22, scores="absolute", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(-2, 6, by=1), font=2, las=2, pos=-0.04)
mtext(side=2, line=2.5, text="Specificity score", font=2, cex=1.4)

legend("topright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("topleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")

## 1.4.1.2. Percentage of BL 
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 140))

plotPrePostLines(data=byBinDataLateVEH[[5]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataLateAP5[[5]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

#plotPrePostLines(data=byBinDataLateVEH[[2]], color=colindx[1], pch=22, scores="percentBL", jitter=0.015) #VEH group, S-. It's confusing so I'm not plotting it
#plotPrePostLines(data=byBinDataLateAP5[[2]], color=colindx[2], pch=22, scores="percentBL", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 140, by=20), font=2, las=2, pos=-0.05)
mtext(side=2, line=3.2, text="% of BL specificity", font=2, cex=1.4)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n", cex=1.2)


### 1.4.2.  Cue specificity: S+ and S- specificity by bin on test day in AP5 vs. VEH

## 1.4.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(-2, 7))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=-2, ytop=6, col="gray95", border="white")
abline(h=0, lty=3)

#lapply(seq(1, length(ratsLateVEH)), function(x) {lines(byBinDataLateVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsLateAP5)), function(x) {lines(byBinDataLateAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataLateVEH[[5]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[5]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(-2, 6, by=2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.4)
mtext(side=2, line=1, text="S+ specificity (s)", font=2, cex=1.4)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)




## 1.4.2.1.B. s- Cued specificity. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(-2, 7))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=-2, ytop=6, col="gray95", border="white")
abline(h=0, lty=3)

matVEH <- do.call("rbind", byBinDataLateVEH[[6]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[6]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(-2, 6, by=2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="S- specificity (s)", font=2, cex=1.2)

legend("topright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)



## 1.4.2.2. Percentage of BL
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 150))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=140, col="gray95", border="white")
abline(h=100, lty=3)

#Get data ready
matVEH <- do.call("rbind", byBinDataLateVEH[[5]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataLateAP5[[5]])

VEHbl <- colMeans(matVEH[,BLbinIndex], na.rm=T); AP5bl <- colMeans(matAP5[,BLbinIndex], na.rm=T)

matVEHperc <- (matVEH/VEHbl)*100
matAP5perc <- (matAP5/AP5bl)*100

#Plot
lines(colMeans(matVEHperc), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEHperc), err=colSds(matVEHperc)/sqrt(nrow(matVEHperc)), color=colindx[1])
points(colMeans(matVEHperc), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5perc), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5perc), err=colSds(matAP5perc)/sqrt(nrow(matAP5perc)), color=colindx[2])
points(colMeans(matAP5perc), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 140, by=20), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="% of baseline S+ specificity", font=2, cex=1.2)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.2)


### 1.4.3. Cued specificity: barplots of pre and post infusion, S+ vs S- and VEH vs AP5
plot.new()
plot.window(xlim=c(0, 6), ylim=c(-2, 7))

#S+ both groups pre and post infusion
plotPrePostBars(data=byBinDataLateVEH[[5]], color=colindx[1], xmiddle=1, barwidth=0.5, colLabel = "white", labelY = 0.5)
plotPrePostBars(data=byBinDataLateAP5[[5]], color=colindx[2], xmiddle=2, barwidth=0.5, colLabel = "white", labelY = 0.5)

#S- both groups pre and post infusion
plotPrePostBars(data=byBinDataLateVEH[[6]], color=colindxB[1], xmiddle=3.25, barwidth=0.5, colLabel = "black", labelY = 0.5)
plotPrePostBars(data=byBinDataLateAP5[[6]], color=colindxB[2], xmiddle=4.25, barwidth=0.5, colLabel = "black", labelY = 0.5)

#Axis and labels
axis(side=1, tick = F, at=c(1.5, 3.75), labels=c("S+", "S-"), cex.axis=1.4, font=2)
axis(side=2, at=seq(-2, 6, 2), cex.axis=1, font=2, las=2)
mtext(side=2, line=2.5, text="Specificity score", cex=1.4, font=2)
rect(xleft=3, xright=3.5, ybottom=5.5, ytop=6, col=colindx[1], border="white")
rect(xleft=3.5, xright=4, ybottom=5.5, ytop=6, col=colindxB[1], border="white")
rect(xleft=3, xright=3.5, ybottom=4.5, ytop=5, col=colindx[2], border="white")
rect(xleft=3.5, xright=4, ybottom=4.5, ytop=5, col=colindxB[2], border="white")
text(x=4.5, y=5.8, labels="VEH", cex=1.5)
text(x=4.5, y=4.8, labels="AP5", cex=1.5)







######################################################
######################################################
### STATISTICAL ANALYSES                           ###
######################################################
######################################################
# I'll conduct two-way (VEH vs. AP5 and Pre vs. Post infusion as within-subjects factors) repeated measures ANOVA for each index. I'll take the 30min before infusion vs. 30min after infusion.

Late_LongFormat #This is our object of reference

indexes <- unique(Late_LongFormat$Index)


## S+ Response ratio
DSRR <- subset(x=Late_LongFormat, Late_LongFormat$Index==indexes[1])
ezANOVA(data=DSRR, dv=Performance, within=.(Infusion, Drug), wid=Rat, type=3)
#Nothing was significant
# $ANOVA
#          Effect DFn DFd         F          p p<.05       ges
# 2      Infusion   1   4 5.9428499 0.07137787       0.2994483
# 3          Drug   1   4 0.4490463 0.53948830       0.0454966
# 4 Infusion:Drug   1   4 0.4490463 0.53948830       0.0454966


## S- Response ratio
NSRR <- subset(x=Late_LongFormat, Late_LongFormat$Index==indexes[2])
ezANOVA(data=NSRR, dv=Performance, within=.(Infusion, Drug), wid=Rat, type=3) 

#Nothing was significant
# $ANOVA
#          Effect DFn DFd          F         p p<.05          ges
# 2      Infusion   1   4 0.85982604 0.4062723       0.0164609238
# 3          Drug   1   4 0.06168540 0.8160816       0.0061376527
# 4 Infusion:Drug   1   4 0.03332762 0.8640235       0.0003347305

## S+ latency
DSlat <- subset(x=Late_LongFormat, Late_LongFormat$Index==indexes[3])
ezANOVA(data=DSlat, dv=Performance, within=.(Infusion, Drug), wid=Rat, type=3)

# Nothing was significant
# $ANOVA
#          Effect DFn DFd         F         p p<.05        ges
# 2      Infusion   1   4 4.1531576 0.1112051       0.13331434
# 3          Drug   1   4 0.4475211 0.5401433       0.05564145
# 4 Infusion:Drug   1   4 1.3249969 0.3138235       0.05007153


## S- latency
NSlat <- subset(x=Late_LongFormat, Late_LongFormat$Index==indexes[4])
ezANOVA(data=NSlat, dv=Performance, within=.(Infusion, Drug), wid=Rat, type=3)

#Nothing was significant:
# $ANOVA
#          Effect DFn DFd           F         p p<.05          ges
# 2      Infusion   1   4 1.743198156 0.2572151       0.0261640644
# 3          Drug   1   4 0.001622495 0.9698000       0.0002049124
# 4 Infusion:Drug   1   4 0.116647144 0.7498877       0.0010269948


## ITI latency
ITIlat <- subset(x=Late_LongFormat, Late_LongFormat$Index==indexes[7])
ezANOVA(data=ITIlat, dv=Performance, within=.(Infusion, Drug), wid=Rat, type=3)
# $ANOVA
#          Effect DFn DFd           F         p p<.05          ges
# 2      Infusion   1   4 1.323862618 0.3140043       0.0632833364
# 3          Drug   1   4 0.001660553 0.9694482       0.0001713418
# 4 Infusion:Drug   1   4 2.783263588 0.1705790       0.0236162481


## S+ specificity
DSspec <- subset(x=Late_LongFormat, Late_LongFormat$Index==indexes[5])
ezANOVA(data=DSspec, dv=Performance, within=.(Infusion, Drug), wid=Rat, type=3)

# Nothing was significant
# $ANOVA
#          Effect DFn DFd          F         p p<.05         ges
# 2      Infusion   1   4 0.06054481 0.8177474       0.001005104
# 3          Drug   1   4 0.09550315 0.7727219       0.015299556
# 4 Infusion:Drug   1   4 0.24627803 0.6457405       0.006770484

## S- specificity
NSspec <- subset(x=Late_LongFormat, Late_LongFormat$Index==indexes[6])
ezANOVA(data=NSspec, dv=Performance, within=.(Infusion, Drug), wid=Rat, type=3)

# $ANOVA
#           Effect DFn DFd          F          p p<.05         ges
# 2      Infusion   1   4  0.3150865 0.60453326       0.016660731
# 3          Drug   1   4  0.0147881 0.90907516       0.001544113
# 4 Infusion:Drug   1   4 11.0585371 0.02922977     * 0.082284052

#The interaction was significant. As post-hoc test, I'll split the dataset into the groups and, within each group, use a paired t-test for the pre vs. post
NSspec_VEH <- subset(NSspec, Drug=="VEH")
NSspec_AP5 <- subset(NSspec, Drug=="AP5")

vehtest <- t.test(x=NSspec_VEH$Performance[NSspec_VEH$Infusion=="Pre"], y=NSspec_VEH$Performance[NSspec_VEH$Infusion=="Post"], paired=T, alternative="two.sided") #t(3)= 0.50965, p=0.6454191
ap5test <- t.test(x=NSspec_AP5$Performance[NSspec_AP5$Infusion=="Pre"], y=NSspec_AP5$Performance[NSspec_AP5$Infusion=="Post"], paired=T, alternative="two.sided") #t(3)=-1.2502, p=0.5997133
p.adjust(p=c(vehtest$p.value, ap5test$p.value), method="holm")


##############################

### PERFORMANCE INDEX BY BIN
Late_LongFormatByBin$bin <- as.factor(Late_LongFormatByBin$bin)

bins <- unique(Late_LongFormatByBin$bin)
newBinsIndex <- c(1, 1, 1, 0, 2, 2, 2, 3, 3, 3, 4, 4) #I want to create 30min bins instead of 10min bins because to compare so many bins reduces my p values a lot when adjusting

newBinsVals <- sapply(seq(1, nrow(Late_LongFormatByBin)), function(l){
        sel <- Late_LongFormatByBin$bin[l]
        newBinsIndex[sel]
})

Late_LongFormatByBin$Bigbins <- newBinsVals

smallbins.aov <- summary(aov(perf ~ drug * bins + Error(rat/(bins*drug)), data=Late_LongFormatByBin))
bigbins.aov <- summary(aov(perf ~ drug * Bigbins + Error(rat/(Bigbins*drug)), data=Late_LongFormatByBin))

ezANOVA(data=Late_LongFormatByBin, dv=perf, within=.(bin, drug), wid=rat, type=1)
# $`ANOVA`
#     Effect DFn DFd         F         p p<.05        ges
# 1      bin  11  44 1.3622909 0.2246165       0.13232317
# 2     drug   1   4 2.2513012 0.2078935       0.01801055
# 3 bin:drug  11  44 0.6863665 0.7438122       0.08186467


newBinsVals <- sapply(seq(1, nrow(LateVEH_DStaskAccByBin_LongFormat)), function(l){
        sel <- LateVEH_DStaskAccByBin_LongFormat$bin[l]
        newBinsIndex[sel]
})

LateVEH_DStaskAccByBin_LongFormat$Bigbins <- newBinsVals

newBinsVals <- sapply(seq(1, nrow(LateAP5_DStaskAccByBin_LongFormat)), function(l){
        sel <- LateAP5_DStaskAccByBin_LongFormat$bin[l]
        newBinsIndex[sel]
})

LateAP5_DStaskAccByBin_LongFormat$Bigbins <- newBinsVals


ttestPerBigBin <- do.call("rbind", lapply(seq(1, length(unique(Late_LongFormatByBin$Bigbins))), function(m){
        bindex <- unique(Late_LongFormatByBin$Bigbins)[m]
        tst <- t.test(x=LateVEH_DStaskAccByBin_LongFormat$perf[LateVEH_DStaskAccByBin_LongFormat$Bigbins==bindex], 
                      y=LateAP5_DStaskAccByBin_LongFormat$perf[LateAP5_DStaskAccByBin_LongFormat$Bigbins==bindex], paired=T, alternative="greater") 
        
        data.frame(Bigbins=bindex, t=tst$statistic, df=tst$parameter, p=tst$p.value)
})
)

#Adjust the t test p values and also the p values of the ANOVA (using the big bins)
padjusted <- p.adjust(p=c(0.2246165, 0.2078935, 0.7438122, ttestPerBigBin$p), method="holm")


ttestPerBigBin$p.adjusted <- p.adjust(p=ttestPerBigBin$p, method="holm")
#    Bigbins           t df         p p.adjusted
# t        1  0.18791913 14 0.4268183          1
# t1       0  0.86050153  4 0.2190169          1
# t2       2  0.68442297 14 0.2524406          1
# t3       3  0.26590056 14 0.3970958          1
# t4       4 -0.05136734  9 0.5199226          1


#################################################################
####### NEURONAL FIRING                                ##########
#################################################################

###############################################################
#1. PSTH pre and post infusion
###############################################################

#########
## VEH ##
#########

# S+ Onset
psthInf(formatDat="Zscores", group="VEH", event="S+", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), stimulus="cue", imgFormat="pdf")       
psthInf(formatDat="raw", group="VEH", event="S+", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=26, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), stimulus="cue", imgFormat="pdf")       

#S- Onset
psthInf(formatDat="Zscores", group="VEH", event="S-", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsNSLateVEHPreInf, allNeuronsNSLateVEHPostInf), stimulus="cue", imgFormat="pdf")       
psthInf(formatDat="raw", group="VEH", event="S-", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=26, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsNSLateVEHPreInf, allNeuronsNSLateVEHPostInf), stimulus="cue", imgFormat="pdf")       


#S+ Entry
psthInf(formatDat="Zscores", group="VEH", event="S+ Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryDSLateVEHPreInf, allNeuronsEntryDSLateVEHPostInf), stimulus="entry",  BLNeuData=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="VEH", event="S+ Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryDSLateVEHPreInf, allNeuronsEntryDSLateVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), imgFormat="pdf")       

#S- Entry
psthInf(formatDat="Zscores", group="VEH", event="S- Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryNSLateVEHPreInf, allNeuronsEntryNSLateVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="VEH", event="S- Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryNSLateVEHPreInf, allNeuronsEntryNSLateVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), imgFormat="pdf")       

#ITI Entry
psthInf(formatDat="Zscores", group="VEH", event="ITI Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryITILateVEHPreInf, allNeuronsEntryITILateVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="VEH", event="ITI Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryITILateVEHPreInf, allNeuronsEntryITILateVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), imgFormat="pdf")       


########
#AP5
########

#S+ Onset
psthInf(formatDat="Zscores", group="AP5", event="S+", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), stimulus="cue", imgFormat="pdf")       
psthInf(formatDat="raw", group="AP5", event="S+", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=26, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), stimulus="cue", imgFormat="pdf")       

#S- Onset
psthInf(formatDat="Zscores", group="AP5", event="S-", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsNSLateAP5PreInf, allNeuronsNSLateAP5PostInf), stimulus="cue", imgFormat="pdf")       
psthInf(formatDat="raw", group="AP5", event="S-", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=26, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsNSLateAP5PreInf, allNeuronsNSLateAP5PostInf), stimulus="cue", imgFormat="pdf")       


#S+ Entry
psthInf(formatDat="Zscores", group="AP5", event="S+ Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryDSLateAP5PreInf, allNeuronsEntryDSLateAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="AP5", event="S+ Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryDSLateAP5PreInf, allNeuronsEntryDSLateAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), imgFormat="pdf")       


#S- Entry
psthInf(formatDat="Zscores", group="AP5", event="S- Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryNSLateAP5PreInf, allNeuronsEntryNSLateAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="AP5", event="S- Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryNSLateAP5PreInf, allNeuronsEntryNSLateAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), imgFormat="pdf")       


#ITI Entry
psthInf(formatDat="Zscores", group="AP5", event="ITI Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryITILateAP5PreInf, allNeuronsEntryITILateAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="AP5", event="ITI Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Late", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryITILateAP5PreInf, allNeuronsEntryITILateAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), imgFormat="pdf")       



###############################################################
#2. POINTS pre and post infusion around time of cue
###############################################################

dotplot(neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf, allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf),
        expName="Late", dot="Medians", Lines=T, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400, graphFolder = neuGraphFolder)
dotplot(neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf, allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf),
        expName="Late", dot="Means", Lines=T, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400)
dotplot(neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf, allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf),
        expName="Late", dot="Medians", Lines=F, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400)
dotplot(neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf, allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf),
        expName="Late", dot="Means", Lines=F, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400)


#Same but with boxplot instead of dotplot
#Z scores
dotplot(boxplot=T, neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf, allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf),
        expName="Late 100-400", Lines=T, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400, ytop=12, comp=c("VEH", "AP5"))
dotplot(boxplot=T, neudata=list(allNeuronsNSLateVEHPreInf, allNeuronsNSLateVEHPostInf, allNeuronsNSLateAP5PreInf, allNeuronsNSLateAP5PostInf),
        expName="Late 100-400", Lines=T, ytop=12, ybottom=-2, col=colindx, plotWidth=0.3, event="S-", winmin=100, winmax=400, comp=c("VEH", "AP5"))
dotplot(boxplot=T, neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf, allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf),
        expName="Late 100-400", Lines=F, col=colindx, ybottom=-2, plotWidth=0.3, event="S+", winmin=100, winmax=400, ytop=9, comp=c("VEH", "AP5"))
dotplot(boxplot=T, neudata=list(allNeuronsNSLateVEHPreInf, allNeuronsNSLateVEHPostInf, allNeuronsNSLateAP5PreInf, allNeuronsNSLateAP5PostInf),
        expName="Late 100-400", Lines=F, ytop=9, ybottom=-2, col=colindx, plotWidth=0.3, event="S-", winmin=100, winmax=400, comp=c("VEH", "AP5"))

#Raw
dotplot(boxplot=T, neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf, allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf),
        expName="Late 100-400", Lines=T, formatDat="Raw", ytop=20, ybottom=0, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400, comp=c("VEH", "AP5"))
dotplot(boxplot=T, neudata=list(allNeuronsNSLateVEHPreInf, allNeuronsNSLateVEHPostInf, allNeuronsNSLateAP5PreInf, allNeuronsNSLateAP5PostInf),
        expName="Late 100-400", Lines=T, formatDat="Raw", ytop=20, ybottom=0, col=colindx, plotWidth=0.3, event="S-", winmin=100, winmax=400, comp=c("VEH", "AP5"))
dotplot(boxplot=T, neudata=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf, allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf),
        expName="Late 100-400", Lines=F, formatDat="Raw", ytop=20, ybottom=0, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400, comp=c("VEH", "AP5"))
dotplot(boxplot=T, neudata=list(allNeuronsNSLateVEHPreInf, allNeuronsNSLateVEHPostInf, allNeuronsNSLateAP5PreInf, allNeuronsNSLateAP5PostInf),
        expName="Late 100-400", Lines=F, formatDat="Raw", ytop=20, ybottom=0, col=colindx, plotWidth=0.3, event="S-", winmin=100, winmax=400, comp=c("VEH", "AP5"))


###############################################################
#3. SCATTERPLOT pre and post infusion around time of cue
###############################################################

prePostInf_FR(data=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), 
              dataformat="Raw", BLwdw=2, winmin=100, winmax=400, col_labels="purple",
              comparison="Late VEH Pre vs. Post", graphfolder=preVsPostFRFolder)

prePostInf_FR(data=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), 
              dataformat="Raw", BLwdw=2, winmin=100, winmax=400, col_labels="purple",
              comparison="Late AP5 Pre vs. Post", graphfolder=preVsPostFRFolder)


#Baseline
prePostInf_FR(data=list(allNeuronsDSLateVEHPreInf, allNeuronsDSLateVEHPostInf), 
              dataformat="Raw", BLwdw=2, winmin=-2000, winmax=0, col_labels="gray65",
              comparison="Baseline_Late VEH Pre vs. Post", graphfolder=preVsPostFRFolder)

prePostInf_FR(data=list(allNeuronsDSLateAP5PreInf, allNeuronsDSLateAP5PostInf), 
              dataformat="Raw", BLwdw=2, winmin=-2000, winmax=0, col_labels="gray65",
              comparison="Baseline_Late AP5 Pre vs. Post", graphfolder=preVsPostFRFolder)



### STATISTICAL TEST
#Comparison of cue-evoked firing rate (100ms-300ms window) pre vs post injection for both groups. I'll use a Wilcoxon paired test separately for each group
#To get the dotplotDataVEH or dotplotDataAP5 objects, I need to run the code inside the dotplot function separately. I should fix this to make it easier.
#Z scores with 2s BL

load(file=paste(dataForRdir, "dotPlotByGroup.rdat", sep=""))
dotPlotByGroupLate <- dotPlotByGroup
dotplotDataVEH_Late <- dotPlotByGroupLate$VEH
dotplotDataAP5_Late <- dotPlotByGroupLate$AP5

wilcox.test(x=dotplotDataVEH_Late[,1], y=dotplotDataVEH_Late[,2], paired=T) #W= 257, p=0.000862; pcorrected=0.86340/ #Raw scores: V = 375, p-value = 0.04537
wilcox.test(x=dotplotDataAP5_Late[,1], y=dotplotDataAP5_Late[,2], paired=T) #W = 1213, p= 0.01344; pcorrected=1/ #Raw scores: V = 1172, p-value = 0.03058

p.adjust(p=c(0.64, 1.824e-05, 0.000862, 0.01344), method="holm")
#6.400e-01 7.296e-05 2.340e-03, 2.688e-01

#Raw scores
p.adjust(p=c(0.46, 0.0011, 0.045, 0.03), method="holm")
# 0.4600 0.0044 0.0900 0.0900


################################################
### EXCITATION AND INHIBITION BY BIN
###############################################
postInfTargetWdw <- 1800+12*60+30*60 #For the post infusion window, I'll choose the period between the end of the infusion +30'.

#Matrix in which rows are 50ms bins after the cue, columns are individual neurons and the values indicate if the neuron was EXCITED (ExcBins) or INHIBITED (InhBins) on that bin
NEXfiles <- "E:/Dropbox/NMDA/EXP1_Performance/Late VEH/NEX files/"
LateVEHPreInf_ExcBins <- KC.sigbins(path=NEXfiles, startt=0, endt=1800, event=1, BLwdw=2, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
LateVEHPostInf_ExcBins <- KC.sigbins(path=NEXfiles, startt=2520, endt=postInfTargetWdw, event=1, BLwdw=2, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
LateVEHPreInf_InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=0, endt=1800, event=1, BLwdw=2, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
LateVEHPostInf_InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=2520, endt=postInfTargetWdw, event=1, BLwdw=2, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)

# In neuralhist, I flagged neurons as CUE-EXCITED if they were excited (>99.9% confidence interval of a Poisson distribution given by BL firing) for 3 consecutive 10ms bins in the 500ms window after the cue. I used the 2s precue window as baseline to define my Poisson distribution.
LateVEH_ExcUnits <- unlist(allNeuronsDSLateVEHPreInf$cueexidx) #Index of cue-excited units

#Redefine NEXfiles now so that it sends the function to the AP5 files and repeat
NEXfiles <- "E:/Dropbox/NMDA/EXP1_Performance/Late AP5/NEX files/"
LateAP5PreInf_ExcBins <- KC.sigbins(path=NEXfiles, startt=0, endt=1800, event=1, BLwdw=2, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
LateAP5PostInf_ExcBins <- KC.sigbins(path=NEXfiles, startt=2520, endt=postInfTargetWdw, event=1, BLwdw=2, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
LateAP5PreInf_InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=0, endt=1800, event=1, BLwdw=2, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
LateAP5PostInf_InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=2520, endt=postInfTargetWdw, event=1, BLwdw=2, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)

# In neuralhist, I flagged neurons as CUE-EXCITED if they were excited (>99.9% confidence interval of a Poisson distribution given by BL firing) for 3 consecutive 10ms bins in the 500ms window after the cue. I used the 2s precue window as baseline to define my Poisson distribution.
LateAP5_ExcUnits <- unlist(allNeuronsDSLateAP5PreInf$cueexidx) #Index of cue-excited units


#Save these files
save(LateVEHPreInf_ExcBins, file=paste(dataForRdir, "LateVEHPreInf_ExcBins.rdat", sep=""))
save(LateVEHPostInf_ExcBins, file=paste(dataForRdir, "LateVEHPostInf_ExcBins.rdat", sep=""))
save(LateVEHPreInf_InhBins, file=paste(dataForRdir, "LateVEHPreInf_InhBins.rdat", sep=""))
save(LateVEHPostInf_InhBins, file=paste(dataForRdir, "LateVEHPostInf_InhBins.rdat", sep=""))

save(LateAP5PreInf_ExcBins, file=paste(dataForRdir, "LateAP5PreInf_ExcBins.rdat", sep=""))
save(LateAP5PostInf_ExcBins, file=paste(dataForRdir, "LateAP5PostInf_ExcBins.rdat", sep=""))
save(LateAP5PreInf_InhBins, file=paste(dataForRdir, "LateAP5PreInf_InhBins.rdat", sep=""))
save(LateAP5PostInf_InhBins, file=paste(dataForRdir, "LateAP5PostInf_InhBins.rdat", sep=""))


#Plot % bins excited/inhibited before and after infusion of VEH or AP5

#Function to calculate the percentage of units exc/inh to apply on the objects that I created with KC.sigbins.R and KC.inhib.sigbins.R
PercBins <- function(sigBinData){
        sapply(seq(1, nrow(sigBinData)), function (x){
                sum(sigBinData[x,])/ncol(sigBinData)
        })
}

#Late VEH
plot.new()
plot.window(xlim = c(0, nrow(LateVEHPreInf_ExcBins)), ylim=c(0, 1))
abline(h=seq(-1, 1, by=0.25), col="gray90")

lines(x=seq(1, nrow(LateVEHPreInf_ExcBins)), y=PercBins(LateVEHPreInf_ExcBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(LateVEHPostInf_ExcBins)), y=PercBins(LateVEHPostInf_ExcBins), col="blue", lwd=2)

lines(x=seq(1, nrow(LateVEHPreInf_InhBins)), y=-PercBins(LateVEHPreInf_InhBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(LateVEHPostInf_InhBins)), y=-PercBins(LateVEHPostInf_InhBins), col="blue", lwd=2)

axis(side=1, at=seq(0, nrow(LateVEHPreInf_ExcBins), by=10), labels=seq(0, 1, by=0.5), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, 1, by=0.5), labels=seq(0, 100, 50), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, -1, by=-0.5), labels=seq(0, 100, 50), cex.axis=1.4)
mtext(side=1, text="Time from S+ onset (s)", font=2, cex=1.5, line=2.5)
mtext(side=2, text="% Excited", at=0.5, font=2, cex = 1.5, line=2.5)
mtext(side=2, text="% Inhibited", at=-0.5, font=2, cex = 1.5, line=2.5)




#Late AP5
plot.new()
plot.window(xlim = c(0, nrow(LateAP5PreInf_ExcBins)), ylim=c(-1, 1))
abline(h=seq(-1, 1, by=0.25), col="gray90")

lines(x=seq(1, nrow(LateAP5PreInf_ExcBins)), y=PercBins(LateAP5PreInf_ExcBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(LateAP5PostInf_ExcBins)), y=PercBins(LateAP5PostInf_ExcBins), col="red", lwd=2)

lines(x=seq(1, nrow(LateAP5PreInf_InhBins)), y=-PercBins(LateAP5PreInf_InhBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(LateAP5PostInf_InhBins)), y=-PercBins(LateAP5PostInf_InhBins), col="red", lwd=2)

axis(side=1, at=seq(0, nrow(LateVEHPreInf_ExcBins), by=10), labels=seq(0, 1.5, by=0.5), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, 1, by=0.5), labels=seq(0, 100, 50), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, -1, by=-0.5), labels=seq(0, 100, 50), cex.axis=1.4)
mtext(side=1, text="Time from S+ onset (s)", font=2, cex=1.5, line=2.5)
mtext(side=2, text="% Excited", at=0.5, font=2, cex = 1.5, line=2.5)
mtext(side=2, text="% Inhibited", at=-0.5, font=2, cex = 1.5, line=2.5)



### COMPARE FR AND % CUE-EXC BINS BETWEEN EARLY VS LATE GROUPS. I'LL COMBINE THE PREINFUSION DATA OF VEH AND AP5 GROUPS
Early_PreInf_ExcBins <- cbind(EarlyVEHPreInf_ExcBins, EarlyAP5PreInf_ExcBins)
Late_PreInf_ExcBins <- cbind(LateVEHPreInf_ExcBins, LateAP5PreInf_ExcBins)

plot.new()
plot.window(xlim = c(0, 20), ylim=c(0, 1)) #Plot 1s (20 bins of 50ms)
abline(h=seq(-1, 1, by=0.25), col="gray90")

lines(x=seq(1, nrow(Early_PreInf_ExcBins)), y=PercBins(Early_PreInf_ExcBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(Late_PreInf_ExcBins)), y=PercBins(Late_PreInf_ExcBins), col="black", lwd=2)

axis(side=1, at=seq(0, nrow(LateVEHPreInf_ExcBins), by=10), labels=seq(0, 1.5, by=0.5), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, 1, by=0.5), labels=seq(0, 100, 50), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, -1, by=-0.5), labels=seq(0, 100, 50), cex.axis=1.4)
mtext(side=1, text="Time from S+ onset (s)", font=2, cex=1.5, line=2.5)
mtext(side=2, text="% Excited", at=0.5, font=2, cex = 1.5, line=2.5)
mtext(side=2, text="% Inhibited", at=-0.5, font=2, cex = 1.5, line=2.5)

legend(x=15, y=0.99, col=c("gray30", 'black'), legend=c("Early", "Late"), lty=1, lwd=2, cex=1.5, bty="n")

# % of CUE-EXCITED UNITS 
# In neuralhist, I flagged neurons as CUE-EXCITED if they were excited (>99.9% confidence interval of a Poisson distribution given by BL firing) for 3 consecutive 10ms bins in the 500ms window after the cue. I used the 2s precue window as baseline to define my Poisson distribution.
# In neuralhist, I flagged neurons as CUE-EXCITED if they were excited (>99.9% confidence interval of a Poisson distribution given by BL firing) for 3 consecutive 10ms bins in the 500ms window after the cue. I used the 2s precue window as baseline to define my Poisson distribution.
EarlyVEHPreInf_ExcUnits <- unlist(allNeuronsDSEarlyVEHPreInf$cueexidx)
EarlyAP5PreInf_ExcUnits <- unlist(allNeuronsDSEarlyAP5PreInf$cueexidx)
EarlyVEHPostInf_ExcUnits <- unlist(allNeuronsDSEarlyVEHPostInf$cueexidx)
EarlyAP5PostInf_ExcUnits <- unlist(allNeuronsDSEarlyAP5PostInf$cueexidx)

LateVEHPreInf_ExcUnits <- unlist(allNeuronsDSLateVEHPreInf$cueexidx)
LateAP5PreInf_ExcUnits <- unlist(allNeuronsDSLateAP5PreInf$cueexidx)
LateVEHPostInf_ExcUnits <- unlist(allNeuronsDSLateVEHPostInf$cueexidx)
LateAP5PostInf_ExcUnits <- unlist(allNeuronsDSLateAP5PostInf$cueexidx)




################################################
### EXCITATION AND INHIBITION BY BIN
################################################

#Late VEH
plot.new()
par(mar=c(2, 6, 2, 2))
plot.window(xlim=c(0, 2), ylim=c(0, 1))
rect(xleft=0, xright=1, ybottom=0, ytop=sum(LateVEHPreInf_ExcUnits)/length(LateVEHPreInf_ExcUnits), col="gray30", border = F)
rect(xleft=1, xright=2, ybottom=0, ytop=sum(LateVEHPostInf_ExcUnits)/length(LateVEHPostInf_ExcUnits), col=colindx[1], border= F)
axis(side=1, at=c(0.5, 1.5), tick = F, labels=c("Pre", "Post"), cex.axis=1.5, font=2)
axis(side=2, at=seq(0, 1, 0.25), cex.axis=1.4, las=2)
mtext(side=2, line=4, text="Proportion", cex=1.5, font=2)


#Late AP5
plot.new()
par(mar=c(2, 6, 2, 2))
plot.window(xlim=c(0, 2), ylim=c(0, 1))
rect(xleft=0, xright=1, ybottom=0, ytop=sum(LateAP5PreInf_ExcUnits)/length(LateAP5PreInf_ExcUnits), col="gray30", border = F)
rect(xleft=1, xright=2, ybottom=0, ytop=sum(LateAP5PostInf_ExcUnits)/length(LateAP5PostInf_ExcUnits), col=colindx[2], border= F)
axis(side=1, at=c(0.5, 1.5), tick = F, labels=c("Pre", "Post"), cex.axis=1.5, font=2)
axis(side=2, at=seq(0, 1, 0.25), cex.axis=1.4, las=2)
mtext(side=2, line=4, text="Proportion", cex=1.5, font=2)


Early_PreInf_ExcBins_ExcUnits <- cbind(EarlyVEHPreInf_ExcBins, EarlyAP5PreInf_ExcBins)
Late_PreInf_ExcBins <- cbind(LateVEHPreInf_ExcBins, LateAP5PreInf_ExcBins)


sum(Early_PreInf_ExcBins)/length(Early_PreInf_ExcBins)
sum(Late_PreInf_ExcBins)/length(Late_PreInf_ExcBins)





#Dot plot (or boxplot) that says: of the cue-excited neurons, during what % of bins were those units excited before and after injection 
LateVEHPre_cueExcOnly_ExcPerBin <- LateVEHPreInf_ExcBins[,LateVEH_ExcUnits]
LateVEHPost_cueExcOnly_ExcPerBin <- LateVEHPostInf_ExcBins[,LateVEH_ExcUnits]

LateVEHPre_ExcDotPlot <- colSums(LateVEHPre_cueExcOnly_ExcPerBin)/nrow(LateVEHPre_cueExcOnly_ExcPerBin)
LateVEHPost_ExcDotPlot <- colSums(LateVEHPost_cueExcOnly_ExcPerBin)/nrow(LateVEHPost_cueExcOnly_ExcPerBin)


LateAP5Pre_cueExcOnly_ExcPerBin <- LateAP5PreInf_ExcBins[,LateAP5_ExcUnits]
LateAP5Post_cueExcOnly_ExcPerBin <- LateAP5PostInf_ExcBins[,LateAP5_ExcUnits]

LateAP5Pre_ExcDotPlot <- colSums(LateAP5Pre_cueExcOnly_ExcPerBin)/nrow(LateAP5Pre_cueExcOnly_ExcPerBin)
LateAP5Post_ExcDotPlot <- colSums(LateAP5Post_cueExcOnly_ExcPerBin)/nrow(LateAP5Post_cueExcOnly_ExcPerBin)


#Make dotplot
makeBoxPlot <- function(data, xmin, xmax, color){
        rect(xleft=xmin, xright=xmax, ybottom=summary(data)[2], ytop=summary(data)[5], border=color, lwd=2) #IQR
        segments(x0=xmin, x1=xmax, y0=summary(data)[3], col=color, lwd=2) #Median
        segments(x0=xmin, x1=xmax, y0=summary(data)[4], col="black", lwd=2) #Mean
}

plot.new()
plot.window(xlim=c(0, 2), ylim=c(0, 0.5))

makeBoxPlot(data=LateVEHPre_ExcDotPlot, xmin=0, xmax=0.3, color=colindx[1])
makeBoxPlot(data=LateVEHPost_ExcDotPlot, xmin=0.4, xmax=0.7, color=colindx[1])

makeBoxPlot(data=LateAP5Pre_ExcDotPlot, xmin=1, xmax=1.3, color=colindx[2])
makeBoxPlot(data=LateAP5Post_ExcDotPlot, xmin=1.4, xmax=1.7, color=colindx[2])


axis(side=1, at=c(0.15, 0.55), labels =c("Pre", "Post"), cex.axis=1.5, font=2)
axis(side=1, at=c(1.15, 1.55), labels =c("Pre", "Post"), cex.axis=1.5, font=2)
axis(side=2, las=2, cex.axis=1.4)
mtext(side=2, line=3, text="% excited bins", cex=1.5, font=2)

wilcox.test(LateVEHPre_ExcDotPlot, LateVEHPost_ExcDotPlot, paired=T) #V = 57.5, p-value = 0.7769
wilcox.test(LateAP5Pre_ExcDotPlot, LateAP5Post_ExcDotPlot, paired=T) #V = 439, p-value = 0.0001806

