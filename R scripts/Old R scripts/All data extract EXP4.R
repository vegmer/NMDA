

#############################################################
### EXPERIMENT 3: NO INFUSIONS, SWITCH SHORT TO LONG ITI  ###
#############################################################


########### BEHAVIORAL DATA ##############################################################
### EXTRACT BASIC BEHAVIOR OBJECTS
funcdirect <- "E:/Dropbox/NMDA/R functions/"
datafolder <- "E:/Dropbox/NMDA/EXP3_NAc FR acquisition/MedPC files/"
dataForRdir <- "E:/Dropbox/NMDA/EXP3_NAc FR acquisition/Data for R/"
dataForRCumulative <- "E:/Dropbox/NMDA/EXP3_NAc FR acquisition/Data for R cumulative/"
behGraphFolder <- "E:/Dropbox/NMDA/EXP3_NAc FR acquisition/Graphs/Behavior/"
CPfuncFolder <- paste(funcdirect, 'Change_Point-master/', sep="")
CPGraphFolder <- "E:/Dropbox/NMDA/EXP3_NAc FR acquisition/Graphs/Behavior/Change point/"
load(file=paste(funcdirect, "MedPCextract.r", sep=""))
load(file=paste(funcdirect, "mpcextract_blockSingle.Rfunc", sep=""))
load(file=paste(funcdirect, "CPextract.R", sep=""))
load(file=paste(funcdirect, "cumulativeIndGraphs.R", sep=""))
load(file=paste(funcdirect, "PerformanceFromCP.r", sep=""))


#### GENERATE IMPORTANT DATA OBJECTS
MedPCextract(funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative, cuelength=10)

# Load important behavior-related objects
files <- paste(dataForRdir, list.files(dataForRdir), sep="")
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep="")
for(i in 1:length(files)){load(files[[i]])}
for(i in 1:length(filesCum)){load(filesCum[[i]])}

##### THIS IS JUST TO RUN THE HIDDEN MARKOV MODEL, WHICH I DECIDED NOT TO DO
#### SAVE IN MATLAB'S FRIENDLY FORMAT: save one object per rat as csv. Column 1: binary measure (resp. ratio); Column 2: continuous measure (DStaskAcc inverted and positive, this is, smaller numbers indicate faster latency during cue than during ITI and I added 10 to make all numbers positive)
# lapply(seq(1, length(rats)), function(x){
#         filename <- paste(dataForRCumulative, paste("CumDat", rats[[x]], ".csv", sep=""), sep="")
#         invPosDStaskAcc <- (DStaskAcc[[x]]*-1)+10 #Transform DStaskAcc into its inverse (small numbers are better than high numbers)
#         write.csv(cbind(DSrespAll[[x]], invPosDStaskAcc), file=filename)
#         
# })


#### BEHAVIORAL GRAPHS

# RASTERS
load(file=paste(funcdirect, "MAKERASTER.r", sep=""))
MAKERASTER(i=6, data=alldata, idxdata=csacqidx)
MAKERASTER(i=4, data=alldata, idxdata=csacqidx)

#AVERAGE PERFORMANCE S+ VS S-. 5 trial bins.



### CHANGE POINT /R_code/
CPdata <- CPextract(GallCrit=1.3, minSlope=0, plot=F, CPfuncFolder=CPfuncFolder, idx=idx, CPGraphFolder=CPGraphFolder, dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative)

CPdataExp3 <- CPdata
alldataExp3 <- alldata
csacqidxExp3 <- csacqidx
ratsExp3 <- rats
idxExp3 <- idx

save(csacqidxExp3, file=paste(dataForRdir, "csacqidxExp3.ridx", sep=""))
save(alldataExp3, file=paste(dataForRdir, "alldataExp3.rdat", sep=""))
save(ratsExp3, file=paste(dataForRdir, "ratsExp3.rdat", sep=""))
save(idxExp3, file=paste(dataForRdir, "idxExp3.rdat", sep=""))

#CUMULATIVE INDIVIDUAL PERFORMANCE
IndcumPerfGraphFolder <- "E:/Dropbox/NMDA/EXP3_NAc FR acquisition/Graphs/Behavior/Cumulative perf/"
cumulativeIndGraphs(numSess=7, sessLines=F, dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf")

#MEAN PERFORMANCE AS A FUNCTION OF DISTANCE FROM CP
PerfRelToCPFolder <- "E:/Dropbox/NMDA/EXP3_NAc FR acquisition/Graphs/Behavior/Beh rel to CP/"
PerformanceFromCP(relTrialMin=-100, relTrialMax=100, trialBinSize=5, typegraph="both", imgFormat="pdf", 
                  dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=PerfRelToCPFolder, CPdata=CPdata, 
                  idx=idx, rats = rats, alldata = alldata)

#PRE CP VS. POST CP PERFORMANCE

# S+ SPECIFICITY
plot.new()
plot.window(xlim=c(0, 0.5), ylim=c(-1.5, 4))

rect(xleft=0, xright=0.25, ybottom=0, ytop=mean(CPdata$slopePre, na.rm=T), border = "white", col="gray60")
rect(xleft=0.25, xright=0.5, ybottom=0, ytop=mean(CPdata$slopePost, na.rm=T), border = "white", col="gray60")

sapply(seq(1, nrow(CPdata)), function(x){lines(x=c(0.125, 0.375), y=c(CPdata$slopePre[x], CPdata$slopePost[x]), lwd=2, col="gray20")})

abline(h=0, lty=3)
axis(side=1, at=c(0.125, 0.375), labels=c("Before CP", "After CP"), cex.axis=1.4, tick = FALSE)
axis(side=2, at=seq(-1, 4, 1), las=2, cex.axis=1.4)
mtext(side=2, line=2.5, text="S+ specificity", cex=1.5, font=2)

t.test(x=CPdata$slopePre, y=CPdata$slopePost, paired=T) #t(6)=-8.115, p=0.0001879. No correction.

PrevsPostFolder <- paste(PerfRelToCPFolder, "Average Pre vs Post/", sep="")
PrePostCP_Perf(data=DSrespAll, CP, y_axis_label="S+ response ratio", graphFolder=PrevsPostFolder, plot=T) #t(6)=-5.3517, p=0.001742
PrePostCP_Perf(data=NSrespAll, CP, y_axis_label="S- response ratio", graphFolder=PrevsPostFolder, plot=T) #t(6)= -0.36,  p=0.73
PrePostCP_Perf(data=DSlatency, CP, y_axis_label="S+ latency", graphFolder=PrevsPostFolder, plot=T) #t(6)=6.4922, p=0.0006353
PrePostCP_Perf(data=NSlatency, CP, y_axis_label="S- latency", graphFolder=PrevsPostFolder, plot=T) #t(6)= -0.049108,  p=0.9624
PrePostCP_Perf(data=DStaskAcc, CP, y_axis_label="S+ specificity", graphFolder=PrevsPostFolder, plot=T) #t(6)=-7.927, p=0.0002141
PrePostCP_Perf(data=NStaskAcc, CP, y_axis_label="S- specificity", graphFolder=PrevsPostFolder, plot=T) #t(6)= 0.88921,  p=0.4081
PrePostCP_Perf(data=ITIlatency, CP, y_axis_label="ITI latency", graphFolder=PrevsPostFolder, plot=T) #t(6)= 1.3057,  p=0.2395



############ NEURONAL DATA #################################################3
### Before running the function, define where the NEX files are
NEXfiles <- paste("E:/Dropbox/DISSERTATION/Experiment 4/Neuronal data/NEX files/", sep="")

### Load necessary functions
load(file = paste(funcdirect, "neuralhist.r", sep=""))

allNeuronsDS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T)

allNeuronsDSresponded = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T)
allNeuronsNSresponded = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T)

allNeuronsDSmissed = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T)
allNeuronsNSmissed = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T)

#sAVE THESE OBJECTS
save(allNeuronsDS, file=paste(dataForRdir, 'allNeuronsDS.rdat', sep=""))
save(allNeuronsNS, file=paste(dataForRdir, 'allNeuronsNS.rdat', sep=""))
save(allNeuronsDSresponded, file=paste(dataForRdir, 'allNeuronsDSresponded.rdat', sep=""))
save(allNeuronsNSresponded, file=paste(dataForRdir, 'allNeuronsNSresponded.rdat', sep=""))
save(allNeuronsDSmissed, file=paste(dataForRdir, 'allNeuronsDSmissed.rdat', sep=""))
save(allNeuronsNSmissed, file=paste(dataForRdir, 'allNeuronsNSmissed.rdat', sep=""))

#Load files again      
files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}

#BUILD A DATAFRAME WITH FIRING RATE OF EACH NEURON ON EACH TRIAL WITH BEH INFO ABOUT EACH TRIAL
load(file=paste(funcdirect, "FRbyNEURONbyBINcue.r", sep=""))
load(file=paste(funcdirect, "plotFRandCP.r", sep=""))
load(file=paste(funcdirect, "errBars.r", sep=""))
load(file=paste(funcdirect, "errCloud.r", sep=""))
load(file=paste(funcdirect, "plotFRandCPhistogram.r", sep=""))
load(file=paste(funcdirect, "masterDFsummary.r", sep=""))

masterDF <- FRbyNEURONbyBINcue(neudata=allNeuronsDS, funcdirect=funcdirect, dataForRdir=dataForRdir, format="Zsc", BLduration=psthmin)

masterDFsumm <- masterDFsummary(masterDF=masterDF, neudata=allNeuronsDS, format="Zscores", postCueWdw=c(0, 0.4), postSpikeWdw=c(0.4, 5), orEntry=T, cuelength=10)


### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point
MixedGraphFolder <- 'E:/Dropbox/DISSERTATION/Experiment 4/Graphs/Mixed graphs/'

plotFRandCP(experiment="Exp 4", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="raw", correctOnly=FALSE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 7, yAxMaxRaw = 10, neudata=allNeuronsDS)
plotFRandCP(experiment="Exp 4", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="Zscores", correctOnly=FALSE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 7, yAxMaxRaw = 10, neudata=allNeuronsDS)
plotFRandCP(experiment="Exp 4", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="raw", correctOnly=TRUE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 7, yAxMaxRaw = 10, neudata=allNeuronsDS)
plotFRandCP(experiment="Exp 4", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, postCueWdw=400, dataProcess="Zscores", correctOnly=TRUE, color="black", capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 7, yAxMaxRaw = 10, neudata=allNeuronsDS)


### PLOT FR IN THE POSTCUE WINDOW AND TAIL PERIOD BY RAT BY TRIAL BY NEURON
PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostCue", threshold=2, MixedGraphFolder=MixedGraphFolder)
PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostSpike", threshold=2, MixedGraphFolder=MixedGraphFolder)


### PLOT FR AROUND THE CUE (PTSH) as a function of distance to change point
MixedGraphFolder <- 'E:/Dropbox/DISSERTATION/Experiment 4/Graphs/Mixed graphs/'

plotFRandCPhistogram(experiment="Exp 4", masterDF=list(masterDF), graphFolder=MixedGraphFolder, trialBinSize=5, dataProcess="raw", correctOnly=TRUE, color="black", capped=T, capValue = c(-30, 30), yAxMinZ = -1, yAxMaxZ = 7, yAxMaxRaw = 10, imgFormat="png")

   





#COMPARE WITH OTHER GROUPS
compareCPs(data=list(CPdataExp4, CPdataExp10tr, CPdataHYBunilAP5, CPdataHYBbilAP5) , imgFormat="png", expNames=c("Task1", "Task2", "UnilAP5", "BilAP5"), colindx=c("gray30", "gray45", "gray60", "gray75"), graphFolder=behGraphFolder, minSess=5, graphs=c("Change point", "Non learners", "Dynamic Interval", "Asymptote"))


#RAW BEH GRAPH BY ANIMAL
