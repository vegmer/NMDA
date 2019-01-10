FRbyNEURONbyBINcue <- function(neudata, CPvector=CPdata$CP, sessionCPperRat=sessionCPperRat, funcdirect=funcdirect, SpikeEnd=0.4, dataForRdir=dataForRdir, BLduration=2, cueExcOnly=F, cueKind=1){
        
        #cueKind=1 is DS; cueKind=2 is NS
        
        #rats <- unique(sapply(seq(1, length(allNeuronsDSVEH$nexdata)), function(x){allNeuronsDSVEH$nexdata[[x]]$ratname}))
        ### DEFINE OBJECTS THAT WILL BE USEFUL LATER ON
        binw <- neudata$parameters$binw
        psthmin <- neudata$parameters$psthmin
        psthmax <- neudata$parameters$psthmax
        binIdx <- seq((-psthmin), psthmax, by=binw/1000)
        binNo <- 1:length(binIdx)
        CueBin <- (psthmin/(binw/1000))+1
        
        #Calculate BL firing rate and sd for each unit (I used avg and sd firing rate of that unit on a window before the cue throughout the session)
        load(file=paste(funcdirect, file='calculateBLfiring.r', sep=""))
        BLdata <- calculateBLfiring(data=neudata, dataForRdir = dataForRdir, BLduration=BLduration)
        uniqSessIdx <- paste(csacqidx$subject, csacqidx$session)
        
        zscore <- function (x, avg, sd){(x-avg)/sd}
        
        neuronIdx <- c()
        for(i in 1:length(neudata$neurons)){
                a <- rep(i, length(neudata$neurons[[i]]))
                neuronIdx <- c(neuronIdx, a)
        }
        uniqNeuronIdx <- 1:length(neuronIdx)
        
      
        ### This function makes the actual master data frame
        masterDF <- do.call(rbind, lapply(seq(1, length(neudata$masterlist)), function(x){
                neuData <- neudata$masterlist[[x]] #Object with firing of each neuron on each trial w respect to the event (cue)
                nexratname <- neudata$nexdata[[x]]$ratname
                ratIdx <- rats==nexratname
                #The function I used to extract NEX data, transforms any 3digit rat name into a 2 digit one. I have a bunch of rats numbered above 100, so add a 1 in front of their NEX-extracted name. Others are named, for ex. MV06 in NEX but MV6 in MedPC. Fix that.
                if(sum(ratIdx)==0){
                        
                        if(grepl("MV0", nexratname)==T){ratIdx <- rats==gsub("MV0", "MV", nexratname)}
                        else {ratIdx <- rats==gsub("MV", "MV1", nexratname)}
                }
                sessNo <- as.numeric(neudata$nexdata[[x]]$expt); if(is.na(sessNo)){sessNo <- 1}
                sessAbsIdx <- (1:length(uniqSessIdx))[uniqSessIdx==paste(rats[ratIdx], sessNo)]
                
                behData <- alldata[[sessAbsIdx]]
                
                #Create objects for processing trial info
                kindOfTrial <- behData$orderCues
                trialAll <- (1:length(kindOfTrial))[kindOfTrial==cueKind] #Number of trial taking into account all kinds of cues
                trialIdx <- 1:length(trialAll) #Number of trial counting only that kind of cue
                cumSumTrl <- cumsum(sapply((1:sessNo), function(p){length(alldata[[idx[ratIdx][[1]][p]]]$CSpluscue)})) #Cumulative sum of trials of that kind throughout the sessions
                if(sessNo==1){trialCum <- trialIdx} else {trialCum <- trialIdx+cumSumTrl[sessNo-1]} #Number of trial of that kind of cue taking into account previous sessions
                
                ratsCP <- CPvector[rats==behData$ratname]
                ratsCPsess <- sessionCPperRat[rats==behData$ratname]
                
                #When the 1st trial comes after 5s and it's a long cue, I discarded trial for the Cue specificity calculation. Replace now with an NA if that's the case
                CSplusSpecif=behData$CSplusTA
                if(length(CSplusSpecif)<length(behData$CSpluscue)){CSplusSpecif <- c(NA, CSplusSpecif)}
                
                #Data frame with behavioral data for each trial
                BEH <- data.frame(SessidxNEX=x,
                                  SessidxBEH=sessAbsIdx,
                                  rat=behData$ratname,
                                  session=sessNo,
                                  group=behData$group,
                                  trialAll=trialAll,
                                  trialIdx=trialIdx,
                                  trialCum=trialCum,
                                  trialfromCP=trialCum-ratsCP,
                                  sessfromCPsess=sessNo-ratsCPsess,
                                  CSplusTime=behData$CSpluscue,
                                  CSplusResponse=behData$CSplusresponse,
                                  CSplusLat=behData$CSplusLat,
                                  CSplusSpecif=CSplusSpecif,
                                  Cuelength=behData$CueLength[kindOfTrial==1],
                                  ITIlatency=behData$ITIlatency[kindOfTrial==1],
                                  CueBin= CueBin, #First bin after cue onset
                                  PostSpikeBin=CueBin+SpikeEnd/(binw/1000), #Bin that demarcates the end of my cue-evoked ROI for the spike
                                  FiveSecBin= CueBin+5/(binw/1000), #Bin 5s after cue onset (for those rats that don't make an entry before that)
                                  CSplusEndBin= round((CueBin+(behData$CSplusLat/(binw/1000))), 0) #Bin in which the cue ended (by entry or miss)
                                  
                )
                
                BEH$uniqTrial <- paste(BEH$SessidxBEH, BEH$trialIdx)
                
                #List with firing frequency (Hz) of each unit per bin around time of event (cue onset)
                freqAllUnits <- lapply(seq(1, length(neuData)), function(y){ #For each unit on this session
                        lapply(seq(1, length(neuData[[y]])), function(z){ #For each trial
                                corresp <- findInterval(neuData[[y]][[z]], binIdx)  #Correspondencia spikes to bins
                                sapply(binNo, function(a){
                                        sum(binNo[a]==corresp)/(binw/1000)
                                })
                        })
                })
                
                
                #Index indicating if units were cue excited or not (criteria for flagging as cue-excidted is defined in the script neurohist.r)
                cueExIdx <- unlist(neudata$cueexidx[[x]])   
                BLavg <- BLdata[neuronIdx==x,1]
                BLsd <- BLdata[neuronIdx==x,2]
                
                
                freqDF <- do.call(rbind, lapply(seq(1, length(freqAllUnits)), function(y){
                        cueex <- cueExIdx[y]
                        allNeuronIdx <- uniqNeuronIdx[neuronIdx==x]
                        do.call(rbind, lapply(seq(1, length(freqAllUnits[[y]])), function(z){
                                UnitTrlFR <- freqAllUnits[[y]][[z]]
                                c(paste(sessAbsIdx, z), x, z, allNeuronIdx[y], y, cueex, round(BLavg[y], 3), round(BLsd[y], 3), UnitTrlFR)
                                
                        }))
                }))
                
                freqDF <- as.data.frame(freqDF); colnames(freqDF) <- c("uniqTrial", "SessIdxNEX", "TrialIdx", "allUnitIdx", "UnitIdx", "CueExcited", "BLavg", "BLsd",  seq(1, length(binIdx)))
                
                
                #Expand behavioral data frame to include one row per neuron
                BehAndFR <- merge(BEH, freqDF, by="uniqTrial")
                
                if(cueExcOnly==T){
                        cueExUnitsIdx <- as.logical(BehAndFR$CueExcited)
                        BehAndFR <- BehAndFR[cueExUnitsIdx,]
                } else {BehAndFR}
                
                return(BehAndFR)
        }))
        
        save(masterDF, file=paste(dataForRdir, 'masterDF.rdat', sep=""))
        return(masterDF)
}




save(FRbyNEURONbyBINcue, file="E:/Dropbox/NMDA/R functions/FRbyNEURONbyBINcue.R")
save(FRbyNEURONbyBINcue, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R functions/FRbyNEURONbyBINcue.R")

