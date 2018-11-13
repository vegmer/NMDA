cueExcByBin <- function(masterDF=masterDF_DS, neudata=allNeuronsDS, capValue,
                        trialBinSize, WdwStart, WdwEnd, binw=50, BLforPoisson=5000){
        
        #Index of units on neudata (the number in masterDF was assigned in the order they appear in neudata)

        nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
        trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
        
        correspBins <- findInterval(masterDF$trialfromCP, trialBins)
        masterDF$correspBins <- correspBins
        
        FRcols <- (1:ncol(masterDF))[is.element(colnames(masterDF), 1:ncol(masterDF))]
        CueBin <- masterDF$CueBin[1]
        TargetStartBin <- CueBin+ (WdwStart/binw)
        TargetEndBin <- CueBin + (WdwEnd/binw)
        
        BLStartBin <- CueBin- (BLforPoisson/binw)
        BLEndBin <- CueBin-1
        
        datperbin <- do.call("rbind", lapply(seq(1, nDivisions, by=1), function(i){

                if(sum(correspBins==i, na.rm=T)>0){

                        dataSel <- filter(masterDF, correspBins==i)

                        uniqUnits <- unique(dataSel$allUnitIdx)

                        byBin <- do.call("rbind", lapply(seq(1, length(uniqUnits)), function(u){

                                unitData <- dataSel[dataSel$allUnitIdx==uniqUnits[u], ]
                                unitFRinHz <- apply(unitData[, FRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric

                                unitFRinCounts <- (unitFRinHz*binw)/1000
                                
                                

                                #Excitations
                                thresholdExc = 3
                                if(is.null(dim(unitFRinCounts))){
                                        baselinecounts <- unitFRinCounts[BLStartBin:BLEndBin]
                                        critwin <- unitFRinCounts[TargetStartBin:TargetEndBin]} else {
                                        baselinecounts <- colSums(unitFRinCounts[ , BLStartBin:BLEndBin])
                                        critwin = colSums(unitFRinCounts[, TargetStartBin:TargetEndBin])    
                                }
                                
                                #This is the window (500 ms) in which we look for excitation
                                critval = poisson.test(x = sum(baselinecounts), T = length(baselinecounts), conf.level = 0.999)$conf.int[2]   #computes the upper limit of the confindence interval based on the baseline

                                diffs = diff(which(critwin > critval))       #computes the differences in indices for bins exceeding the critical value (to check whether excited bins are consecutive or not)
                                cueex = F
                                if(length(which(rle(diffs)$values == 1 & rle(diffs)$lengths >= thresholdExc))>0) (cueex = T)   #looks for consecutive bins (diff equal to 1) that are at least 2 bins long -or whatever #of bins specified by the argument "threshold" (rle length of at least 2)

                                #Inhibitions
                                thresholdInh = 2
                                critvalInh = poisson.test(x = sum(baselinecounts), T = length(baselinecounts), conf.level = 0.999)$conf.int[1]   #computes the lower limit of the confindence interval based on the baseline
                                diffsInh = diff(which(critwin < critval))
                                cueinh = F
                                if(length(which(rle(diffsInh)$values == 1 & rle(diffsInh)$lengths >= thresholdInh))>0) (cueinh = T)   #looks for consecutive bins (diff equal to 1) that are at least 2 bins long -or whatever #of bins specified by the argument "threshold" (rle length of at least 2)


                                #For output
                                rat <- unitData$rat[1]
                                allUnitIdx <- uniqUnits[u]
                                correspBins <- unitData$correspBins[1]

                                data.frame(rat=rat, allUnitIdx=allUnitIdx, correspBins=correspBins, DSexc=cueex, DSinh=cueinh)

                        })
                        )

                        byBin
                }
        })
        )

        
        ####Using the data in masterlist (allNeuronsDS) just to double check that I get the same results
        # datperbin <- do.call("rbind", lapply(seq(1, nDivisions, by=1), function(i){
        #         
        #         if(sum(correspBins==i, na.rm=T)>0){
        #                 
        #                 dataSel <- filter(masterDF, correspBins==i)
        #                 
        #                 uniqUnits <- unique(dataSel$allUnitIdx)
        #                 
        #                 byBin <- do.call("rbind", lapply(seq(1, length(uniqUnits)), function(u){
        #                         
        #                         unitData <- dataSel[dataSel$allUnitIdx==uniqUnits[u], ]
        #                         unitDistr <- c(0, cumsum(neudata$nNeurons)) #Distribution of unit indexes in NEX files 
        #                         unitName <- as.numeric(as.character(uniqUnits[u]))
        # 
        #                         NeudataSess <- findInterval(unitName, unitDistr, all.inside=TRUE, left.open = TRUE)
        #                         unitPosInSess <- unitName-as.numeric(unitDistr[NeudataSess])
        #                         
        #                         masterlistUnit <- neudata$masterlist[[NeudataSess]][[unitPosInSess]]
        #                         
        #                         trialsToSelect <- unitData$trialIdx
        #                         
        #                         ML_TrialBin <- masterlistUnit[c(trialsToSelect)]
        #                         
        #                         
        #                         #Excitations
        #                         thresholdExc = 3
        #                         pbin <- binw/1000
        #                         BLwdw=BLforPoisson/1000
        #                         
        #                         allvals = unlist(ML_TrialBin)
        #                         hcounts = hist(allvals[which(allvals >= -BLwdw & allvals <= .5)], breaks = seq(-BLwdw, .5, pbin), plot = F)$counts 
        #                         
        #                         
        #                         baselinecounts=hcounts[1:(BLwdw/pbin)]
        #                         critwin = hcounts[(BLwdw/pbin + 1):((BLwdw+.5)/pbin)]         #this is the window (500 ms) in which we look for excitation
        #                         critvalexc = poisson.test(x = sum(baselinecounts), T = length(baselinecounts), conf.level = 0.999)$conf.int[2]   #computes the upper limit of the confindence interval based on the baseline
        #                         
        #                         diffs = diff(which(critwin > critvalexc))       #computes the differences in indices for bins exceeding the critical value (to check whether excited bins are consecutive or not)
        #                         cueex = F
        #                         if(length(which(rle(diffs)$values == 1 & rle(diffs)$lengths >= thresholdExc))>0) (cueex = T)   #looks for consecutive bins (diff equal to 1) that are at least 2 bins long -or whatever #of bins specified by the argument "threshold" (rle length of at least 2)
        #                         
        #                         
        #                         #Inhibitions
        #                         thresholdInh = 2
        #                         
        #                         critvalInh = poisson.test(x = sum(baselinecounts), T = length(baselinecounts), conf.level = 0.999)$conf.int[1]   #computes the lower limit of the confindence interval based on the baseline
        #                         
        #                         diffsInh = diff(which(critwin > critvalInh))       #computes the differences in indices for bins exceeding the critical value (to check whether excited bins are consecutive or not)
        #                         cueinh = F
        #                         if(length(which(rle(diffsInh)$values == 1 & rle(diffsInh)$lengths >= thresholdInh))>0) (cueinh = T)   #looks for consecutive bins (diff equal to 1) that are at least 2 bins long -or whatever #of bins specified by the argument "threshold" (rle length of at least 2)
        #                         
        #                         
        #                         #For output
        #                         rat <- unitData$rat[1]
        #                         expt <- unitData$session[1]
        #                         correspBins <- unitData$correspBins[1]
        #                         
        #                         data.frame(rat=rat, expt=expt, allUnitIdx=unitName, correspBins=correspBins, DSexc=cueex, DSinh=cueinh)
        #                         
        #                 })
        #                 )
        #                 
        #                 byBin
        #         }
        # })
        # )
        
        datperbin
        
}


save(cueExcByBin, file=paste(funcdirect, "cueExcByBin.R", sep=""))
