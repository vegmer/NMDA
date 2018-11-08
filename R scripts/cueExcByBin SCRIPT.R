cueExcByBin <- function(masterDF=masterDF_DS, neudata=allNeuronsDS, capValue=c(-160, 160),
                        trialBinSize=40, WdwStart=100, WdwEnd=400, BLforPoisson=5000){
        
        binw <- neudata$parameters$binw
        
        #Index of units on neudata
        neudata$nNeurons
        
        
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
                                
                                baselinecounts <- colSums(unitFRinCounts[ , BLStartBin:BLEndBin])
                                
                                #Excitations
                                thresholdExc = 3
                                critwin = colSums(unitFRinCounts[, TargetStartBin:TargetEndBin])         #this is the window (500 ms) in which we look for excitation
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
        
        
        
        datperbin <- do.call("rbind", lapply(seq(1, nDivisions, by=1), function(i){
                
                if(sum(correspBins==i, na.rm=T)>0){
                        
                        dataSel <- filter(masterDF, correspBins==i)
                        
                        uniqUnits <- unique(dataSel$allUnitIdx)
                        
                        byBin <- do.call("rbind", lapply(seq(1, length(uniqUnits)), function(u){
                                
                                unitData <- dataSel[dataSel$allUnitIdx==uniqUnits[u], ]
                                
                                
                                
                                
                                
                                unitFRinHz <- apply(unitData[, FRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                
                                unitFRinCounts <- (unitFRinHz*binw)/1000
                                
                                baselinecounts <- colSums(unitFRinCounts[ , BLStartBin:BLEndBin])
                                
                                #Excitations
                                thresholdExc = 3
                                critwin = colSums(unitFRinCounts[, TargetStartBin:TargetEndBin])         #this is the window (500 ms) in which we look for excitation
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
        
}