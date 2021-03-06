compareDSvsNSfromCP <- function(masterDF=list(masterDF_DS, masterDF_NS), trialBinSize=15, paired=T, event="cue", 
                                cueExcOnly=F, correctOnly=F, WdwStart=100, WdwEnd=400, capped=T, capValue=c(-120, 90), 
                                dataProcess="Zscores"){
        
        output <- lapply(seq(1, length(masterDF)), function(c){
                
                if(capped==T){
                        nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                        trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
                        
                } else {
                        nDivisions <- round(((max(masterDF[[c]]$trialfromCP, na.rm=T)-min(masterDF[[c]]$trialfromCP, na.rm=T))/trialBinSize), 0)
                        trialBins <- seq(min(masterDF[[c]]$trialfromCP, na.rm=T), max(masterDF[[c]]$trialfromCP, na.rm=T), by=trialBinSize)
                }
                
                correspBins <- findInterval(masterDF[[c]]$trialfromCP, trialBins)
                masterDF[[c]]$correspBins <- correspBins
                
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                
                
                #Define the window of interest in terms of bins
                if(event=="cue"){
                        CueBin <- masterDF[[c]]$CueBin[1]
                        TargetStartBin <- CueBin+ (WdwStart/binw)
                        TargetEndBin <- CueBin + (WdwEnd/binw)
                }
                
                #I need to figure out how to define TargetStartBin and TargetEndBin around entries
                if(event=="entry"){
                        EntryBin <- masterDF[[c]]$EntryBin[1]
                        TargetStartBin <- EntryBin + (WdwStart/binw)
                        TargetEndBin <- EntryBin + (WdwEnd/binw)
                }
                
                ### ALL TRIALS (RESPONDED TO AND MISSED)
                
                dat <- do.call("rbind", lapply(seq(1, nDivisions), function(i){
                        
                        if(sum(correspBins==i, na.rm=T)>0){
                                
                                dataSel <- filter(masterDF[[c]], correspBins==i)
                                
                                if(cueExcOnly==T){
                                        dataSel <- filter(dataSel, CueExcited==T)
                                }
                                
                                if(correctOnly==T){
                                        dataSel <- filter(dataSel, !is.na(CSplusresponse))
                                }
                                
                                #When examining the tail of the excitations, discard trials in which the animal entered the port in the window under scrutiny
                                if(event=="cue" & WdwEnd>500){
                                        dataSel <- filter(dataSel, (CueLat*1000)>=WdwEnd)
                                }
                                
                                if(event=="entry" & WdwStart<0){
                                        dataSel <- filter(dataSel, (CueLat*1000)>=abs(WdwStart))
                                }
                                
                                #Calculate FR based on the dataProcess parameter
                                if(dataProcess=="Zscores"){
                                        
                                        numericDF <- apply(dataSel[, FRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        numericDFTargetWdw <- numericDF[,TargetStartBin:TargetEndBin]
                                        BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                        BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                        Zscores <- (rowMeans(numericDFTargetWdw, na.rm=T)-BLaverage)/BLsd
                                        
                                        return(data.frame(dataSel[,-FRcols], FR=Zscores))
                                }
                                
                                if(dataProcess=="raw"){
                                        
                                        numericDF <- apply(dataSel, MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        
                                        return(data.frame(dataSel[,-FRcols], FR=numericDF))
                                        
                                }
                                
                        }
                })
                )
                
        })
        
        #Compare DS and NS by bin
        
        DSdat <- output[[1]]
        NSdat <- output[[2]]
        
        save(DSdat, file=paste(dataForRdir, "DSdat.rdat", sep=""))
        save(NSdat, file=paste(dataForRdir, "NSdat.rdat", sep=""))
        
        #I calculated this inside the previous function, so I need to rerun it
        if(capped==T){
                nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
        } else {
                nDivisions <- round(((max(masterDF[[1]]$trialfromCP, na.rm=T)-min(masterDF[[1]]$trialfromCP, na.rm=T))/trialBinSize), 0)
                trialBins <- seq(min(masterDF[[1]]$trialfromCP, na.rm=T), max(masterDF[[1]]$trialfromCP, na.rm=T), by=trialBinSize)
        }
        
        do.call("rbind", lapply(seq(1, nDivisions), function(i){
                 
              selDat_DS <- DSdat[DSdat$correspBins==i, ]
              selDat_NS <- NSdat[NSdat$correspBins==i, ]
              
              unitsDS <- unique(selDat_DS$allUnitIdx)
              unitsNS <- unique(selDat_NS$allUnitIdx)
              
              #Only analyze data from neurons that contribute to both DS and NS data on that particular bin of trials
              unitIdx <- unitsDS[unitsDS %in% unitsNS]
              
              summ <- do.call("rbind", lapply(seq(1, length(unitIdx)), function(j){
                      unitDS <- selDat_DS[selDat_DS$allUnitIdx==as.numeric(unitIdx[j]), ]
                      unitNS <- selDat_NS[selDat_NS$allUnitIdx==as.numeric(unitIdx[j]), ]
                      
                      data.frame(correspBins=i, unit=unitIdx[j], DS_FR=mean(unitDS$FR, na.rm=T), NS_FR=mean(unitNS$FR, na.rm=T))
                      
              })
              )
              
              #A few neurons give me NAs when I look at entry related data (bc there may have been no entries in that window of trials in teh session in which the units were recorded), remove
              NAidx <- (1:nrow(summ))[is.na(summ$DS_FR) | is.na(summ$NS_FR)]
              if(length(NAidx)>0){summ <- summ[-NAidx, ]} 
              
              INFidx <- (1:nrow(summ))[is.infinite(summ$DS_FR) | is.infinite(summ$NS_FR)]
              if(length(INFidx)>0){summ <- summ[-INFidx, ]} 
              
              if(mean(summ$DS_FR, na.rm=T)>mean(summ$NS_FR, na.rm=T)){alt.pick="greater"}
              if(mean(summ$DS_FR, na.rm=T)<mean(summ$NS_FR, na.rm=T)){alt.pick="less"}
              if(mean(summ$DS_FR, na.rm=T)==mean(summ$NS_FR, na.rm=T)){alt.pick="two.sided"}
              
              #In one rare case, there's one neuron whose FR was the same in both conditions. Remove
              if(sum((summ$DS_FR-summ$NS_FR)==0, na.rm=T)>=1){
                      summ <- summ[-(1:nrow(summ))[(summ$DS_FR-summ$NS_FR)==0], ]
              } else {summ <- summ}
              
              #test <- t.test(x=summ$DS_FR, y=summ$NS_FR, paired = T, alternative = alt.pick)
              #I can't do t tests because it's FR data and it violates assumptions required for the t test. I have to use wilcoxon sign rank test (for paired samples).
              test <- wilcox.test(x=summ$DS_FR, y=summ$NS_FR, paired = paired, alternative = alt.pick)
              
              data.frame(bin=paste(trialBins[i], "to", trialBins[i+1]), V=test[1]$statistic, p=test[3]$p.value, n=nrow(summ))
                       
        })
        )

}



save(compareDSvsNSfromCP, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/compareDSvsNSfromCP.r")
save(compareDSvsNSfromCP, file="E:/Dropbox/NMDA/R Functions/compareDSvsNSfromCP.r")
save(compareDSvsNSfromCP, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R Functions/compareDSvsNSfromCP.r")

