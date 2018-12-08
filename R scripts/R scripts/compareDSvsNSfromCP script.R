compareDSvsNSfromCP <- function(masterDF=list(masterDF_DS, masterDF_NS), trialBinSize=15, event="cue", correctOnly=F, WdwStart=0, WdwEnd=400, capped=T, capValue=c(-90, 90), dataProcess="Zscores"){
        
        output <- lapply(seq(1, length(masterDF)), function(c){
                
                if(capped==T){
                        nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                        trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
                } else {
                        nDivisions <- round(((max(masterDF[[c]]$trialfromCP, na.rm=T)-min(masterDF[[c]]$trialfromCP, na.rm=T))/trialBinSize), 0)
                        trialBins <- seq(min(masterDF[[c]]$trialfromCP, na.rm=T), max(masterDF[[c]]$trialfromCP, na.rm=T), by=trialBinSize)
                }
                
                #Assign trial to the right bin with respect to CP
                correspBins <- findInterval(masterDF[[c]]$trialfromCP, trialBins)
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                masterDF[[c]]$correspBins <- correspBins
                
                
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
                                        dataSel <- filter(dataSel, (CueLatency*1000)>=WdwEnd)
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
              
              if(mean(summ$DS_FR, na.rm=T)>mean(summ$NS_FR, na.rm=T)){alt.pick="greater"}
              if(mean(summ$DS_FR, na.rm=T)<mean(summ$NS_FR, na.rm=T)){alt.pick="less"}
              if(mean(summ$DS_FR, na.rm=T)==mean(summ$NS_FR, na.rm=T)){alt.pick="two.sided"}
              
              test <- t.test(x=summ$DS_FR, y=summ$NS_FR, paired = T, alternative = alt.pick)
              
              data.frame(bin=paste(trialBins[i], "to", trialBins[i+1]), t=test[1]$statistic, df=test[2]$parameter, p=test[3]$p.value)
                       
        })
        )

}



save(compareDSvsNSfromCP, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/compareDSvsNSfromCP.r")
save(compareDSvsNSfromCP, file="E:/Dropbox/NMDA/R Functions/compareDSvsNSfromCP.r")

