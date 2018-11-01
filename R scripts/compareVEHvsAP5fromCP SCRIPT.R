compareVEHvsAP5fromCP <- function(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), trialBinSize=15, 
                                  paired=T, event="cue", correctOnly=F, WdwStart=0, WdwEnd=400, cueExcOnly=F,
                                  capped=T, capValue=c(-90, 90), dataProcess="Zscores"){
        
        if(capped==T){
                nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
        } else {
                nDivisions <- round(((max(masterDF[[1]]$trialfromCP, na.rm=T)-min(masterDF[[1]]$trialfromCP, na.rm=T))/trialBinSize), 0)
                trialBins <- seq(min(masterDF[[1]]$trialfromCP, na.rm=T), max(masterDF[[1]]$trialfromCP, na.rm=T), by=trialBinSize)
        }
        
        output <- lapply(seq(1, length(masterDF)), function(c){
                
                
                
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
                                        dataSel <- filter(dataSel, (CueLat*1000)>=WdwEnd)
                                }
                                
                                if(event=="entry" & WdwStart<0){
                                        dataSel <- filter(dataSel, (CueLat*1000)<=abs(WdwStart))
                                }
                                
                                #Calculate FR based on the dataProcess parameter
                                if(dataProcess=="raw"){
                                        
                                        numericDF <- apply(dataSel[, FRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        
                                        return(data.frame(dataSel[,-FRcols], FR=numericDF))
                                        
                                }
                                if(dataProcess=="Zscores"){
                                        
                                        numericDF <- apply(dataSel[, FRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        numericDFTargetWdw <- numericDF[,TargetStartBin:TargetEndBin]
                                        BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                        BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                        Zscores <- (rowMeans(numericDFTargetWdw, na.rm=T)-BLaverage)/BLsd
                                        
                                        w <- data.frame(dataSel[,-FRcols], FR=Zscores)
                                        
                                        # if(adjust==T){
                                        #         w <- w[!is.nan(w$FR), ] #When FR is NaN, it's bc for that rat I have units on the other side but not on this side. Get rid of that
                                        #         ratsInBin <- unique(w$rat)
                                        #         sapply(seq(1, length(ratsInBin)), function(r){
                                        #                 ratW <- w[w$rat==ratsInBin[r], ] #All the units recorded for that rat on this bin. Now I want to test whether they belong to the same session or not.
                                        #                 uniqSess <- unique(ratW$session)
                                        #                 if(length(uniqSess)>1){
                                        #                         table(ratW$session, ratW$trialfromCP)
                                        #                         idxMost <- (rle(ratW$session)$length)==max(rle(ratW$session)$length)
                                        #                         sessMost <- uniqSess[idxMost] #Identify the session from which more trials come in this bin. Choose that session
                                        #                         ratW[ratW$session==sessMost, ]
                                        # 
                                        #                 }
                                        #         })
                                        # }

                                        return(w)
                                }
                                
                        }
                })
                )
                
        })
        
        #Compare DS and NS by bin
        
        VEHdat <- output[[1]]
        AP5dat <- output[[2]]
        
        save(VEHdat, file=paste(dataForRdir, "VEHdat.rdat", sep=""))
        save(AP5dat, file=paste(dataForRdir, "AP5dat.rdat", sep=""))
        
        #I calculated this inside the previous function, so I need to redefine it
        nDivisions <- length(unique(VEHdat$correspBins))
        
        sig.test <- do.call("rbind", lapply(seq(1, nDivisions), function(i){
                
                selDat_VEH <- VEHdat[VEHdat$correspBins==i, ]
                selDat_AP5 <- AP5dat[AP5dat$correspBins==i, ]
                
                unitsVEH <- unique(selDat_VEH$allUnitIdx)
                unitsAP5 <- unique(selDat_AP5$allUnitIdx)
                
                ratsVEH <- unique(selDat_VEH$rat)
                ratsAP5 <- unique(selDat_AP5$rat)
                
                #If I only want to compare VEH vs. AP5 side within individual rats
                if(paired==TRUE){
                        #Only analyze data from rats that contribute to both VEH and AP5 data on that particular bin of trials
                        ratsIdx <- ratsVEH[ratsVEH %in% ratsAP5]
                        
                        if(length(ratsIdx)==0){wilcox.test.results<- data.frame(bin=NA, V=NA, p=NA)}
                        
                        if(length(ratsIdx)!=0){
                                
                                #"summ" is a data frame where each row is FR of neurons in VEH side (VEH_FR) or AP5 side (AP5_FR) of a single animal in each one of the trials in the bin under scrutiny (bin). Also, the number of units on each side (VEH_units and AP5_units)
                                summ <- do.call("rbind", lapply(seq(1, length(ratsIdx)), function(j){
                                        ratVEH <- selDat_VEH[selDat_VEH$rat==ratsIdx[j], ]
                                        ratAP5 <- selDat_AP5[selDat_AP5$rat==ratsIdx[j], ]
                                        
                                        NunitsVEH <- length(unique(ratVEH$allUnitIdx))
                                        NunitsAP5 <- length(unique(ratAP5$allUnitIdx))
                                        
                                        uniqueTrialFromCP <- unique(ratVEH$trialfromCP)
                                        
                                        VEHfrPerTrial <- sapply(seq(1, length(uniqueTrialFromCP)), function(t){
                                                selTrial <- ratVEH[ratVEH$trialfromCP==uniqueTrialFromCP[t], ]
                                                mean(selTrial$FR, na.rm=T)
                                        })
                                        
                                        AP5frPerTrial <- sapply(seq(1, length(uniqueTrialFromCP)), function(t){
                                                selTrial <- ratAP5[ratAP5$trialfromCP==uniqueTrialFromCP[t], ]
                                                mean(selTrial$FR, na.rm=T)
                                        })
                                        
                                        data.frame(correspBins=i, rat=ratsIdx[j], VEH_FR=VEHfrPerTrial, AP5_FR=AP5frPerTrial, 
                                                   VEH_units=NunitsVEH, AP5_units=NunitsAP5)
                                        
                                })
                                )
                                
                                #I get NAs in some cases because there are rats that didn't have units on one of the sides on that particular set of trials
                                NAidx <- (1:nrow(summ))[is.na(summ$VEH_FR) | is.na(summ$AP5_FR)]
                                if(length(NAidx)>0){summ <- summ[-NAidx, ]} 
                                
                                
                                
                                if(mean(summ$VEH_FR, na.rm=T)>mean(summ$AP5_FR, na.rm=T)){alt.pick="greater"}
                                if(mean(summ$VEH_FR, na.rm=T)<mean(summ$AP5_FR, na.rm=T)){alt.pick="less"}
                                if(mean(summ$VEH_FR, na.rm=T)==mean(summ$AP5_FR, na.rm=T)){alt.pick="two.sided"}
                                
                                #In one rare case, there's one neuron whose FR was the same in both conditions. Remove
                                if(sum((summ$VEH_FR-summ$AP5_FR)==0)>=1){
                                        summ <- summ[-(1:nrow(summ))[(summ$VEH_FR-summ$AP5_FR)==0], ]
                                } else {summ <- summ}
                                
                                #test <- t.test(x=summ$DS_FR, y=summ$NS_FR, paired = T, alternative = alt.pick)
                                #I can't do t tests because it's FR data! I have to use wilcoxon sign rank test (for paired samples).
                                test <- wilcox.test(x=summ$VEH_FR, y=summ$AP5_FR, paired = paired, alternative = alt.pick)
                                
                                wilcox.test.results <- data.frame(bin=paste(trialBins[i], "to", trialBins[i+1]), V=test[1]$statistic, p=test[3]$p.value)
                        }
                        
                        wilcox.test.results
                        
                }
                
                #If I want to compare VEH vs. AP5 sides regardless of what rat the data came from
                if(paired==FALSE){
                        
                        VEH <- sapply(seq(1, length(unitsVEH)), function(u){
                                unitSel <- selDat_VEH[selDat_VEH$allUnitIdx==unitsVEH[u], ]
                                mean(unitSel$FR, na.rm=T)
                        })
                        
                        AP5 <- sapply(seq(1, length(unitsAP5)), function(u){
                                unitSel <- selDat_AP5[selDat_AP5$allUnitIdx==unitsAP5[u], ]
                                mean(unitSel$FR, na.rm=T)
                        })
                        
                        if(mean(VEH, na.rm=T) > mean(AP5, na.rm = T)){alt.pick="greater"}
                        if(mean(VEH, na.rm=T) < mean(AP5, na.rm=T)){alt.pick="less"}
                        if(mean(VEH, na.rm=T) == mean(AP5, na.rm=T)){alt.pick="two.sided"}
                        
                        test <- wilcox.test(x=VEH, y=AP5, paired=FALSE, alternative = alt.pick)  
                        
                        wilcox.test.results<- data.frame(bin=paste(trialBins[i], "to", trialBins[i+1]), 
                                   n=paste(length(VEH), " vs. ", length(AP5), sep=""), 
                                   W=test[1]$statistic, p=test[3]$p.value)
                        
                        wilcox.test.results
                }
                
                return(wilcox.test.results)
                           
                }) 
                )
        
        sig.test
        
}



save(compareVEHvsAP5fromCP, file="E:/Dropbox/NMDA/R Functions/compareVEHvsAP5fromCP.R")
