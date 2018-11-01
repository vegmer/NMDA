##### PLOT FIRING RATE POST CUE AS A FUNCTION OF DISTANCE FROM CHANGE POINT

plotFRandCP_Entry <- function(cue="S+", wdwLabel="Post S+", experiment, masterDF, graphFolder=MixedGraphFolder, 
                        trialBinSize=5, WdwStart, WdwEnd, dataProcess="Zscores", correctOnly=FALSE, 
                        colindx=colindx, legLabels=c("VEH side", "AP5 side"), yAxMinZ=-1, yAxMaxZ=1, yAxMaxRaw=7, 
                        capped=T, capValue=c(-100, 100), cueExcOnly=F, neudata=allNeuronsDS) {
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        
        
        binw <- neudata$parameters$binw
        
        #IF I WANT TO PLOT MORE THAN ONE GROUP PER GRAPH, MAKE A LIST WITH THE DIFFERNT masterDF OF THE DIFFERENT GROUPS
        
        filename=paste(graphFolder, experiment, "FR from CP", trialBinSize, "bins",  dataProcess, wdwLabel, WdwStart, WdwEnd, ".pdf", sep=" ")
        
        pdf(file=filename)
        
               
                plot.new()
                
                lapply(seq(1, length(masterDF)), function(c){
                        dotColIdx <- colindx
                        trialBinSize <- trialBinSize
                        if(capped==T){
                                nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                                trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
                                } else {
                                nDivisions <- round(((max(masterDF[[c]]$trialfromCP, na.rm=T)-min(masterDF[[c]]$trialfromCP, na.rm=T))/trialBinSize), 0)
                                trialBins <- seq(min(masterDF[[c]]$trialfromCP, na.rm=T), max(masterDF[[c]]$trialfromCP, na.rm=T), by=trialBinSize)
                        }
                        
                        
                        correspBins <- findInterval(masterDF[[c]]$trialfromCP, trialBins)
                        FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                        masterDF[[c]]$correspBins <- correspBins
                        EntryBin <- masterDF[[c]]$EntryBin[1]
                        TargetStartBin <- EntryBin + (WdwStart/binw)
                        TargetEndBin <- EntryBin + (WdwEnd/binw)
                        
                        ### ALL TRIALS (RESPONDED TO AND MISSED)
                        sapply(seq(1, nDivisions), function(i){
                                        if(sum(correspBins==i, na.rm=T)>0){
                                                if(dataProcess=="Zscores"){
                                                        plot.window(xlim=c(1-0.5, nDivisions), ylim=c(yAxMinZ, yAxMaxZ))
                                                        
                                                        dataSel <- filter(masterDF[[c]], correspBins==i)
                                                        
                                                        if(cueExcOnly==T){
                                                                dataSel <- filter(dataSel, CueExcited==T)
                                                        }
                                                        
                                                        #If I'm going to be examining the period before the entry, exclude those trials in which the window includes the cue
                                                        if(WdwStart<0){
                                                                dataSel <- filter(dataSel, (CueLat*1000)>=abs(WdwStart))
                                                        }
                                                        
                                                        
                                                        
                                                        numericDF <- apply(dataSel[, FRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                                        numericDFTargetWdw <- numericDF[,TargetStartBin:TargetEndBin]
                                                        BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                                        BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                                        #If BLsd=0 (it happens in a couple of units), when I calculate Zsc it gives me an Inf value. In those cases, treat BLsd as 1 (Zsc would just be the difference between the value and the mean of BL)
                                                        BLsd[BLsd==0] <- 1
                                                        Zscores <- (rowMeans(numericDFTargetWdw, na.rm=T)-BLaverage)/BLsd
                                                        FRtargetPeriod <- mean(Zscores, na.rm=T)
                                                        sdtargetPeriod <- sd(Zscores, na.rm=T)/sqrt(length(Zscores))
                                                        #Legend
                                                        if(length(masterDF)>1){legend(x=2, y=yAxMaxZ-c/2, legend=legLabels[c], fill=dotColIdx[c], border='white', bty = 'n')}
                                                        
                                                        #Plot the actual values of FR in selected bins of trials with respect to CP
                                                        errBars(x=i-0.5, y=FRtargetPeriod, err=sdtargetPeriod, hatLength = 0, jitter=0, vertical=T)
                                                        
                                                        # if(length(cue)>1){
                                                        #         cueorig <- cue
                                                        #         if(c==1){cue <- cueorig[1]}
                                                        #         if(c==2){cue <- cueorig[2]}
                                                        # }
                                                        
                                                        if(cue[c]=="S-"){points(x=i-0.5, y=FRtargetPeriod, col=dotColIdx[c], pch=21, bg="white", cex=2)}
                                                        if(cue[c]=="S+"){points(x=i-0.5, y=FRtargetPeriod, col=dotColIdx[c], pch=19, cex=2)}
                                                       
                                                }
                                                
                                                if(dataProcess=="raw"){
                                                        plot.window(xlim=c(1-0.5, nDivisions), ylim=c(0, yAxMaxRaw))
                                                        dataSel <- filter(masterDF[[c]], correspBins==i)[,FRcols]
                                                        if(cueExcOnly==T){
                                                                dataSel <- filter(dataSel, CueExcited==T)
                                                        }
                                                        
                                                        #If I'm going to be examining the period before the entry, exclude those trials in which the window includes the cue
                                                        if(WdwStart<0){
                                                                dataSel <- filter(dataSel, (CueLat*1000)>=abs(WdwStart))
                                                        }
                                                        
                                                        numericDF <- apply(dataSel, MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                                        FRtargetPeriod <- mean(numericDF[,TargetStartBin:TargetEndBin], na.rm=T)
                                                        sdtargetPeriod <- sd(numericDF[,TargetStartBin:TargetEndBin], na.rm=T)/sqrt(nrow(numericDF))
                                                        #Legend
                                                        if(length(masterDF)>1){legend(x=2, y=yAxMaxRaw-c/2, legend=legLabels[c], fill=dotColIdx[c], border='white', bty = 'n')}
                                                        
                                                        #Plot the actual values of FR in selected bins of trials with respect to CP
                                                        errBars(x=i-0.5, y=FRtargetPeriod, err=sdtargetPeriod, hatLength = 0, jitter=0, vertical=T)
                                                        
                                                        # if(length(cue)>1){
                                                        #         cueorig <- cue
                                                        #         if(c==1){cue <- cueorig[1]}
                                                        #         if(c==2){cue <- cueorig[2]}
                                                        # }
                                                        
                                                        if(cue[c]=="S-"){points(x=i-0.5, y=FRtargetPeriod, col=dotColIdx[c], pch=21, bg="white", cex=2)}
                                                        if(cue[c]=="S+"){points(x=i-0.5, y=FRtargetPeriod, col=dotColIdx[c], pch=19, cex=2)}
                                                        
                                                }
                                                
                                                       
                                                # if(dataProcess=="PercCueExc"){
                                                #         
                                                #         plot.window(xlim=c(0, nDivisions+1), ylim=c(0, 40))
                                                #         
                                                #         dataSel <- filter(masterDF[[c]], correspBins==i)
                                                #         selUnits <- unique(dataSel$allUnitIdx)
                                                #         BLwdw <- 2 #in s
                                                #         CueBin <- masterDF[[c]]$CueBin[1]
                                                #         BLwdwStart <- CueBin- (BLwdw*1000)/binw
                                                #         BLwdwEnd <- CueBin-1
                                                #         
                                                #         cueExcUnits <- sapply(seq(1, length(selUnits)), function(m){
                                                #                 #Find out the Poisson distribution of FR of that unit during BL throughout the session in which it was recorded
                                                #                 unitSession <- filter(masterDF[[c]], masterDF[[c]]$allUnitIdx==selUnits[m])
                                                #                 FRdf <- unitSession[, FRcols]
                                                #                 numericDF <- apply(FRdf, MARGIN=2, as.numeric)
                                                #                 BLfr <- mean(colSums(numericDF[, BLwdwStart:BLwdwEnd])/nrow(numericDF))
                                                #                 
                                                #                 #Find out the FR in the target period ONLY in the trials that belong to the defined bin
                                                #                 unitSelTrials <- filter(dataSel, allUnitIdx==selUnits[m])
                                                #                 FRdf <- apply(unitSelTrials[, FRcols], MARGIN=2, as.numeric)
                                                #                 if(dim(unitSelTrials)[1]==1){                             #I have to do this bc if there's only one trial for that unit, I can't subset columns
                                                #                         
                                                #                         TargetWdw <- FRdf[TargetStartBin:TargetEndBin]
                                                #                         meanFRTargetWdw <- mean(TargetWdw)
                                                #                         } else { 
                                                #                         TargetWdw <- FRdf[, TargetStartBin:TargetEndBin]
                                                #                         meanFRTargetWdw <- colMeans(TargetWdw)
                                                #                 }
                                                #                
                                                #                 #Given the Poisson distribution of FR during BL throughout the session, what's the probability that the unit is excited at different BINS in the subset of trials we're examining (in the target window)
                                                #                 probVect <- ppois(meanFRTargetWdw, BLfr)
                                                #                 
                                                #                 sigbins = rep(F, length(probVect))
                                                #                 sigbins[which(probVect >= .999)] = T #If the probability of observing the count value exceeds a 99.9% confidence interval of our Poisson distribution (defined in ppois using the 2s pre-cue window), then flag that neuron as excited on that bin
                                                #                 
                                                #                 cueex=F
                                                #                 
                                                #                 if(length(which(rle(sigbins)$values == 1 & rle(sigbins)$lengths >= threshold))>0) (cueex = T) #looks for consecutive exc bins (sigbins equal to 1) that are at least 3 bins long (rle length of at least 3)
                                                #                 
                                                #                 return(cueex)
                                                #                 
                                                #         })
                                                #         
                                                #        PercCueExc <- (sum(cueExcUnits)/length(cueExcUnits))*100
                                                #        if(PercCueExc==0){borderCol=colindx[c]} else {borderCol="white"}
                                                #        rect(xleft = i-1, xright=i, ybottom =0, ytop=PercCueExc, col=colindx[c], border=borderCol)
                                                #    
                                                # }
                                                
                                                
                                                
                                        }
                                })
                
                        #Axis
                        if(capped==T){trialfromCPvalues <- capValue[1]:capValue[2]} else {
                                trialfromCPvalues <- min(masterDF[[c]]$trialfromCP, na.rm=TRUE):max(masterDF[[c]]$trialfromCP, na.rm=TRUE)   
                        }
                        
                        trialsToLabel <- seq(-1000, 1000, by=10)
                        lngthScreenUnit <- nDivisions/length(trialfromCPvalues)
                        screenMarks <- cumsum(rep(lngthScreenUnit, length(trialfromCPvalues)))
                        atVals <- screenMarks[trialfromCPvalues %in% trialsToLabel]
                        labels <- trialfromCPvalues[trialfromCPvalues %in% trialsToLabel]
                        axis(side=1, at=atVals, labels=labels, cex.axis=1.8)
                        
                        abline(v=screenMarks[which(trialfromCPvalues==0)])
                        axis(side=2, las=2, cex.axis=1.8)
                        
                        #Labels
                        if(dataProcess=="Zscores"){mtext("Firing rate (Zsc)", side=2, cex=1.2, line=2.5)}  
                        if(dataProcess=="raw"){mtext("Firing rate", side=2, cex=1.2, line=2.5)}
                        if(dataProcess=="PercCueExc"){mtext("% cue excited", side=2, cex=1.2, line=2.5)}
                        mtext("Trial from behavioral change point", side=1, cex=1.2, line=2.5)
                        mtext(paste("Bin size=", trialBinSize, "trials", sep=" "), side=3, adj=0.05, line=-2, font=3)
                        
                        #Title
                        titleLabel <- paste(wdwLabel,  "(", WdwStart, "-", WdwEnd,  "ms) firing rate with respect to behavioral change", sep="")
                        mtext(titleLabel, side=3, line=2.5, cex=1.2, font=2)
                        mtext(subLabel, side=3, line=1, cex=1, font=3)
                        
                })
                
                dev.off()
                
      
}



save(plotFRandCP_Entry, file="E:/Dropbox/DISSERTATION/R functions/plotFRandCP_Entry.r")
save(plotFRandCP_Entry, file="E:/Dropbox/NMDA/EXP1_Performance/R Functions/plotFRandCP_Entry.R")
save(plotFRandCP_Entry, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/plotFRandCP_Entry.R")
save(plotFRandCP_Entry, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R Functions/plotFRandCP_Entry.R")
save(plotFRandCP_Entry, file="E:/Dropbox/NMDA/R Functions/plotFRandCP_Entry.R")


        
