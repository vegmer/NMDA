##### PLOT FIRING RATE POST CUE AS A FUNCTION OF DISTANCE FROM CHANGE POINT

plotFRandCP <- function(cue="S+", experiment, masterDF, graphFolder=MixedGraphFolder, trialBinSize=5, 
                        WdwStart=0, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                        colindx="black", legLabels=c("VEH side", "AP5 side"), yAxMinZ=-1, yAxMaxZ=1, yAxMaxRaw=7, 
                        capped=T, capValue=c(-105, 105), cueExcOnly=F, neudata=allNeuronsDS) {
        
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        
        
        binw <- neudata$parameters$binw
        
        #IF I WANT TO PLOT MORE THAN ONE GROUP PER GRAPH, MAKE A LIST WITH THE DIFFERNT masterDF OF THE DIFFERENT GROUPS
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        if(cueExcOnly==TRUE){unitSel="CueExcOnly"} else {unitSel="all units"}
        
        
        filename=paste(graphFolder, experiment, "FR by trial from CP", "bin size", trialBinSize, unitSel, WdwStart, WdwEnd, dataProcess, trialSel, ".pdf", sep="_")
        
        pdf(file=filename)
        
        plot.new()
        
        allDat <- do.call("list", lapply(seq(1, length(masterDF)), function(c){
                
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
                CueBin <- masterDF[[c]]$CueBin[1]
                TargetStartBin <- CueBin+ (WdwStart/binw)
                TargetEndBin <- CueBin + (WdwEnd/binw)
                
                
                ### ALL TRIALS (RESPONDED TO AND MISSED)
                
                datperbin <- do.call("rbind", lapply(seq(1, nDivisions, by=1), function(i){
                        if(sum(correspBins==i, na.rm=T)>0){
                                
                                
                                dataSel <- filter(masterDF[[c]], correspBins==i)
                                
                                if(cueExcOnly==T){
                                        dataSel <- filter(dataSel, CueExcited==T)
                                }
                                
                                #When examining the tail of the excitations, discard trials in which the animal entered the port in the window under scrutiny
                                if(WdwEnd>500){
                                        dataSel <- filter(dataSel, (CueLatency*1000)>=WdwEnd)
                                }
                                
                                if(correctOnly==T){
                                        
                                        dataSel <- filter(dataSel, !is.na(dataSel$CueResponse))
                                }
                                
                                if(dataProcess=="Zscores"){
                                        plot.window(xlim=c(1, nDivisions), ylim=c(yAxMinZ, yAxMaxZ))
                                        
                                        BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                        BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                        
                                        numericDF <- apply(dataSel[, FRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        #If I only have one unit in that bin I can't treat it as a DF
                                        if(sum(correspBins==i, na.rm=T)==1){
                                                numericDFTargetWdw <- numericDF[TargetStartBin:TargetEndBin]
                                                Zscores <- (mean(numericDFTargetWdw, na.rm=T)-BLaverage)/BLsd
                                        }else {
                                                numericDFTargetWdw <- numericDF[,TargetStartBin:TargetEndBin]
                                                Zscores <- (rowMeans(numericDFTargetWdw, na.rm=T)-BLaverage)/BLsd
                                        }
                                        
                                        FRtargetPeriod <- mean(Zscores, na.rm=T)
                                        sdtargetPeriod <- sd(Zscores, na.rm=T)/sqrt(length(Zscores))
                                        
                                        #Plot the actual values of FR in selected bins of trials with respect to CP
                                        if(cue[c]=="S-"){pchsel=21}
                                        if(cue[c]=="S+"){pchsel=19}
                                        
                                        errBars(x=i, y=FRtargetPeriod, err=sdtargetPeriod, hatLength = 0, jitter=0, vertical=T)
                                        points(x=i, y=FRtargetPeriod, col=colindx[c], pch=pchsel, cex=2, bg="white")
                                        
                                        #Legend
                                        if(length(masterDF)>1){legend(x=2, y=yAxMaxZ-c/2, legend=legLabels[c], fill=colindx[c], pch=pchsel, border='white', bty = 'n')}
                                        
                                }
                                
                                if(dataProcess=="raw"){
                                        
                                        plot.window(xlim=c(1, nDivisions), ylim=c(0, yAxMaxRaw))
                                        
                                        numericDF <- apply(dataSel, MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        FRtargetPeriod <- mean(numericDF[,TargetStartBin:TargetEndBin], na.rm=T)
                                        sdtargetPeriod <- sd(numericDF[,TargetStartBin:TargetEndBin], na.rm=T)/sqrt(nrow(numericDF))
                                        #Legend
                                        if(length(masterDF)>1){legend(x=2, y=yAxMaxZ-c/2, legend=legLabels[c], col=colindx[c], pch=pchsel[c], border='white', bty = 'n')}
                                        
                                        #Plot the actual values of FR in selected bins of trials with respect to CP
                                        if(cue[c]=="S-"){pchsel=21}
                                        if(cue[c]=="S+"){pchsel=19}
                                        
                                        errBars(x=i, y=FRtargetPeriod, err=sdtargetPeriod, hatLength = 0, jitter=0, vertical=T)
                                        points(x=i, y=FRtargetPeriod, col=colindx[c], pch=pchsel, cex=2, bg="white")
                                        
                                        #Legend
                                        if(length(masterDF)>1){legend(x=2, y=yAxMaxZ-c/2, legend=legLabels[c], col=colindx[c], pch=pchsel, border='white', bty = 'n')}
                                        
                                }
                                
                                
                                if(dataProcess=="PercCueExc"){
                                        
                                        selUnits <- unique(dataSel$allUnitIdx)
                                        
                                        cueExcUnits <- as.logical(sapply(seq(1, length(selUnits)), function(m){
                                                unitSession <- filter(masterDF[[c]], masterDF[[c]]$allUnitIdx==selUnits[m])
                                                exc <- unitSession$CueExcited[1]
                                                exc
                                        }))
                                        
                                        cueInhUnits <- as.logical(sapply(seq(1, length(selUnits)), function(m){
                                                unitSession <- filter(masterDF[[c]], masterDF[[c]]$allUnitIdx==selUnits[m])
                                                inh <- unitSession$CueInhibited[1]
                                                inh
                                        }))
                                        
                                        sumCueExc <- sum(cueExcUnits)
                                        sumNotCueExc <- length(cueExcUnits)-sumCueExc
                                        sumCueInh <- sum(cueInhUnits)
                                        sumNotCueInh <- length(cueInhUnits)-sumCueInh
                                        
                                        PercCueExc <- (sumCueExc/length(cueExcUnits))*100
                                        PercCueInh <- (sumCueInh/length(cueInhUnits))*100
                                        
                                        if(PercCueExc==0){exc_borderCol=colindx[c]} else {exc_borderCol="white"}
                                        if(PercCueInh==0){inh_borderCol=colindx[c]} else {inh_borderCol="white"}
                                        
                                        #Plot
                                        plot.window(xlim=c(0, nDivisions+1), ylim=c(-100, 100))
                                        
                                        #Excitations (top plot)
                                        rect(xleft = i-1, xright=i, ybottom =0, ytop=PercCueExc, col=colindx[c], border=exc_borderCol, lwd=2)
                                        
                                        #Inhibitions (bottom plot)
                                        rect(xleft = i-1, xright=i, ybottom=-PercCueInh, ytop=0, col=colindx[c], border=inh_borderCol, lwd=2)
                                        
                                        return(data.frame(bin=i, trialBins=trialBins[i], CueEx=sumCueExc, notCueExc=sumNotCueExc,
                                                          CueInh=sumCueInh, notCueInh=sumNotCueInh))
                                        
                                }
                                
                        }
                })
                )
                
                
                
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
                if(dataProcess=="PercCueExc"){
                        mtext("% cue excited", side=2, at=40, cex=1.2, line=2.5)
                        mtext("% cue inhibited", side=2, at=-5, cex=1.2, line=2.5)
                }
                mtext("Trial from behavioral change point", side=1, cex=1.2, line=2.5)
                mtext(paste("Bin size=", trialBinSize, "trials", sep=" "), side=3, adj=0.05, line=-2, font=3)
                
                return(datperbin)
                
        })
        )
        
        #Title
        titleLabel <- paste("Post ", cue, " (", WdwStart, "-", WdwEnd,  "ms) firing rate with respect to behavioral change", sep="")
        if(correctOnly==T){subLabel <- "Correct trials only"} else {subLabel="All S+ trials"}
        mtext(titleLabel, side=3, line=2.5, cex=1.2, font=2)
        mtext(subLabel, side=3, line=1, cex=1, font=3)
        
        return(allDat)
        
        dev.off()
        
        
}



save(plotFRandCP, file="E:/Dropbox/DISSERTATION/R functions/plotFRandCP.r")
save(plotFRandCP, file="E:/Dropbox/NMDA/EXP1_Performance/R Functions/plotFRandCP.R")
save(plotFRandCP, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/plotFRandCP.R")
save(plotFRandCP, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R Functions/plotFRandCP.R")
save(plotFRandCP, file="E:/Dropbox/NMDA/R Functions/plotFRandCPs.R")