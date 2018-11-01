

plotFRandCPhistogram_Entry <- function(experiment="Exp 3", masterDF=list(masterDF_DS, masterDF_NS), 
                                       graphFolder=MixedGraphFolder, trialBinSize=15, dataProcess="Zscores", 
                                       correctOnly=FALSE, color="black", capped=T, capValue = c(-90, 90), 
                                       yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, WdwStart=-2, WdwEnd=5, 
                                       imgFormat="pdf", neudata=allNeuronsDS, eventName="Entry"){
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        if(!require(matrixStats)){install.packages("matrixStats")}
        library(matrixStats)
        
        
        if(imgFormat=="pdf"){pdf(file=paste(graphFolder, experiment, "FR from CP PSTH", "bin size", trialBinSize, dataProcess, WdwStart, WdwEnd , ".pdf", sep=" "))}
        if(imgFormat=="png"){png(filename=paste(graphFolder, experiment, "FR by trial from CP PSTH", "bin size", trialBinSize, dataProcess, WdwStart, WdwEnd, ".png", sep=" "))}
        
        binw <- neudata$parameters$binw
        
        plot.new()
        
        ZscoreCalc <- function(x, avg, sd){(x-avg)/sd}
        
        
        sapply(seq(1, length(masterDF)), function(c){
                
                if(capped==T){
                        nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                        trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
                } else {
                        nDivisions <- round(((max(masterDF[[c]]$trialfromCP, na.rm=T)-min(masterDF[[c]]$trialfromCP, na.rm=T))/trialBinSize), 0)
                        trialBins <- seq(min(masterDF[[c]]$trialfromCP, na.rm=T), max(masterDF[[c]]$trialfromCP, na.rm=T), by=trialBinSize)
                }
                correspBins <- findInterval(masterDF[[c]]$trialfromCP, trialBins)
                
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                psthWins <- ((WdwStart)/binw):((WdwEnd)/binw)
                psthWinMin <- (WdwStart)/binw
                psthWinMax <- (WdwEnd)/binw
                subFRcolNames <- (unique(masterDF[[c]]$EntryBin)+psthWinMin):(unique(masterDF[[c]]$EntryBin)+psthWinMax) #Select bins to be plotted
                subFRcols <- colnames(masterDF[[c]]) %in% as.character(subFRcolNames)
                masterDF[[c]]$correspBins <- correspBins

                sapply(seq(1, nDivisions, by=1), function(i){
                        
                        if(sum(correspBins==i, na.rm=T)>0){
                                
                                if(dataProcess=="Zscores"){
                                        plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, 3*(nDivisions+1)*yAxMaxZ))
                                        dataSel <- filter(masterDF[[c]], correspBins==i)
                                        
                                        if(cueExcOnly==T){
                                                dataSel <- filter(dataSel, CueExcited==T)
                                        }
                                        
                                        #If I'm going to be examining the period before the entry, exclude those trials in which the window includes the cue
                                        if(WdwStart<0){
                                                dataSel <- filter(dataSel, (CueLat*1000)>=abs(WdwStart))
                                        }
                                        
                                        numericDF <- apply(dataSel[, subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                        BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                        
                                        #Sometimes I get units with BLsd=0, and it gives me Inf values when calculating Zscores. Make it 1 so that it doesn't yield an INF
                                        BLsd[BLsd==0] <- 1
                                        numericDFZsc <- sapply(seq(1, nrow(numericDF)), function(w){
                                                sapply(seq(1, length(numericDF[w,])), function(v){
                                                        ZscoreCalc(x=numericDF[w,][v], avg=BLaverage[w], sd=BLsd[w])})
                                                })
                                        numericDFZsc <- t(numericDFZsc)
                                        
                                        MeanByBin <- colMeans(numericDFZsc, na.rm=T)
                                        SEMByBin <- colSds(numericDFZsc, na.rm=T)/sqrt(nrow(numericDFZsc))
                                        
                                        yAxMax <- yAxMaxZ
                                        
                                        # ### TEST SPECIFIC UNITS IN SPECIFIC BINS
                                        # neuronidx <- unique(dataSel$allUnitIdx)
                                        # plot.new()
                                        # plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, length(neuronidx)*40+30))
                                        # test <- sapply(seq(1, length(neuronidx)), function(j){
                                        #         numericDF <- apply(dataSel[dataSel$allUnitIdx==neuronidx[j], subFRcols], MARGIN=2, as.numeric)
                                        #         BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                        #         BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                        #         numericDFZsc <- sapply(seq(1, nrow(numericDF)), function(w){
                                        #                 sapply(seq(1, length(numericDF[w,])), function(v){
                                        #                         ZscoreCalc(x=numericDF[w,][v], avg=BLaverage[w], sd=BLsd[w])})
                                        #         })
                                        #         numericDFZsc <- t(numericDFZsc)
                                        #         
                                        #         MeanByBin <- colMeans(numericDFZsc, na.rm=T)
                                        #         SEMByBin <- colSds(numericDFZsc, na.rm=T)/sqrt(nrow(numericDFZsc))
                                        #         
                                        #         errCloud(col="black", x=1:length(MeanByBin), y=j*40+MeanByBin, err=SEMByBin)
                                        #         lines(x=1:length(MeanByBin), y=j*40+MeanByBin)
                                        # 
                                        # })
                                        # 
                                        # 
                                        # 
                                        
                                } else {
                                        
                                        plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, 3*nDivisions*yAxMaxRaw))
                                        dataSel <- filter(masterDF[[c]], correspBins==i)
                                        if(cueExcOnly==T){
                                                dataSel <- filter(dataSel, CueExcited==T)
                                        }
                                        
                                        #If I'm going to be examining the period before the entry, exclude those trials in which the window includes the cue
                                        if(WdwStart<0){
                                                dataSel <- filter(dataSel, (CueLat*1000)>=abs(WdwStart))
                                        }
                                        
                                        numericDF <- apply(dataSel[,subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        MeanByBin <- colMeans(numericDF, na.rm=T)
                                        SEMByBin <- colSds(numericDFna.rm=T)/sqrt(nrow(numericDF))
                                       
                                }
                                
                                #These lines are the ones that actually plot the graph
                                abline(h=(3*i*yAxMax), col=color, lty=3)
                                errCloud(col=color[c], x=1:ncol(numericDF), y=(3*i*yAxMax)+MeanByBin, err=SEMByBin)
                                lines(x=1:ncol(numericDF), y=(3*i*yAxMax)+MeanByBin, col=color[c], lwd=1.5, type='l')
                                #if(trialBinSize==1){text(paste("Trial", trialBins[i], sep=" "), col="red", x=5, y=(3*i*yAxMax)+yAxMax*1.5)}
                                # else {text(paste("Trials", trialBins[i], "to", trialBins[i]+trialBinSize-1, sep=" "), col="red", x=5, y=(3*i*yAxMax)+yAxMax*1.5)}
                        }
                
                })
                
        })
        
        #I have to do this again outside the function (for the axis)
        if(capped==T){
                nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
        } else {
                nDivisions <- round(((max(masterDF[[1]]$trialfromCP, na.rm=T)-min(masterDF[[1]]$trialfromCP, na.rm=T))/trialBinSize), 0)
                trialBins <- seq(min(masterDF[[1]]$trialfromCP, na.rm=T), max(masterDF[[1]]$trialfromCP, na.rm=T), by=trialBinSize)
        }
        psthWins <- ((WdwStart)/binw):((WdwEnd)/binw)
        psthWinMin <- (WdwStart)/binw
        psthWinMax <- (WdwEnd)/binw
        subFRcolNames <- (unique(masterDF[[c]]$EntryBin)+psthWinMin):(unique(masterDF[[c]]$EntryBin)+psthWinMax) #Select bins to be plotted
        
        if(dataProcess=="Zscores"){yAxMax <- yAxMaxZ; labelLeg <- paste(yAxMax, "(Zsc.)")}
        if(dataProcess=="raw"){yAxMax <- yAxMaxRaw; labelLeg <- paste(yAxMax, "(Hz)")}
       
        steps <- seq(1, length(subFRcolNames), by=10)
        labs <- seq(WdwStart/1000, WdwEnd/1000, length.out = length(steps))
        axis(side=1, at=steps, labels=labs)
        mtext(side=1, text=paste("Time from ",  eventName, " (s)", sep=""), cex=1.2, line=2.5)
        abline(v=seq(1, length(subFRcolNames))[seq(WdwStart, WdwEnd, by=binw)==0], lwd=1.5)
        
        atY <- (1:nDivisions)*3*yAxMax
        axis(side=2, at=atY, labels = trialBins[-length(trialBins)], las=2, cex=1.4)
        mtext(side=2, line=2.5, text="Trial from CP", font=2, cex=1.5)
        
        #abline(h=40, lty=2)
        segments(x0=length(subFRcolNames)-2, x1=length(subFRcolNames), y0=nDivisions*3*yAxMax, y1=nDivisions*3*yAxMax, lwd=1.5)
        segments(x0=length(subFRcolNames), x1=length(subFRcolNames), y0=nDivisions*3*yAxMax, y1=(nDivisions*3*yAxMax)+3*yAxMax, lwd=1.5)
        text(labelLeg, x=length(subFRcolNames)+1, y=(nDivisions*3*yAxMax)+1.5*yAxMax, srt=90, cex=0.8)  

        dev.off()
}




save(plotFRandCPhistogram_Entry, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/plotFRandCPhistogram_Entry.R")
save(plotFRandCPhistogram_Entry, file="E:/Dropbox/NMDA/R Functions/plotFRandCPhistogram_Entry.R")
