

plotFRandCPhistogram <- function(experiment="Exp 3", events=c("S+", "S-"), masterDF=list(masterDF_DS, masterDF_NS), 
                                 graphFolder=MixedGraphFolder, trialBinSize=15, dataProcess="Zscores", 
                                 correctOnly=FALSE, color="black", capped=T, capValue = c(-60, 60), 
                                 yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, 
                                 imgFormat="pdf", neudata=allNeuronsDS, cueExcOnly=FALSE){
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        if(!require(matrixStats)){install.packages("matrixStats")}
        library(matrixStats)
        
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        
        if(imgFormat=="pdf"){pdf(file=paste(graphFolder, experiment, "FR by trial from CP PSTH", "bin size", trialBinSize, dataProcess, trialSel, ".pdf", sep="_"))}
        if(imgFormat=="png"){png(filename=paste(graphFolder, experiment, "FR by trial from CP PSTH", "bin size", trialBinSize, dataProcess, trialSel, ".png", sep="_"))}
        
        binw <- neudata$parameters$binw
        psthWins <- ((psthmin*1000)/binw):((psthmax*1000)/binw)
        psthWinMin <- (psthmin*1000)/binw
        psthWinMax <- (psthmax*1000)/binw
        
        plot.new()
        
        lapply(seq(1, length(masterDF)), function(c){
                
                if(capped==T){
                        nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                        trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
                } else {
                        nDivisions <- round(((max(masterDF[[c]]$trialfromCP, na.rm=T)-min(masterDF[[c]]$trialfromCP, na.rm=T))/trialBinSize), 0)
                        trialBins <- seq(min(masterDF[[c]]$trialfromCP, na.rm=T), max(masterDF[[c]]$trialfromCP, na.rm=T), by=trialBinSize)
                }
                
                correspBins <- findInterval(masterDF[[c]]$trialfromCP, trialBins)
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                
                subFRcolNames <- (unique(masterDF[[c]]$CueBin)-psthWinMin):(unique(masterDF[[c]]$CueBin)+psthWinMax) #Select bins to be plotted
                subFRcols <- colnames(masterDF[[c]]) %in% as.character(subFRcolNames)
                masterDF[[c]]$correspBins <- correspBins
                
                ZscoreCalc <- function(x, avg, sd){(x-avg)/sd}
                
                ### ALL TRIALS (RESPONDED TO AND MISSED)
                if(correctOnly==FALSE){
                        for(i in 1:nDivisions){
                                if(sum(correspBins==i, na.rm=T)>0){
                                        if(dataProcess=="Zscores"){
                                                plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, 3*(nDivisions+1)*yAxMaxZ))
                                                dataSel <- filter(masterDF[[c]], correspBins==i)
                                                
                                                if(cueExcOnly==TRUE){
                                                        #If I select Cue-exc only, in truth I'm interested in S+ excited units (not S- exc units). So if the event used to make the masterDF was a S-, I won't have the selection of S+ units (but of S- units). Fix that
                                                        if(events[c]=="S-"){
                                                                excIndexAll <- unlist(neudata$cueexidx)
                                                                excIndexSel <- excIndexAll[dataSel$allUnitIdx]
                                                                dataSel$CueExcited <- excIndexSel
                                                        }
                                                        
                                                        dataSel <- dataSel[dataSel$CueExcited==TRUE, ]
                                                }
                                                
                                                numericDF <- apply(dataSel[, subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                                BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                                BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                                numericDFZsc <- sapply(seq(1, nrow(numericDF)), function(w){
                                                        sapply(seq(1, length(numericDF[w,])), function(v){
                                                                ZscoreCalc(x=numericDF[w,][v], avg=BLaverage[w], sd=BLsd[w])})
                                                })
                                                numericDFZsc <- t(numericDFZsc)
                                                
                                                MeanByBin <- colMeans(numericDFZsc, na.rm=T)
                                                SEMByBin <- colSds(numericDFZsc, na.rm=T)/sqrt(nrow(numericDFZsc))
                                                yAxMax=yAxMaxZ
                                                labelLeg=paste(yAxMax, "(Zsc.)")
                                                
                                                
                                        } else {
                                                plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, 3*nDivisions*yAxMaxRaw))
                                                dataSel <- filter(masterDF[[c]], correspBins==i)
                                                if(cueExcOnly==TRUE){
                                                        dataSel <- dataSel[dataSel$CueExcited==TRUE, ]
                                                }
                                                numericDF <- apply(dataSel[,subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                                MeanByBin <- colMeans(numericDF, na.rm=T)
                                                SEMByBin <- colSds(numericDFna.rm=T)/sqrt(nrow(numericDF))
                                                yAxMax=yAxMaxRaw
                                                labelLeg=paste(yAxMax, "(Hz)")
                                        }
                                        
                                        abline(h=(3*i*yAxMax)+0, col="gray60", lty=3)
                                        errCloud(col=color[c], x=1:ncol(numericDF), y=(3*i*yAxMax)+MeanByBin, err=SEMByBin)
                                        lines(x=1:ncol(numericDF), y=(3*i*yAxMax)+MeanByBin, col=color[c], lwd=1.5, type='l')
                                        
                                        if(trialBinSize==1){text(paste("Trial", trialBins[i], sep=" "), col="red", x=5, y=(3*i*yAxMax)+yAxMax*1.5)}
                                        else {text(paste("Trials", trialBins[i], "to", trialBins[i]+trialBinSize-1, sep=" "), col="red", x=5, y=(3*i*yAxMax)+yAxMax*1.5)}
                                }
                                
                        }
                        
                        labs <- seq(-psthmin, psthmax, by=binw/100)
                        steps <- seq(1, length(subFRcolNames), by=10)
                        axis(side=1, at=steps, labels=labs)
                        mtext(side=1, "Time from S+ onset (s)", cex=1.2, line=2.5)
                        abline(v=seq(1, length(subFRcolNames))[seq(-psthmin, psthmax, by=binw/1000)==0], lwd=1.5)
                        #abline(h=40, lty=2)
                        segments(x0=length(subFRcolNames)-2, x1=length(subFRcolNames), y0=nDivisions*3*yAxMax, y1=nDivisions*3*yAxMax, lwd=1.5)
                        segments(x0=length(subFRcolNames), x1=length(subFRcolNames), y0=nDivisions*3*yAxMax, y1=(nDivisions*3*yAxMax)+3*yAxMax, lwd=1.5)
                        text(labelLeg, x=length(subFRcolNames)+1, y=(nDivisions*3*yAxMax)+1.5*yAxMax, srt=90, cex=0.8)  
                        
                        atY <- (0.5:(nDivisions+0.5))*3*yAxMax
                        if(length(trialBins)==1){trialBins <- c(trialBins, trialBinSize)}
                        axis(side=2, at=atY, labels = trialBins, las=2, cex=1.4)
                        
                        mtext(side=2, line=2.5, text="Trial from CP", font=2, cex=1.5)
                        
                }
                
                # If I select only correct trials
                if(correctOnly==TRUE){
                        for(i in 1:nDivisions){
                                #dataSelIdx <- !is.na((1:nrow(masterDF[[c]]))[correspBins==i])
                                #TrialsInBinDF <- masterDF[[c]][dataSelIdx,] #DF with rows that correspond to bin i
                                TrialsInBinDF <- filter(masterDF[[c]], correspBins==i)
                                
                                if(cueExcOnly==TRUE){
                                        TrialsInBinDF <- TrialsInBinDF[TrialsInBinDF$CueExcited==TRUE, ]
                                }
                                
                                if(sum(!is.na(TrialsInBinDF$CSplusResponse))>0){
                                        DSrespDF <- TrialsInBinDF[!is.na(TrialsInBinDF$CSplusResponse),] #Index of rows of correct DS trials
                                        numericDF <- apply(DSrespDF[, subFRcols], MARGIN=2, as.numeric)
                                        
                                        if(dataProcess=="Zscores"){
                                                plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, 3*(nDivisions+1)*yAxMaxZ))
                                                BLaverage <- as.numeric(format(DSrespDF$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                                BLsd <- as.numeric(format(DSrespDF$BLsd, digits=2))
                                                numericDFZsc <- sapply(seq(1, nrow(numericDF)), function(w){
                                                        sapply(seq(1, length(numericDF[w,])), function(v){
                                                                ZscoreCalc(x=numericDF[w,][v], avg=BLaverage[w], sd=BLsd[w])})
                                                })
                                                numericDFZsc <- t(numericDFZsc)
                                                
                                                MeanByBin <- colMeans(numericDFZsc)
                                                SEMByBin <- colSds(numericDFZsc)/sqrt(nrow(numericDFZsc))
                                                yAxMax=yAxMaxZ
                                                labelLeg=paste(yAxMaxZ, "(Zsc.)")
                                                
                                        } else {
                                                plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, 3*nDivisions*yAxMaxRaw))
                                                MeanByBin <- colMeans(numericDF)
                                                SEMByBin <- colSds(numericDF)/sqrt(nrow(numericDF))
                                                yAxMax=yAxMaxRaw
                                                labelLeg=paste(yAxMax, "(Hz)")}
                                        
                                        errCloud(col=color, x=1:ncol(numericDF), y=(3*i*yAxMax)+MeanByBin, err=SEMByBin)
                                        lines(x=1:ncol(numericDF), y=(3*i*yAxMax)+MeanByBin, col=color, lwd=1.5, type='l')
                                        #if(trialBinSize==1){text(paste("Trial", trialBins[i], sep=" "), col="red", x=10, y=(3*i*yAxMax)+yAxMax*1.5)}
                                        #else {text(paste("Trials", trialBins[i], "to", trialBins[i]+trialBinSize-1, sep=" "), col="red", x=10, y=(3*i*yAxMax)+yAxMax*1.5)}
                                        
                                        
                                }
                                mtext(side=1, "Time from S+ onset (s)", cex=1.2, line=2.5)
                                abline(v=seq(1, length(subFRcolNames))[seq(psthmin, psthmax, by=binw/1000)==0])
                                segments(x0=length(subFRcolNames)-2, x1=length(subFRcolNames), y0=nDivisions*3*yAxMax, y1=nDivisions*3*yAxMax)
                                segments(x0=length(subFRcolNames), x1=length(subFRcolNames), y0=nDivisions*3*yAxMax, y1=(nDivisions*3*yAxMax)+3*yAxMax)
                                text(labelLeg, x=length(subFRcolNames)+1, y=(nDivisions*3*yAxMax)+1.5*yAxMax, srt=90, cex=0.8)
                                
                                atY <- (1:nDivisions)*3*yAxMax
                                axis(side=2, at=atY, labels = trialBins[-length(trialBins)], las=2, cex=1.4)
                                mtext(side=2, line=2.5, text="Trial from CP", font=2, cex=1.5)
                                
                                #abline(h=40, lty=2)
                                segments(x0=length(subFRcolNames)-2, x1=length(subFRcolNames), y0=nDivisions*3*yAxMax, y1=nDivisions*3*yAxMax, lwd=1.5)
                                segments(x0=length(subFRcolNames), x1=length(subFRcolNames), y0=nDivisions*3*yAxMax, y1=(nDivisions*3*yAxMax)+3*yAxMax, lwd=1.5)
                                text(labelLeg, x=length(subFRcolNames)+1, y=(nDivisions*3*yAxMax)+1.5*yAxMax, srt=90, cex=0.8)  
                                
                        }
                        
                }
                
                
        })
        
        
        dev.off()
}




save(plotFRandCPhistogram, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/plotFRandCPhistogram.R")
save(plotFRandCPhistogram, file="E:/Dropbox/NMDA/R Functions/plotFRandCPhistogram.R")
save(plotFRandCPhistogram, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R Functions/plotFRandCPhistogram.R")
