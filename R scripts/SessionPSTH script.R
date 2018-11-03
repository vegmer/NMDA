

SessionPSTH <- function(experiment="Exp 4 Unil AP5 EXT", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
                        graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
                        yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, psthmin=-1, psthmax=1, imgFormat="pdf", 
                        neudata=allNeuronsDS){
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
                
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                
                #Subset the columns with the bins around the event that I'm interested in (defined by psthmin and psthmax)
                subFRcolNames <- (unique(masterDF[[c]]$CueBin)-psthWinMin):(unique(masterDF[[c]]$CueBin)+psthWinMax) #Select bins to be plotted
                subFRcols <- colnames(masterDF[[c]]) %in% as.character(subFRcolNames)
                
                #Subset the trials (correct only or all) and the units (cue excited or not) to be used for the PSTH
                dataSel <- masterDF[[c]]
                if(correctOnly==TRUE){
                        dataSel <- dataSel[!is.na(dataSel$CSplusResponse), ] 
                }
                if(cueExcOnly==TRUE){
                        dataSel <- dataSel[dataSel$CueExcited==TRUE, ]
                }
                
                ZscoreCalc <- function(x, avg, sd){(x-avg)/sd}
                
                #Calculate the mean average rate per bin around the time of the event of interest (in Z scores or raw scores)
                if(dataProcess=="Zscores"){
                        plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(yAxMinZ, yAxMaxZ))
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
                        plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, yAxMaxRaw))
                        numericDF <- apply(dataSel[,subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                        MeanByBin <- colMeans(numericDF, na.rm=T)
                        SEMByBin <- colSds(numericDFna.rm=T)/sqrt(nrow(numericDF))
                        yAxMax=yAxMaxRaw
                        labelLeg=paste(yAxMax, "(Hz)")
                }
                                        
                abline(h=0, col="gray30", lty=3)
                
                #Plot the actual PSTH
                errCloud(col=color[c], x=1:ncol(numericDF), y=MeanByBin, err=SEMByBin)
                lines(x=1:ncol(numericDF), y=MeanByBin, col=color[c], lwd=1.5, type='l')
                                
                })
                
        labs <- seq(psthmin, psthmax, by=binw/100) #Labels of psthmin to psthmax by 0.5s (50ms/100)
        steps <- seq(1, length(subFRcolNames), length.out = length(labs)) #Location of the labels
        axis(side=1, at=steps, labels=labs)
        axis(side=2, las=2)
        mtext(side=1, "Time from S+ onset (s)", cex=1.2, line=2.5, font=2)
        mtext(side=2, "Firing rate", line=2.5, font=2, cex=1.2)
        abline(v=steps[labs==0], lwd=1.5)
        #abline(v=seq(1, length(subFRcolNames))[seq(-psthmin, psthmax, by=binw/1000)==0], lwd=1.5)
        #abline(h=40, lty=2)
                        
          
        
        dev.off()
}


save(SessionPSTH, file="E:/Dropbox/NMDA/R Functions/SessionPSTH.R")
