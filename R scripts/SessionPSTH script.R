

SessionPSTH <- function(experiment="Exp 4 Unil AP5 EXT", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
                        graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
                        yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, psthmin=-1, psthmax=1, imgFormat="pdf", 
                        capped=F, capValue=c(1, 15),
                        neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=FALSE, legendLabels=c("S+ VEH side L", "S+ AP5 side L")){
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        if(!require(matrixStats)){install.packages("matrixStats")}
        library(matrixStats)
        
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        if(cueExcOnly==TRUE){unitSel="CueExcOnly"} else {unitSel="all units"}
        
        if(imgFormat=="pdf"){pdf(file=paste(graphFolder, experiment, "PSTH", dataProcess, trialSel, unitSel, ".pdf", sep="_"))}
        if(imgFormat=="png"){pnf(filename=paste(graphFolder, experiment, "PSTH", dataProcess, trialSel, unitSel, ".pdf", sep="_"))}
        
        binw <- neudata[[1]]$parameters$binw
        if(is.null(binw)){binw <- neudata$parameters$binw} #Preventive troubleshooting: In case I define neudata not as a list (as it should), but as a single object
        
        psthWins <- ((psthmin*1000)/binw):((psthmax*1000)/binw)
        psthWinMin <- (psthmin*1000)/binw
        psthWinMax <- (psthmax*1000)/binw
        
        plot.new()
        
        #Subset the columns with the bins around the event that I'm interested in (defined by psthmin and psthmax)
        FRcols <- (1:ncol(masterDF[[1]]))[is.element(colnames(masterDF[[1]]), 1:ncol(masterDF[[1]]))]
        subFRcolNames <- (unique(masterDF[[1]]$CueBin)-psthWinMin):(unique(masterDF[[1]]$CueBin)+psthWinMax) #Select bins to be plotted
        subFRcols <- colnames(masterDF[[1]]) %in% as.character(subFRcolNames)
        
        
        byBinbyUnit <- do.call("list", lapply(seq(1, length(masterDF)), function(c){
                
                dataSel <- masterDF[[c]]
                
                #Cue excited units taking into account the whole session
                cueexidx <- neudata[[c]]$cueexidx
                propcueex <- sum(unlist(cueexidx))/length(unlist(cueexidx))
                
                if(capped==TRUE){
                        dataSel <- dataSel[dataSel$trialIdx>=capValue[1] & dataSel$trialIdx<=capValue[2], ]
                        
                        masterlist <- neudata[[c]]$masterlist
                        
                        #Cue exc units taking into account those specific trials only
                        cueexidx = lapply(seq(1, length(masterlist)), function(x) { #For each rat
                                lapply(seq(1, length(masterlist[[x]])), function(y) { #For each unit
                                        
                                        trialSel <- capValue[1]:capValue[2]
                                        pbin = .05
                                        threshold = 3
                                        BLwdw=2
                                        
                                        allvals = unlist(masterlist[[x]][[y]][trialSel])
                                        hcounts = hist(allvals[which(allvals >= -BLwdw & allvals <= .5)], breaks = seq(-BLwdw, .5, pbin), plot = F)$counts 
                                        
                                        freq = (hcounts/pbin)/length(masterlist[[x]][[y]])  #just for comparison's sake, to see how counts compare to actual frequency
                                        baseline = mean(freq[1:(BLwdw/pbin)])   
                                        
                                        baselinecounts=hcounts[1:(BLwdw/pbin)]
                                        critwin = hcounts[(BLwdw/pbin + 1):((BLwdw+.5)/pbin)]         #this is the window (500 ms) in which we look for excitation
                                        
                                        critval = poisson.test(x = sum(baselinecounts), T = length(baselinecounts), conf.level = 0.999)$conf.int[2]   #computes the upper limit of the confindence interval based on the baseline
                                        
                                        diffs = diff(which(critwin > critval))       #computes the differences in indices for bins exceeding the critical value (to check whether excited bins are consecutive or not)
                                        cueex = F
                                        if(length(which(rle(diffs)$values == 1 & rle(diffs)$lengths >= threshold))>0) (cueex = T)   #looks for consecutive bins (diff equal to 1) that are at least 2 bins long -or whatever #of bins specified by the argument "threshold" (rle length of at least 2)
                                        
                                        return(cueex)  
                                })
                        })
                        
                        
                        # propcueex <- lapply(seq(1, length(cueexidx)), function(x){
                        #         length(which(cueexidx[[x]]==T))/length(cueexidx[[x]])
                        # })
                        
                        propcueex <- sum(unlist(cueexidx))/length(unlist(cueexidx))
                }
                
                #Subset the trials (correct only or all) and the units (cue excited or not) to be used for the PSTH
                if(correctOnly==TRUE){
                        dataSel <- dataSel[!is.na(dataSel$CueResponse), ] 
                }
                
                
                if(cueExcOnly==TRUE){
                        
                        cueExc <- unlist(cueexidx)[dataSel$allUnitIdx]
                        dataSel$CueExcited <- cueExc
                        
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
                        
                        
                } else {
                        plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, yAxMaxRaw))
                        numericDF <- apply(dataSel[,subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                        MeanByBin <- colMeans(numericDF, na.rm=T)
                        SEMByBin <- colSds(numericDFna.rm=T)/sqrt(nrow(numericDF))
                        yAxMax=yAxMaxRaw
                        labelLeg=paste(yAxMax, "(Hz)")
                }
                                        

                #Plot the actual PSTH
                errCloud(col=color[c], x=1:ncol(numericDF), y=MeanByBin, err=SEMByBin)
                lines(x=1:ncol(numericDF), y=MeanByBin, col=color[c], lwd=1.5, type='l')
                
                results <- list(numericDF)
                names(results) <- legendLabels[c]
                
                nUnits <- length(unique(dataSel$allUnitIdx))
                text(x=10, y=5-c, labels=paste("n=", nUnits, "; ", "Prop Cue Exc=", round(propcueex, 2)), col=colindx[c])
                
                return(results)              
                })
        )
        
        abline(h=0, col="gray30", lty=3)
        labs <- seq(-psthmin, psthmax, by=binw/100) #Labels of psthmin to psthmax by 0.5s (50ms/100)
        steps <- seq(1, length(subFRcolNames), length.out = length(labs)) #Location of the labels
        axis(side=1, at=steps, labels=labs)
        axis(side=2, las=2)
        mtext(side=1, "Time from S+ onset (s)", cex=1.2, line=2.5, font=2)
        mtext(side=2, "Firing rate", line=2.5, font=2, cex=1.2)
        abline(v=steps[labs==0], lwd=1.5)
        #abline(v=seq(1, length(subFRcolNames))[seq(-psthmin, psthmax, by=binw/1000)==0], lwd=1.5)
        #abline(h=40, lty=2)
        legend("topleft", legend=legendLabels, col=color, lty=1)   
        
        dev.off()
        
        byBinbyUnit
}


save(SessionPSTH, file="E:/Dropbox/NMDA/R Functions/SessionPSTH.R")
