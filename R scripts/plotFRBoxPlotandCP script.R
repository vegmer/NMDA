##### PLOT FIRING RATE POST CUE AS A FUNCTION OF DISTANCE FROM CHANGE POINT

plotFRBoxPlotandCP <- function(cue="S+", experiment, masterDF, graphFolder=MixedGraphFolder, trialBinSize=40, 
                        WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, points=TRUE, lines=F,
                        color="black", legLabels=c("VEH side", "AP5 side"), yAxMinZ=-1, yAxMaxZ=1, yAxMaxRaw=7, 
                        capped=T, capValue=c(-120, 120), cueExcOnly=F, neudata=allNeuronsDS) {
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        
        
        binw <- neudata$parameters$binw
        
        #IF I WANT TO PLOT MORE THAN ONE GROUP PER GRAPH, MAKE A LIST WITH THE DIFFERNT masterDF OF THE DIFFERENT GROUPS
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        
        filename=paste(graphFolder, experiment, "BOXPLOT FR trial from CP", "bin size", trialBinSize, WdwStart, WdwEnd, dataProcess, trialSel, ".pdf", sep="_")
        
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
                                        dataSel <- filter(dataSel, (CueLat*1000)>=WdwEnd)
                                }
                                
                                if(correctOnly==T){
                                        dataSel <- filter(dataSel, !is.na(dataSel$CueResponse))
                                }
                                
                                #It's possible that after subsetting by all those criteria (cue excited, etc.), we're left with no rows. In that case, skip the rest of the codde
                                if(nrow(dataSel)>0){
                                        screenRange <- c(0, nDivisions+1)
                                        
                                        plot.window(xlim=screenRange, ylim=c(yAxMinZ, yAxMaxZ))
                                        
                                        BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                        BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                        
                                        numericDF <- apply(dataSel[, FRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        
                                        if(nrow(numericDF)>1){numericDFTargetWdw <- numericDF[, TargetStartBin:TargetEndBin]} else {
                                                numericDFTargetWdw <- numericDF[TargetStartBin:TargetEndBin]
                                        }
                                        
                                        
                                        binUnits <- unique(dataSel$allUnitIdx)
                                        
                                        byBinbyUnit <- do.call("rbind", lapply(seq(1, length(binUnits)), function(u){
                                                
                                                unitSel <- dataSel$allUnitIdx==binUnits[u]
                                                byUnitbyTrialFR <- numericDFTargetWdw[unitSel, ]
                                                
                                                byUnitFR <-  mean(byUnitbyTrialFR)
                                                
                                                if(dataProcess=="Zscores"){
                                                        BLaverage <- BLaverage[unitSel][1]
                                                        BLsd <- BLsd[unitSel][1]
                                                        byUnitFR <- (byUnitFR-BLaverage)/BLsd  
                                                }
                                                
                                                sessInfo <- dataSel[unitSel, c(2:8, 28)][1, ]
                                                
                                                #Summary of performance index for the trials recorded on this unit on this session
                                                perfData <- dataSel[unitSel, c(17, 16, 15, 19) ]
                                                perfData$CueResponse <- !is.na(perfData$CueResponse)
                                                
                                                perfData$CueSpecif <- as.numeric(as.character(perfData$CueSpecif))
                                                perfData$CueLat <- as.numeric(as.character(perfData$CueLat))
                                                perfData$CueLat[which(perfData$CueLat==10)] <- NA
                                                perfData$CueResponse <- as.logical(perfData$CueResponse)
                                                perfData$ITIlatency <- as.numeric(as.character(perfData$ITIlatency))
                                                
                                                
                                                perfSumm <- sapply(1:length(perfData), function(j){mean(perfData[, j], na.rm=T)})
                                                names(perfSumm) <- colnames(perfData)
                                                perfSumm <- as.list(perfSumm)
                                                
                                                return(data.frame(bin=trialBins[i], unitIdx=binUnits[u], sessInfo, perfSumm, byUnitFR, cue=cue[c]))
                                                
                                        })
                                        )
                                        
                                        #Plot the actual values of FR in selected bins of trials with respect to CP
                                        forBoxplot <- summary(byBinbyUnit$byUnitFR)
                                        
                                        Q1 <- forBoxplot[2]; Q3 <- forBoxplot[5]
                                        MEDIAN <- forBoxplot[3]; AVG <- forBoxplot[4]
                                        MIN <- forBoxplot[1]; MAX <- forBoxplot[6]
                                        
                                        if(length(masterDF)==1){xleft=-0.3; xright=0.3}
                                        
                                        if(length(masterDF) >1){
                                                minBoxArea <- 0.3
                                                maxBoxArea <- 0.3
                                                SubBoxArea <- (maxBoxArea+minBoxArea)/length(masterDF)
                                                xleftSeq <- seq(-minBoxArea, maxBoxArea-SubBoxArea, length.out = length(masterDF))
                                                xrightSeq <- seq(-minBoxArea+SubBoxArea, maxBoxArea, length.out = length(masterDF))
                                                
                                                xleft=xleftSeq[c]
                                                xright=xrightSeq[c]
                                                
                                                legend(x=length(trialBins)+0.5, y=yAxMaxZ-c/2, legend=legLabels[c], 
                                                       fill=color[c], pch=19, border='white', bty = 'n')
                                        }
                                        
                                        center <- i #Because I want my boxplots to be between the bin cuts (and plot.window said X goes from 0 to 6 and I have 6 bins)
                                        
                                        rect(xleft=center+xleft, xright=center+xright, ybottom=Q1, ytop=Q3, col=color[c], border = "white")
                                        segments(x0=center+xleft, x1=center+xright, y0=MEDIAN, y1=MEDIAN, col="white", lwd=2)
                                        segments(x0=center+xleft, x1=center+xright, y0=AVG, y1=AVG, col="black", lwd=2)
                                        
                                        if(points==TRUE){
                                                segments(x0=center+xleft, x1=center+xright, y0=MAX, y1=MAX, col="black", lwd=1.5)
                                                segments(x0=center+xleft, x1=center+xright, y0=MIN, y1=MIN, col="black", lwd=1.5)
                                                
                                                halfBox <- SubBoxArea/2
                                                points(x=rep(center+xleft+halfBox, nrow(byBinbyUnit)), y=byBinbyUnit$byUnitFR, pch=19, cex=0.2, col="black")
                                                
                                        }
                                        
                                        
                                        
                                        byBinbyUnit
                                }
                               
                        }
                })
                )
                
                
                
                #Axis
                if(capped==T){
                        trialfromCPvalues <- capValue[1]:capValue[2]
                } else {
                        trialfromCPvalues <- min(masterDF[[c]]$trialfromCP, na.rm=TRUE):max(masterDF[[c]]$trialfromCP, na.rm=TRUE)   
                }
                
                
                windowLandmarks <- c(0, 1:nDivisions, nDivisions+1) #Based on the X parameter of plot.window()
                atVals <- (windowLandmarks-0.5)[-1] #Location of the plots
                
                axis(side=1, at=atVals[1:length(trialBins)], labels=trialBins, cex.axis=1.8)
                
                abline(v=atVals[which(trialBins==0)])
                axis(side=2, las=2, cex.axis=1.8)
                
                mtext("Trial from behavioral change point", side=1, cex=1.2, line=2.5)
                mtext(paste("Bin size=", trialBinSize, "trials", sep=" "), side=3, adj=0.05, line=-2, font=3)
                
                return(datperbin)
                
        })
        )
        
        
        if(length(masterDF) >1 & lines==TRUE){ #It makes lines between the same neurons on the S+ trials vs S- trials. The DFs have to be consecutive as indicated by the "cue" parameter
                
                DS_DF <- (1:length(cue))[cue %in% "S+"]
                NS_DF <- (1:length(cue))[cue %in% "S-"]
                
                minBoxArea <- 0.3
                maxBoxArea <- 0.3
                SubBoxArea <- (maxBoxArea+minBoxArea)/length(masterDF)
                halfBoxArea <- SubBoxArea/2
                xleftSeq <- seq(-minBoxArea, maxBoxArea-SubBoxArea, length.out = length(masterDF))
                
                
                sapply(seq(1, length(masterDF)/2), function(c){
                        
                      DS_DFidx <- DS_DF[c] 
                      NS_DFidx <- NS_DF[c]
                      
                      DStrials <- allDat[[DS_DFidx]]
                      NStrials <- allDat[[NS_DFidx]]
                      
                      xleft=xleftSeq[DS_DFidx]+halfBoxArea
                      xright=xleftSeq[NS_DFidx]+halfBoxArea
                      
                      bins <- unique(DStrials$bin)
                      
                      lapply(seq(1, length(bins)), function(i){ #For each bin
                              
                              DStrialsPerBin <- DStrials[DStrials$bin==bins[i], ]
                              NStrialsPerBin <- NStrials[NStrials$bin==bins[i], ]
                              
                              sapply(seq(1, length(DStrialsPerBin$byUnitFR)), function(l){ #For each unit
                                      
                                      if(DStrialsPerBin$modality[l]=="98"){lineCol <- "black"} else {lineCol <- "gray30"}
                                      
                                      lines(x=c(i+xleft, i+xright), y=c(DStrialsPerBin$byUnitFR[l], NStrialsPerBin$byUnitFR[l]), col=lineCol)
                                      
                              })
                              
                              
                      })
                       
                })
                
                legend("topright", col=c("black", "gray"), legend=c("S+ tone", "S+ light"), lty=1, lwd=2)
                
        }
        
        #Title
        titleLabel <- paste("Post ", cue, " (", WdwStart, "-", WdwEnd,  "ms) firing rate with respect to behavioral change", sep="")
        if(correctOnly==T){subLabel <- "Correct trials only"} else {subLabel="All S+ trials"}
        mtext(titleLabel, side=3, line=2.5, cex=1.2, font=2)
        mtext(subLabel, side=3, line=1, cex=1, font=3)
        
        abline(h=0, lty=2)
        
        dev.off()
        
        
        return(allDat)
        
        
}



save(plotFRBoxPlotandCP, file="E:/Dropbox/DISSERTATION/R functions/plotFRBoxPlotandCP.r")
save(plotFRBoxPlotandCP, file="E:/Dropbox/NMDA/EXP1_Performance/R Functions/plotFRBoxPlotandCP.R")
save(plotFRBoxPlotandCP, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/plotFRBoxPlotandCP.R")
save(plotFRBoxPlotandCP, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R Functions/plotFRBoxPlotandCP.R")
save(plotFRBoxPlotandCP, file="E:/Dropbox/NMDA/R Functions/plotFRBoxPlotandCP.R")



