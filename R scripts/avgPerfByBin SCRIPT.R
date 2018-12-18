avgPerfByBin <- function(binsize=10, colors=c("gray20", "gray50"), YMIN=NA, YMAX=NA, 
                         data=list(DStaskAcc, NStaskAcc), cues=c("S+", "S-"), 
                         index="Performance index", legendLocation="topleft", y_axis_label, 
                         behGraphFolder=behGraphFolder, limit=240, plot=T){
        
        forPlotDF <- lapply(seq(1, length(data)), function(x){ #For each index I am interested in
                
                dataSel <- data[[x]]
        
                #If my data object is ITI latency, that object includes ITIs before NS and DS. I'll choose the DS only
                anyITI <- cues=="ITI"
                
                if(sum(anyITI)>0 & cues[x]=="ITI"){dataSel <- lapply(seq(1, length(dataSel)), function(k){
                        ratSess <- idx[[k]]
                        ITIcueKind <- do.call("c", lapply(seq(1, length(ratSess)), function(l){
                                alldata[[ratSess[l]]]$orderCues
                        }))
                        
                        #Select the ITI latency values of the DS trials only
                        dataSel[[k]] <- dataSel[[k]][ITIcueKind==1]
                        
                })}
                
                #In case I only want to plot a subsection of trials
                if(!is.na(limit)){dataSel <- lapply(seq(1, length(dataSel)),  function(m){dataSel[[m]][1:limit]})}
                
                
                byRatPerBin <- lapply(seq(1, length(dataSel)), function(y){ #For each rat
                        
                       
                        
                        trialIdx <- 1:length(dataSel[[y]])
                        bincuts <- seq(1, max(trialIdx, na.rm=T), by=binsize)
                        binIdx <- findInterval(trialIdx, bincuts)
                        sapply(seq(1, max(binIdx)), function(z){ #For each bin
                                mean(dataSel[[y]][binIdx==z])
                        })
                })
                
                #Some rats have less bins than others (because they received less trials (if they lost the cap, or different sessino protocols w different #trials...). Fill in with NAs the bins that they don't have
                maxBin <- max(sapply(byRatPerBin, length), na.rm=T)
                
                forPlot <- as.data.frame(sapply(seq(1, length(byRatPerBin)), function(y){
                        missingBins <- maxBin-length(byRatPerBin[[y]])
                        c(byRatPerBin[[y]], rep(NA, missingBins))
                })
                )
                
                meanPerBin <- rowMeans(forPlot, na.rm=T)
                sdPerBin <- sapply(seq(1, nrow(forPlot)), function(z){sd(forPlot[z,], na.rm=T)})
                SEMPerBin <- sdPerBin/sqrt(ncol(forPlot))
                
                forPlotDF <- cbind(meanPerBin, SEMPerBin)
        })
        
        forPlotDF <- lapply(forPlotDF, as.data.frame) #I have to do this to get forPlotDF to have the properties of a data frame
        
        if(plot==T){
                
                filename=paste(behGraphFolder, paste(cues, index), " by ", binsize, " trial bin.pdf", sep="")
                
                pdf(file = filename)
                
                if(is.na(YMIN)){ymin=floor(min(unlist(forPlotDF), na.rm = T))} else {ymin=YMIN}
                if(is.na(YMAX)){ymax=ceiling(max(unlist(forPlotDF), na.rm = T))} else {ymax=YMAX}
                
                yrange <- ymax-ymin
                
                nBins=max(sapply(forPlotDF, nrow))
                
                plot.new()
                plot.window(xlim=c(0, nBins), ylim=c(ymin, ymax))
                
                if(ymin<0){abline(h=0, lty=3)} #Mark 0 if the y axis spans negative to postive numbers
                
                lapply(seq(1, length(forPlotDF)), function(x){
                        lines(x=seq(1, nBins), y=forPlotDF[[x]][,1], lwd=2, col=colors[x])
                        if(x>1){j=x*0.005} else {j=0}
                        errBars(x=seq(1, nBins), y=forPlotDF[[x]][,1], err=forPlotDF[[x]][,2], jitter=j, color = colors[x])
                })
                
                xaxisLab <- binsize * (1:nBins)
                axis(side=1, at=seq(1, nBins), labels = xaxisLab, cex.axis=1.4)
                
                if(yrange<=1){atVals=seq(0, 1, 0.2)} else {atVals <- seq(ymin, ymax)} #In response ratio graphs, R automatically marks only 0 and 1 in the y axis. I want more breaks
                axis(side=2, las=2, at=atVals, cex.axis=1.4)
                mtext(side=1, line=3, text="Trial", cex = 1.5, font=2)
                mtext(side=2, line=3, text=y_axis_label, cex=1.5, font=2)
                
                
                legend(x=legendLocation, legend=cues, lty=1, lwd=2, col=colors, cex=1.5, bty="n")
                
                legend("bottomleft", legend=paste("Bin size = ", binsize), cex=1.2, bty="n")
                
                dev.off()
        }
        
        print(forPlotDF)
}


save(avgPerfByBin, file="E:/Dropbox/NMDA/R functions/avgPerfByBin.R")
save(avgPerfByBin, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R functions/avgPerfByBin.R")
save(avgPerfByBin, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R functions/avgPerfByBin.R")

