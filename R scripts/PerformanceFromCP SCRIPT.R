PerformanceFromCP <- function(relTrialMin=-100, relTrialMax=100, trialBinSize=5, typegraph="both", 
                              numSess=6, numCues=35, imgFormat="pdf", dataForRCumulative=dataForRCumulative, 
                              dataForRdir=dataForRdir, graphFolder=PerfRelToCPFolder, CPdata=CPdata, 
                              csacqidx=csacqidx, idx=idx, rats=rats, alldata=alldata, color="black"){
        
        #Names for graph titles. Make sure they are in the same order that these objects are loaded from the folder in the next line
        longNames <- c("S+ latency", "S+ resp ratio", "S+ specificity", "S+ time to spare", "ITI latency", "ITI response ratio", "S- latency", "S- resp ratio", "S- specificity", "S- time to spare")
        
        # Load relevant objects and create object 'cumData' with the labels for each object
        filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); cumData <- sapply(seq(1, length(filesCum)), function(i){load(filesCum[[i]])})
        
        # Plot
        lapply(seq(1, length(filesCum)), function(j){
                
                data <- get(cumData[j])
                
                totalCues <- max(sapply(data, length))
                
                
                #Make a list of vectors (one per animal) with trial number relative to change point
                
                #If my data object is ITI latency or ITI response ratio, that object includes ITIs before NS and DS. I'll choose the DS only
                if(j==5 | j==6){data <- lapply(seq(1, length(data)), function(k){
                        ratSess <- idx[[k]]
                        ITIcueKind <- do.call("c", lapply(seq(1, length(ratSess)), function(l){
                                alldata[[ratSess[l]]]$orderCues
                        }))
                        
                        #Select the ITI latency values of the DS trials only
                        data[[k]] <- data[[k]][ITIcueKind==1]
                        
                })}
                
                trialNumbers <- lapply(seq(1, length(data)), function(k){
                        
                        1:length(data[[k]])
                        
                })
                
                
                relTrialNumbers <- lapply(seq(1, length(data)), function(k){
                        trialNumbers[[k]]-CPdata$CP[k]
                })
                
                #Window of trials around CP
                trialVals <- relTrialMin:relTrialMax
                
                #if(imgFormat=="png"){filetitle=paste(graphFolder, paste(longNames[j]), "relative to CP",".png", sep=""); png(filename = filetitle)}
                #if(imgFormat=="pdf"){filetitle=paste(graphFolder, paste(longNames[j]), "relative to CP" ,".pdf", sep=""); pdf(file = filetitle)}
                
                perfRelCPInd <- sapply(seq(1, length(trialVals)), function(l){
                        
                        targetTrial <- trialVals[l]
                        
                        sapply(seq(1, length(data)), function(k){
                                
                                if(sum(!is.na(relTrialNumbers[[k]]))==0){perfVal <- NA
                                } else {
                                        selIdx <- relTrialNumbers[[k]]==targetTrial
                                        if(sum(selIdx)==0){perfVal <- NA} else {perfVal <- data[[k]][selIdx]}
                                }
                        })
                })
                
                perfRelCP <- colMeans(perfRelCPInd, na.rm=T)
                
                nBins=floor(length(trialVals)/trialBinSize)
                binCuts <- seq(relTrialMin, relTrialMax, by=trialBinSize)
                binIdx <- unlist(lapply(binCuts, function(m){rep(m, trialBinSize)}))
                if(length(binIdx>length(perfRelCP))){
                        binCuts[binCuts==max(binCuts)]<- binCuts[length(binCuts)-1] 
                        binIdx <- unlist(lapply(binCuts, function(m){rep(m, trialBinSize)}))
                }
                
                halfBin <- trialBinSize/2
                xcoord <- binCuts+halfBin
                ycoord <-  sapply(binCuts, function(t){
                        mean(perfRelCP[binIdx==t], na.rm=T)
                })
                
                
                ySEM <- sapply(binCuts, function(t){
                        sd(perfRelCP[binIdx==t], na.rm=T)/sqrt(sum(binIdx==t))
                })
                
                plot.new()
                ymin <- floor(min(ycoord))
                ymax <- ceiling(max(ycoord))
                
                plot.window(xlim=c(relTrialMin, relTrialMax), ylim=c(ymin, ymax))
                
                if(typegraph=="lines"){
                        lines(x=xcoord, y=ycoord, lwd=2, col=color)
                        errBars(x=xcoord, y=ycoord, err=ySEM, color=color)} 
                if(typegraph=="points"){
                        points(x=xcoord, y=ycoord, pch=19, cex=1.5, col=color)
                        errBars(x=xcoord, y=ycoord, err=ySEM, color=color)
                }
                if(typegraph=="both"){
                        lines(x=xcoord, y=ycoord, lwd=2, col=color)
                        points(x=xcoord, y=ycoord, pch=19, cex=1.5, col=color)
                        errBars(x=xcoord, y=ycoord, err=ySEM, color=color)
                }
                
                axis(side=1, at=binCuts)
                axis(side=2, las=2)
                abline(h=0, lty=2)
                abline(v=0, lty=1)
                mtext(side=1, line=3, text="Trial from change point", cex=1.5, font=2)
                mtext(side=2, line=3, text=longNames[j], cex=1.5, font=2)
                
                
                mtext(text=paste(longNames[j], "by", trialBinSize, "trial bin relative to change point"), side=3, line=1, cex=1.5, font=2)
                
                dev.off()
                
        })
        
        
        
}


save(PerformanceFromCP, file="E:/Dropbox/NMDA/R functions/PerformanceFromCP.R")

