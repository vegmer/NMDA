### First I need to load the objects: alldata, csacqidx, idx and rats that are generated with teh function "MEDPCextract".
# Select the rat (by name, e.g. "MV175") and the day (in number e.g. 2, which would be the second session).
# Indicate the limits for the raster plot (e.g. 20 s before and after the cue would be: lowerLimit=20, upperLimit=20)
# If combined=TRUE, one single raster with DS and NS trials is shown
# If combined=FALSE, two separate rasters are shown

behRASTERS <- function(subject, day, lowerLimit=20, upperLimit=20, combined=FALSE,
                         graphFolder=rasterGraphFolder){
        
        ratIndex <- subject==rats
        sessionIndex <- idx[c(ratIndex)][[1]][day]
        
        #Subset the behavioral data of the session of interest
        toPlot <- alldata[[sessionIndex]]
        
        entries <- toPlot$receptacleentries
        exits <- toPlot$receptacleexits
        
        entryData <- data.frame(entries=entries, exits=exits)
        
        if(entries[1] > exits[1]){entries <- c(0, entries)} #If for some reason the first event is an exit and not an entry (perhaps the animal was inside when the session started), add an entry at second 0
        if(last(exits) < last(entries)){exits <- c(exits, last(toPlot$allCueEnds))} #If the session ends and the animal is still inside, no exit is recorded. Assign the end of the session (end of last cue) as the last exit
        
        #Establish the time periods around the cue that will be represented in the rasters
        DSwindows <- data.frame(winStart=toPlot$CSpluscue-lowerLimit, 
                                cue=toPlot$CSpluscue, 
                                winEnd=toPlot$CSpluscue+upperLimit,
                                cueIdx=1:length(toPlot$CSpluscue),
                                AllCueIdx=(1:length(toPlot$orderCues))[toPlot$orderCues==1],
                                latency=toPlot$CSplusLat
                                )
        
        NSwindows <- data.frame(winStart=toPlot$CSminuscue-lowerLimit, 
                                cue=toPlot$CSminuscue, 
                                winEnd=toPlot$CSminuscue+upperLimit,
                                cueIdx=1:length(toPlot$CSminuscue),
                                AllCueIdx=(1:length(toPlot$orderCues))[toPlot$orderCues==2],
                                latency=toPlot$CSminusLat
        )
        
        #Make a list with these two objects: the window data for DS trials and for NS trials
        windows <- list(DSwindows, NSwindows); names(windows) <- c("DS", "NS")
        
        #Assign entries and exits to each window
        
        forRaster <- lapply(seq(1, length(windows)), function(c){
                
                #Select the DS or NS data
                
                windowsSel <- windows[[c]]
                
                do.call("rbind", lapply(seq(1, nrow(entryData)), function(x){
                        
                        entry <- entryData[x,1]
                        exit <- entryData[x,2]
                        
                        #Does that entry take place in any of the depicted windows
                        anyWindow <- entry > windowsSel$winStart & entry < windowsSel$winEnd
                        
                        if(sum(anyWindow)>0){
                                a <- windowsSel[anyWindow, ]
                                
                                if(c==1){cueDuration <- a$latency}
                                if(c==2){cueDuration <- 10}
                                
                                #Did the rat exit before the end of that window?
                                sameWindowExit <- exit < a$winEnd
                                
                                if(sum(sameWindowExit)>0){
                                        output <- cbind(a, data.frame(entry=entry, exit=exit, 
                                                                      cueRaster=a$cue-a$cue, cueEndRaster=cueDuration,
                                                                      entryRaster=entry-a$cue, exitRaster=exit-a$cue))
                                } else {
                                        output <- cbind(a, data.frame(entry=entry, exit=exit,
                                                                      cueRaster=a$cue-a$cue, cueEndRaster=cueDuration,
                                                                      entryRaster=entry-a$cue, exitRaster=a$winEnd-a$cue))
                                }
                                
                        } else {
                                #It is possible that the entry happens before the beginning of the window but that the rat is still inside when the window starts. 
                                #In that case, check if the exit of an entry that is NOT inside the windows may have an exit inside the window
                                
                                sameWindowExit <- exit > windowsSel$winStart & exit < windowsSel$winEnd
                                
                                a <- windowsSel[sameWindowExit, ]
                        
                                if(c==1){cueDuration <- a$latency}
                                if(c==2){cueDuration <- 10}
                                
                                if(sum(sameWindowExit)>0){
                                        output <- cbind(a, data.frame(entry=entry, exit=exit,
                                                                      cueRaster=a$cue-a$cue, cueEndRaster=cueDuration,
                                                                      entryRaster=a$winStart-a$cue, exitRaster=exit-a$cue))
                                } else {output <- c()} #If the entry was made outside of any of the windows depicted in the raster, ignore
                        } 
                        
                        output
                        
                })
                )
                
                
                
        })
                
                
        names(forRaster) <- c("S+", "S-")       
        
        
        #### PLOT THE RASTERS
        
        if(combined==FALSE){
                
                pdf(file=paste(graphFolder, "Beh raster ", subject, "day ", day, " Separate trials.pdf", sep=""))
                
                par(mfrow=c(1, 2))
                
                sapply(seq(1, length(forRaster)), function(c){
                        
                        plot.new()
                        
                        toplot <- forRaster[[c]]
                        
                        plot.window(xlim=c(-lowerLimit, upperLimit), ylim = c(1, nrow(windows[[c]])))
                        
                        abline(v=0, col="red", lwd=2)
                        abline(v=-10, col="red", lty=2)
                        abline(v=10, col="darkred", lwd=1.5)
                        
                        #rect(xleft = toplot$cueRaster,
                        #     xright = toplot$cueEndRaster,
                        #     ybottom = toplot$cueIdx-0.5,
                        #     ytop = toplot$cueIdx+0.5,
                        #     col="#2171b5",
                        #     border="white"
                        #     )
                        
                        segments(x0=toplot$entryRaster,
                                 y0=toplot$cueIdx,
                                 x1=toplot$exitRaster,
                                 y1=toplot$cueIdx)
                        
                        if(c==1){
                                sapply(seq(1, nrow(toplot)), function(i){
                                        if(toplot$cueEndRaster[i] > 0 & toplot$cueEndRaster[i] < 10){
                                                points(x=toplot$cueEndRaster[i], y=toplot$cueIdx[i], pch=16, col="darkblue")
                                        }
                                })
                        }
                        
                        axis(side=1, seq(-lowerLimit, upperLimit, 5), cex=1.2)
                        axis(side=2, seq(0, nrow(toplot), by=5), cex=1.2, las=2)
                        
                        mtext(side=1, line=2.5, cex=1.4, font=2, text="Time from cue onset (s)")
                        mtext(side=2, line=2.5, cex=1.4, font=2, text="Trial")
                        
                        title(main=paste(subject, "Day", day, names(forRaster)[c], sep=" "))
                        
                        
                })
                
                dev.off()
                
        }
        
        
        
        if(combined==TRUE){
                
                pdf(file=paste(graphFolder, "Beh raster ", subject, "day ", day, " ALL trials.pdf", sep=""))
                
                plot.new()
                plot.window(xlim=c(-lowerLimit, upperLimit), ylim = c(1, nrow(windows[[1]])+nrow(windows[[2]])))
                
                abline(v=0, col="red", lwd=2)
                abline(v=-10, col="red", lty=2)
                abline(v=10, col="darkred", lwd=1.5)
          
                
                #colorPick <- c("#2171b5", "darkblue")
                colorPick <- c("black", "gray60") 
                    
                sapply(seq(1, length(forRaster)), function(c){
                        segments(x0=forRaster[[c]]$entryRaster,
                                 y0=forRaster[[c]]$AllCueIdx,
                                 x1=forRaster[[c]]$exitRaster,
                                 y1=forRaster[[c]]$AllCueIdx,
                                 col=colorPick[c])
                        
                        if(c==1){
                                sapply(seq(1, nrow(forRaster[[c]])), function(i){
                                        if(forRaster[[c]]$cueEndRaster[i] > 0 & forRaster[[c]]$cueEndRaster[i] < 10){
                                                points(x=forRaster[[c]]$cueEndRaster[i], y=forRaster[[c]]$AllCueIdx[i], pch=16, col="darkblue")
                                        }
                                })
                        }
                        
                })
                        
                axis(side=1, seq(-lowerLimit, upperLimit, 5), cex=1.2)
                axis(side=2, seq(0, nrow(windows[[1]])+nrow(windows[[2]]), by=5), cex=1.2, las=2)
                        
                mtext(side=1, line=2.5, cex=1.4, font=2, text="Time from cue onset (s)")
                mtext(side=2, line=2.5, cex=1.4, font=2, text="Trial")
                        
                title(main=paste(subject, "Day", day, "S+ and S-", sep=" "))
                        
                legend("topleft", lty=1, col=colorPick, legend=c("S+", "S-"), cex=1)
                
                dev.off()
                
        }
        
        
        
}


save(behRASTERS, file="E:/Dropbox/NMDA/R functions/behRASTERS.R")

