cumulativeIndGraphs <- function(numSess=7, sessLines=F, dataForRCumulative, dataForRdir, graphFolder, imgFormat="pdf"){
        #Names for graph titles. Make sure they are in the same order that these objects are loaded from the folder in the next line
        longNames <- c("S+ latency", "S+ resp ratio", "S+ specificity", "S+ time to spare", "ITI latency", "S- latency", "S- resp ratio", "S- specificity", "S- time to spare")
        
        # Load relevant objects and create object 'cumData' with the labels for each object
        filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); cumData <- sapply(seq(1, length(filesCum)), function(i){load(filesCum[[i]])})
        
        # Plot
        lapply(seq(1, length(filesCum)), function(j){
                
                #Create object with per trial data for each subject for that particular (j) index
                data <- get(cumData[j])
                
                #Establish # of cues per session and in total (generalized from 1 animal, if there are individual differences, check)
                numCues <- length(alldata[[idx[[1]][1]]]$CSpluscue)
                numSess <- numSess
                #totalCues <- numCues*numSess
                totalCues <- max(sapply(data, length))
                
                
                #Cumulative sum
                datacumsum <- lapply(seq(1, length(data)), function(k){cumsum(data[[k]])})
                
                if(imgFormat=="png"){filetitle=paste(graphFolder, paste("Cumulative", longNames[j]), ".png", sep=""); png(filename = filetitle)}
                if(imgFormat=="pdf"){filetitle=paste(graphFolder, paste("Cumulative", longNames[j]), ".pdf", sep=""); pdf(file = filetitle)}
                
                plot.new()
                
                ymin <- min(sapply(datacumsum, min))
                ymax <- max(sapply(datacumsum, max))
                
                #I'm hard coding the y axis after to make graphs for the paper (I looked at the maximum values for each variable)
                if("S+ specificity"==longNames[j]){ymin=-200; ymax=800}
                if("S- specificity"== longNames[j]){ymin=-800; ymax=200}
                if(grepl("resp ratio", longNames[j])==T){ymin=0; ymax=300}
                if(grepl("latency", longNames[j])==T){ymin=0; ymax=2500}
                if(grepl("ITI latency", longNames[j])==T){ymin=0; ymax=6000}
                
                
                plot.window(xlim=c(0, totalCues+20), ylim=c(ymin, ymax+10))
                
                lapply(seq(1, length(data)), function(k){
                        
                        xcoord=c(1:length(data[[k]]))
                        ycoord=c(datacumsum[[k]])
                        lines(x=xcoord, y=ycoord, lwd=2, col='black')
                        text(x=max(xcoord+15), y=last(ycoord+7), labels=rats[[k]], cex=1.1, col="red", font=2)
                        
                        abline(h=0, lty=2)
                        if(sessLines==T){abline(v=seq(numCues, totalCues, by=numCues), lwd=1, col="purple")} #Lines demarcating sessions
                        axis(side=1, at=seq(0, totalCues, by=numCues))
                        axis(side=2, las=2)
                        mtext(side=1, line=3, text="Trial", cex=1.5, font=2)
                        mtext(text=paste("Cumulative", longNames[j]), side=2, line=3, cex=1.5, font=2)
                        
                })
                mtext(side=3, text="Individual cumulative performance throughout training")
                
                dev.off()
                
        })
}

save(cumulativeIndGraphs, file=paste("E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R functions/cumulativeIndGraphs.R"))
save(cumulativeIndGraphs, file=paste("E:/Dropbox/NMDA/R functions/cumulativeIndGraphs.R"))


