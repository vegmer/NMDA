#"smartRat" indicates whether I only want to plot cumulative lines for rats that did learn. If so, smartRat=T. If I only want the rats that didn't lear: smartRat=F. If I want all rats: smartRat="NA".
#CP vector: vector with the value of the trials in which the animals made the change point
#limit: I trained some rats for 7 days or more. I may want to discard those latest trials. Limit defines the last trial I want to include in my graph. If limit=NA, then no limit applies


cumulativeIndGraphs <- function(numSess=6, numCues=35, smartRat=T, CPraw=CPrawdata, CPvector=CPdata$CP, sessLines=F, dataForRCumulative, dataForRdir, graphFolder, imgFormat="pdf", limit=210){
        #Names for graph titles. Make sure they are in the same order that these objects are loaded from the folder in the next line
        longNames <- c("S+ latency", "S+ resp ratio", "S+ specificity", "S+ time to spare", "ITI latency", "S- latency", "S- resp ratio", "S- specificity", "S- time to spare")
        
        # Load relevant objects and create object 'cumData' with the labels for each object
        filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); cumData <- sapply(seq(1, length(filesCum)), function(i){load(filesCum[[i]])})
        
        # Plot
        datacumsum <- lapply(seq(1, length(filesCum)), function(j){
                
                #Create object with per trial data for each subject for that particular (j) index
                data <- get(cumData[j])
                
                if(smartRat==T){ ratSelIdx <-!is.na(CPvector); lab="Smart rats"} 
                if(smartRat==F){ ratSelIdx <- is.na(CPvector); lab="Dumb rats"}
                if(smartRat=="NA"){ ratSelIdx <- rep(TRUE, length(CPvector)); lab="All rats"}
                
                #Subset the group of rats I want (smart, dumb or both)
                data <- data[c(ratSelIdx)]
                
                #Cap the number of trials I want for my graphs
                if(!is.na(limit)){
                        if(j!=5){data <- lapply(seq(1, length(data)), function(q){
                                        dat <- data[[q]][1:limit]
                                        NAs <- is.na(dat)
                                        if(sum(NAs)<=10){
                                              datT <- dat[NAs==F]
                                              s <- sample(datT[(length(datT)-15):length(datT)], sum(NAs), replace=T)
                                              dat[is.na(dat)] <- s
                                        }
                                        data <- dat
                                })
                        }
                        
                        if(j==5){data <- lapply(seq(1, length(data)), function(q){data[[q]][1:(2*limit)]})} #j=5 is ITI latency, which includes DS and NS trials. So in this case the limit should be the double of the limit.
                        
                }
                
                #Establish # of cues per session and in total (generalized from 1 animal, if there are individual differences, check)
                #numCues <- length(alldata[[idx[[1]][1]]]$CSpluscue)
                numSess <- numSess
                #totalCues <- numCues*numSess
                totalCues <- max(sapply(data, length))
                
                
                #Cumulative sum
                datacumsum <- lapply(seq(1, length(data)), function(k){cumsum(data[[k]])})
                
                
                if(imgFormat=="png"){filetitle=paste(graphFolder, paste("Cumulative", longNames[j]), lab, ".png", sep=""); png(filename = filetitle)}
                if(imgFormat=="pdf"){filetitle=paste(graphFolder, paste("Cumulative", longNames[j]), lab, ".pdf", sep=""); pdf(file = filetitle)}
                
                plot.new()
                
                #ymin <- min(sapply(datacumsum, min, na.rm=T), na.rm=T)
                #ymax <- max(sapply(datacumsum, max, na.rm=T), na.rm=T)
                
                #I'm hard coding the y axis after to make graphs for the paper (I looked at the maximum values for each variable)
                if("S+ specificity"==longNames[j]){ymin=-400; ymax=800}
                if("S- specificity"== longNames[j]){ymin=-400; ymax=800}
                if(grepl("resp ratio", longNames[j])==T){ymin=0; ymax=300}
                if(grepl("latency", longNames[j])==T){ymin=0; ymax=2500}
                if(grepl("ITI latency", longNames[j])==T){ymin=0; ymax=6000}
                
                
                plot.window(xlim=c(0, totalCues+numCues), ylim=c(ymin, ymax+10))
                
                lapply(seq(1, length(data)), function(k){
                        
                        xcoord=c(1:length(data[[k]]))[!is.na(data[[k]])]
                        ycoord=c(datacumsum[[k]])[!is.na(data[[k]])]
                        lines(x=xcoord, y=ycoord, lwd=2, col='black')
                        ratCode <- LETTERS[1:length(rats)] #Rename rats with letters (cleaner for graphs)
                        ratCodeSel <- ratCode[ratSelIdx]
                        ratlabel=ratCodeSel[k]
                        text(x=max(xcoord+3), y=last(ycoord+5), labels=ratlabel, cex=0.8, col="red", font=2)
                        
                        
                })
                
                abline(h=0, lty=2)
                if(sessLines==T){abline(v=seq(numCues, totalCues, by=numCues), lwd=1, col="purple")} #Lines demarcating sessions
                axis(side=1, at=seq(0, totalCues, by=numCues))
                axis(side=2, las=2)
                mtext(side=1, line=3, text="Trial", cex=1.5, font=2)
                mtext(text=paste("Cumulative", longNames[j]), side=2, line=3, cex=1.5, font=2)
                
                mtext(side=3, text="Individual cumulative performance throughout training")
        
                dev.off() 
                
                return(datacumsum)
        })
        
        #Now, with the good cumsum data, make graphs for individual graphs. I'm mostly interested in graphs that show S+ specificity with respect to CP data:
        objectIdx <- (1:length(longNames))[longNames=="S+ specificity"]
        
        
        #Plot cumulative performance
        if(smartRat==T){ ratSelIdx <-!is.na(CPvector)} 
        if(smartRat==F){ ratSelIdx <- is.na(CPvector)}
        if(smartRat=="NA"){ ratSelIdx <- rep(TRUE, length(CPvector))}
        
        cumDStaskAcc <- datacumsum[[objectIdx]][c(ratSelIdx)]
        
        
        lapply(seq(1, length(cumDStaskAcc)), function(x){
                
                ratname <- rats[ratSelIdx][x]
                
                filename <- paste(CPGraphFolder, ratname, "S+ spec f.pdf", sep="")
                pdf(file = filename, width = 8.5, height = 3.1)
                
                CPrawRat <- CPraw[[x]][[objectIdx]]
                plot(DStaskAcc[[x]], type="l", lwd=1.5, xlab='Trial', ylab="Cumulative S+ spec.")
                abline(h=0, lwd=1.5)
                abline(v=c(0, cumsum(rep(numCues, numSess))), lwd=0.5, col="gray50")
                CP <- CPvector[x]
                points(x=CPrawRat$CP_trial[-1], y=DStaskAcc[[x]][CPrawRat$CP_trial], col="red", pch=19)
                if(!is.na(CP)){abline(v=CP, col="red", lwd=2)}
                
                title(main=ratname, font=2, cex=1.5)
                
                dev.off()
        })
}

save(cumulativeIndGraphs, file=paste("E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R functions/cumulativeIndGraphs.R"))
save(cumulativeIndGraphs, file=paste("E:/Dropbox/NMDA/R functions/cumulativeIndGraphs.R"))
save(cumulativeIndGraphs, file=paste("E:/Dropbox/NMDA/EXP4_Unilateral AP5/R functions/cumulativeIndGraphs.R"))


