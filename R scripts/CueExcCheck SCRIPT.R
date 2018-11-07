CueExcCheck <- function(experiment, data=toplotTC, graphFolder=MixedGraphFolder, minFR=-5, maxFR=15, autom.FRrange=FALSE, FRbin=1){
        
        # Install and call necessary packages
        if(!require(colorRamps)){install.packages("colorRamps")}
        if(!require(dplyr)){install.packages("dplyr")}
        library(colorRamps)
        library(dplyr)
        
        
        
        if(autom.FRrange==TRUE){
                minFR <- round(min(data$ZDS, na.rm=T), 0)
                maxFR <- round(max(data$ZDS, na.rm=T), 0)
        }
        
        breaks <- seq(minFR, maxFR, by=FRbin)
               
        hcounts <- findInterval(data$ZDS, breaks) 
        
        dataForPlot <- data.frame(BeforeCP=data$BeforeCP, DSexc=data$DSExc, DSinh=data$DSInh, ZDS=data$ZDS)
        dataForPlot$FRbin <- hcounts
        
        dataForPlot <- arrange(dataForPlot, DSexc, FRbin, BeforeCP)
        
        dataForPlot$idx <- 1:nrow(dataForPlot)
        
        
        
        
        ####plot
        
        filename=paste(graphFolder, experiment, "Distribution of Post S+ FR and Cue Exc.pdf", sep="")
        pdf(file = filename)
        
        plot.new()
        plot.window(xlim=c(min(breaks), max(breaks)), ylim=c(0, max(table(hcounts))))
        
        abline(v=0, lwd=2)
        
        chunks <- unique(dataForPlot$Before)
        
        sapply(seq(1, length(unique(chunks))), function(x){
                
                dataSel <- dataForPlot[dataForPlot$BeforeCP==chunks[x], ]
                if(chunks[x]==-1){colpick <- "gray50"}
                if(chunks[x]==0){colpick <- "gray25"}
                if(chunks[x]==1){colpick <- "black"}
                
                uniqFRbin <- unique(dataSel$FRbin)
                
                sapply(seq(1, length(uniqFRbin)), function(y){
                        
                        dataSelBin <- dataSel[dataSel$FRbin==uniqFRbin[y], ]
                        
                        sapply(seq(1, nrow(dataSelBin)), function(z){
                                
                                unitInfo <- dataSelBin[z, ]
                                if(unitInfo$DSexc==TRUE){pt.col=colpick} else {pt.col="white"}
                                
                                forYpos<- dataForPlot[dataForPlot$idx<unitInfo$idx, ]
                                ypos <- sum(forYpos$FRbin==uniqFRbin[y])+1
                                
                                points(x=breaks[unitInfo$FRbin], y=ypos, bg=pt.col, pch=21, cex=1.5)
                                
                                if(unitInfo$DSinh==TRUE){
                                        points(x=breaks[unitInfo$FRbin], y=ypos, pch=19, col="blue", cex=0.5)
                                }
                        })
                        
                })
        })
        
        axis(side=1, at=breaks)
        mtext(side=1, text="Firing rate 100-400 post S+ (Z sc.)", line=2.5, cex=1.4, font=2)
        axis(side=2, las=2)
        mtext(side=2, text="# Units", line=2.5, cex=1.4, font=2)
        
        legend("topright", pch=21, pt.bg=c("white", "gray50", "gray25", "black", "blue"), col=c("black", "black", "black", "black", "blue"),
               legend=c("Non S+ exc.", "S+ exc. Pre CP", "S+ exc. CP sess.", "S+ exc. Post CP", "S+ inh."), pt.cex=c(2, 2, 2, 2, 1))
        
        
        dev.off()
}




save(CueExcCheck, file="E:/Dropbox/NMDA/R Functions/CueExcCheck.R")
