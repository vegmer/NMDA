megaplot <- function(data=list(allNeuronsDSresponded, allNeuronsDSmissed, allNeuronsNSresponded, allNeuronsNSmissed), CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=0, winmax=400,
                     colpalette="Rainbow", minFR=-2, maxFR=6, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir,
                     ZcolLabels=c("ZDSresp", "ZDSmissed", "ZNSresp", "ZNSmissed"), arrangeBy=c("ZDSresp")){
        
        # Install and call necessary packages
        if(!require(colorRamps)){install.packages("colorRamps")}
        if(!require(dplyr)){install.packages("dplyr")}
        library(colorRamps)
        library(dplyr)
        
        
        psthmin <- data[[1]]$parameters$psthmin
        psthmax <- data[[1]]$parameters$psthmax
        binw <- data[[1]]$parameters$binw
        
        CUEbin <-(psthmin*1000/binw)+1 #First bin after the cue (50ms)
        
        CUEmin <- CUEbin + winmin/binw
        CUEmax <- CUEbin + winmax/binw
        
        BLmin <- CUEbin - (BLwdw*1000)/binw
        BLmax <- CUEbin
        
      
        unitIdx <- do.call("rbind", sapply(seq(1, length(data[[1]]$nexdata)), function(y){
                rat <- data[[1]]$nexdata[[y]]$ratname
                CP <- CPdata$CPsess[rat==CPdata$rat]
                expt <- as.numeric(data[[x]]$nexdata[[y]]$expt)
                CPfromSess <- expt-CP
                
                #Is the session before CP?
                        if(CPfromSess>0){a <- 1}
                        if(CPfromSess==0){a <- 0}
                        if(CPfromSess<0){a <- -1}
                
                BeforeCP <- a
                CSplusTA <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$CSplusTaskAcc, 2)
                CSminusTA <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$CSminusTaskAcc, 2)
                CSplusLat <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$LatencyCSplus, 2)
                CSminusLat <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$LatencyCSminus, 2)
                DSRR <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$DSperc, 2)
                NSRR <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$NSperc, 2)
                
                units <- names(data[[x]]$nexdata[[y]])[grepl("sig", names(data[[x]]$nexdata[[y]]))]
                
                cbind(rat, expt, CPfromSess, BeforeCP, CSplusTA, CSminusTA, CSplusLat, CSminusLat, DSRR, NSRR, units)
        })
        )
        
        
        CUEfr_Zsc <- sapply(seq(1, length(data)), function(x){
                
                fmat <- data[[x]]$firingmat
                
                BLfr <- colMeans(fmat[BLmin:BLmax, ])
                BLsd <- colSds(fmat[BLmin:BLmax, ])
                CUEfr <- colMeans(fmat[CUEmin:CUEmax, ])
                
                CUEfr_Zsc <- ZscoreCalc(x=CUEfr, avg=BLfr, sd=BLsd)
                
        })
        
        colnames(CUEfr_Zsc) <- ZcolLabels
        
        toplot <- data.frame(unitIdx, CUEfr_Zsc) #This is the tidy dataset I'll use
        
        save(toplot, file=paste(dataForRdir, "toplot.rdat", sep=""))
        
        ####################
        ### Megaplots    ###
        ####################

        toplot2 <- toplot[order(toplot$BeforeCP, toplot[,colnames(toplot)==arrangeBy], na.last = F),]
        breaks <- seq(minFR, maxFR, 0.05)
        
        filename=paste(graphFolder, "Megaplot", length(data), "events", winmax, "ms postcue.pdf", sep=" ")
        pdf(file = filename, width=9, height=9)
        
        plot.new()
        plot.window(xlim=c(0, 10), ylim=c(-10, nrow(toplot2)+10))
        
        #Color palette
        if(colpalette=='Rainbow'){
                colorScale <- matlab.like(length(breaks))
        } else {
                colorScale <- colorRampPalette(colors=c("darkblue", "orangered"))(length(breaks))
        }
        
        
        ### Plot the columns with FR data
        
        #Select columns that contain firing rate data
        FRcolIdx <- grepl("Z", colnames(toplot2))
        FRdata <- toplot2[ ,FRcolIdx]
        
        #Define where the box with FR data starts and the width of individual columns
        box1left <- 1
        subboxWidth <- 1
        
        #Make the columns with the colors
        sapply(seq(1, ncol(FRdata)), function(x){
                
                sapply(seq(1, nrow(FRdata)), function(y){
                        
                        FRsel <- FRdata[y,x]
                        
                        if(!is.na(FRsel) & FRsel<minFR){FRsel <- minFR}
                        if(!is.na(FRsel) & FRsel>maxFR){FRsel <- maxFR}
                        
                        
                        colIdx <- findInterval(FRsel, breaks) 
                        colpick <- colorScale[colIdx]
                        if(x<3){
                                rect(xleft = box1left+(subboxWidth*(x-1)), xright=box1left+(subboxWidth*(x)),
                                     ybottom = y, ytop=y+1, col=colpick, border = NA)
                        }
                        if(x>=3){
                                rect(xleft = box1left+(subboxWidth*(x-1))+0.2, xright=box1left+(subboxWidth*(x))+0.2,
                                     ybottom = y, ytop=y+1, col=colpick, border = NA)
                        }
                        
                })
        })
        
        ### Plot behavior
        ydivisions <- c(0, cumsum(rle(as.numeric(toplot2$BeforeCP))$lengths))
        ylabels <- c("Before CP", "CP session", "After CP")
        sapply(seq(1, length(ydivisions)), function(q){
                segments(x0=box1left, x1=box1left+5, y0=ydivisions[q]+1, y1=ydivisions[q]+1)
                })
        
        axis(2, pos=0.5, at=ydivisions, las=2)
        mtext(side=2, text="Unit #", cex=1.5, font=2)
        
        if(length(data)==4){
                axis(1, pos=0, at=c(box1left+subboxWidth/2, box1left+subboxWidth/2+1), labels=c("S+ resp.", "S+ missed"), las=2, cex.axis=1.4, font=2)
                axis(1, pos=0, at=c(box1left+subboxWidth/2+2+0.2, box1left+subboxWidth/2+3+0.2), labels=c("S- resp.", "S- missed"), las=2, cex.axis=1.4, font=2)
                
        }
        
        if(length(data)==2){
                axis(1, pos=0, at=c(box1left+subboxWidth/2, box1left+subboxWidth/2+1), labels=c("S+", "S-"), las=2, cex.axis=1.4, font=2)
        }
        
        mtext(side=4, line=-17, text="Before CP", cex=1.2, font=2, at=(ydivisions[2]-ydivisions[1])/2)
        mtext(side=4, line=-17, text="CP session", cex=1.2, font=2, at=ydivisions[2]+(ydivisions[3]-ydivisions[2])/2)
        mtext(side=4, line=-17, text="After CP", cex=1.2, font=2, at=ydivisions[3]+(ydivisions[4]-ydivisions[3])/2)
        
        #Legend
        legendBoxLeft <- box1left+7
        legendBoxRight <- box1left+7.3
        legendsubboxHt <- nrow(toplot2)/length(colorScale)
        
        sapply(seq(1, length(colorScale)), function(c){
                rect(xleft=legendBoxLeft, xright=legendBoxRight, 
                     ybottom=legendsubboxHt*c, ytop=legendsubboxHt*(c+1),
                     col=colorScale[c], border=colorScale[c])
        })
        
        legLabels <- seq(minFR, maxFR, by=1)
        legLabelsFixed <- c(paste("<= ", minFR, sep=""), legLabels[-c(1, length(legLabels))], paste(">= ", maxFR, sep=""))
        stepHt <- nrow(toplot2)/(length(legLabels)-1)
        
        axis(side=4, line=-6.5, at=seq(1, nrow(toplot)+1, by=stepHt), labels=legLabelsFixed, las=2)
        
        dev.off()
        
}




save(megaplot, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/megaplot.r")
save(megaplot, file="E:/Dropbox/NMDA/R Functions/megaplot.r")
