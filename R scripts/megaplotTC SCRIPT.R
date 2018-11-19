megaplotTC <- function(experiment, data=list(allNeuronsDS), graphmin=500, graphmax=10000, cuewinmin=100, cuewinmax=400,
                     colpalette="Rainbow", minFR=-2, maxFR=6, graphFolder=MixedGraphFolder, 
                     ZcolLabels=c("ZDS"), arrangeBy=c("ZDS")){
        
        # Install and call necessary packages
        if(!require(colorRamps)){install.packages("colorRamps")}
        if(!require(dplyr)){install.packages("dplyr")}
        library(colorRamps)
        library(dplyr)
        
        
        binw <- data[[1]]$parameters$binw
        psthmin <- data[[1]]$parameters$psthmin
        psthmax <- data[[1]]$parameters$psthmax
        binw <- data[[1]]$parameters$binw
        
        CUEbin <-(psthmin*1000/binw)+1 #First bin after the cue (50ms)
        
        CUEmin <- CUEbin + winmin/binw #Beginning of the window after the cue
        CUEmax <- CUEbin + winmax/binw #End of the window after the cue
        
        BLmin <- CUEbin - (BLwdw*1000)/binw
        BLmax <- CUEbin
        
        graphMinBin <- CUEbin - graphmin/binw
        graphMaxBin <- CUEbin + graphmax/binw
        
        
        unitIdx <- do.call("rbind", sapply(seq(1, length(data[[1]]$nexdata)), function(y){ #In every NEX session I recorded
                rat <- data[[1]]$nexdata[[y]]$ratname
                down <- data[[1]]$nexdata[[y]]$down
                modality <- data[[1]]$nexdata[[y]]$version
                CP <- CPdata$CPsess[rat==CPdata$rat]
                expt <- as.numeric(data[[1]]$nexdata[[y]]$expt)
                CPfromSess <- expt-CP
                
                #Is the session before CP?
                
                if(length(CPfromSess)==0){
                        a <- NA; 
                        CPfromSess <- NA}  
                
                if(is.na(CPfromSess)){a <- NA}
                
                if(length(CPfromSess)!=0 & !is.na(CPfromSess)){
                        if(CPfromSess>0){a <- 1}
                        if(CPfromSess==0){a <- 0}
                        if(CPfromSess<0){a <- -1}
                }      
                
                
                
                BeforeCP <- a #Assign a -1 to before CP sesions, 0 to the session in which CP took place, and 1 to sessions after CP happened
                CSplusTA <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$CSplusTaskAcc, 2)
                CSminusTA <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$CSminusTaskAcc, 2)
                CSplusLat <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$LatencyCSplus, 2)
                CSminusLat <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$LatencyCSminus, 2)
                DSRR <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$DSperc, 2)
                NSRR <- round(csacqidx[csacqidx$subject==rat & csacqidx$session==expt, ]$NSperc, 2)
                
                #units <- names(data[[1]]$nexdata[[y]])[grepl("sig", names(data[[1]]$nexdata[[y]]))]
                units <- 1:length(data[[1]]$neurons[[y]])
                DSExc <- unlist(data[[1]]$cueexidx[[y]])
                DSInh <- unlist(data[[1]]$cueinhidx[[y]])
                
                
                d <- cbind(rat, expt, CPfromSess, BeforeCP, CSplusTA, CSminusTA, CSplusLat, CSminusLat, DSRR, NSRR, units, DSExc, DSInh, down, modality)
                
                if(dim(d)[2] != 15){
                        CSplusTA <- NA; CSminusTA <- NA; CSplusLat <- NA; CSminusLat <- NA; DSRR <- NA; NSRR <- NA
                        d <- cbind(rat, expt, CPfromSess, BeforeCP, CSplusTA, CSminusTA, CSplusLat, CSminusLat, DSRR, NSRR, units, DSExc, DSInh, down, modality)
                }
                
                return(d)
        })
        
        )
        
        #This makes a data frame with post DS firing rate in the cuewinmin:cuewinmax window but also per bin firing rate around the cue (until the boundaries defined by graphmin:graphmax)
        FR_Zsc_CUEandTC <- lapply(seq(1, length(data)), function(x){
                
                fmat <- data[[x]]$firingmat
                
                BLfr <- colMeans(fmat[BLmin:BLmax, ], na.rm = T)
                BLsd <- colSds(fmat[BLmin:BLmax, ], na.rm=T)
                CUEfr <- colMeans(fmat[CUEmin:CUEmax, ], na.rm = T)
                
                ZPerUnitTC <- do.call("rbind", lapply(seq(1, ncol(fmat)), function(y){ #For each unit
                        perUnitZsc <- ZscoreCalc(x=fmat[, y], avg=BLfr[y], sd=BLsd[y])
                })
                )
                
                forGraphZPerUnitTC <- ZPerUnitTC[ , graphMinBin:(graphMaxBin-1)]
                colnames(forGraphZPerUnitTC) <- graphMinBin:(graphMaxBin-1)
                
                CUEfr_Zsc <- as.data.frame(ZscoreCalc(x=CUEfr, avg=BLfr, sd=BLsd))
                colnames(CUEfr_Zsc) <- ZcolLabels
                
                cbind(CUEfr_Zsc, forGraphZPerUnitTC)
                
        })
        
        
        
        toplotTC <- data.frame(unitIdx, FR_Zsc_CUEandTC) #This is the tidy dataset I'll use
        
        save(toplotTC, file=paste(dataForRdir, "toplotTC.rdat", sep=""))
        
        ####################
        ### Megaplots    ###
        ####################
        
        toplot2 <- toplotTC[order(toplotTC$BeforeCP, toplotTC[,colnames(toplotTC)==arrangeBy], na.last = F),]
        breaks <- seq(minFR, maxFR, 0.05)
        
        toplot2 <- toplot2[!is.na(toplot2$ZDS), ]
        
        
        filename=paste(graphFolder, experiment, "Megaplot timecourse", length(data), "events", winmax, "ms postcue.pdf", sep=" ")
        pdf(file = filename, width=9, height=9)
        
        plot.new()
        plot.window(xlim=c(1, length(graphMinBin:graphMaxBin)+30), ylim=c(-10, nrow(toplot2)+10))
        
        #Color palette
        if(colpalette=='Rainbow'){
                colorScale <- matlab.like(length(breaks))
        } else {
                colorScale <- colorRampPalette(colors=c("darkblue", "orangered"))(length(breaks))
        }
        
        
        ### Plot the columns with FR data by bin
        
        #Select columns that contain firing rate data
        FRcolIdx <- grepl("X", colnames(toplot2))
        FRdata <- toplot2[ ,FRcolIdx]
        
        #Define where the box with FR data starts and the width of individual columns
        box1left <- 1
        subboxWidth <- length(graphMinBin:graphMaxBin)
        
        #Make the columns with the colors
        sapply(seq(1, ncol(FRdata)), function(x){ #for each column (bin)
                
                sapply(seq(1, nrow(FRdata)), function(y){ #for each unit
                        
                        FRsel <- FRdata[y,x]
                        
                        if(!is.na(FRsel) & FRsel<minFR){FRsel <- minFR}
                        if(!is.na(FRsel) & FRsel>maxFR){FRsel <- maxFR}
                        
                        
                        colIdx <- findInterval(FRsel, breaks) 
                        colpick <- colorScale[colIdx]
                        
                        #Plot the actual value
                        rect(xleft = x, xright=x+1, 
                             ybottom = y, ytop=y+1, col=colpick, border = NA)
                        
                })
        })
        
        ### Plot behavior
        
        #When I have rats that never learned, their value for "Before CP" is NA. So make sure these sessions are labeled differently. If no not learners, just make 3 divisions (pre CP, CP, post CP)
        if(sum(is.na(toplot2$BeforeCP))==0){
                ydivisions <- c(0, cumsum(rle(as.numeric(toplot2$BeforeCP))$lengths))
                ylabels <- c("Before CP", "CP session", "After CP")
                
                mtext(side=4, line=-10, text="Before CP", cex=1.2, font=2, at=(ydivisions[2]-ydivisions[1])/2)
                mtext(side=4, line=-10, text="CP session", cex=1.2, font=2, at=ydivisions[2]+(ydivisions[3]-ydivisions[2])/2)
                mtext(side=4, line=-10, text="After CP", cex=1.2, font=2, at=ydivisions[3]+(ydivisions[4]-ydivisions[3])/2)
                
        }
        
        if(sum(is.na(toplot2$BeforeCP))>0){
                ytable <- table(toplot2$BeforeCP, exclude = NULL) #NAs come out last but in my data frame are first, so switch order in next line
                ydivisions <- c(0, as.numeric(cumsum(c(ytable[4], ytable[1:3]))))
                ylabels <- c("Non learners", "Before CP", "CPsession", "After CP")
                
                mtext(side=4, line=-10, text="Non learners", cex=1.2, font=2, at=(ydivisions[2]-ydivisions[1])/2)
                mtext(side=4, line=-10, text="Before CP", cex=1.2, font=2, at=ydivisions[2]+(ydivisions[3]-ydivisions[2])/2)
                mtext(side=4, line=-10, text="CP session", cex=1.2, font=2, at=ydivisions[3]+(ydivisions[4]-ydivisions[3])/2)
                mtext(side=4, line=-10, text="After CP", cex=1.2, font=2, at=ydivisions[4]+(ydivisions[5]-ydivisions[4])/2)
                
        }
        
        
        sapply(seq(1, length(ydivisions)), function(q){
                segments(x0=box1left, x1=subboxWidth, y0=ydivisions[q]+1, y1=ydivisions[q]+1)
        })
        
        axis(2, pos=0.5, at=ydivisions, las=2)
        mtext(side=2, line=2.5,  text="Unit #", cex=1.5, font=2)
        mtext(side=1, at=length(graphMinBin:graphMaxBin)/2, line=2.5, text="Time from cue onset (s)", cex=1.5, font=2)

        
       
        xlabels <- seq(-graphmin/1000, graphmax/1000, by=binw/100)
        xat <- seq(1, length(graphMinBin:graphMaxBin), length.out=length(xlabels))
        axis(side=1, at=xat, labels = xlabels)
        abline(v=xat[xlabels==0])
        
        # #Mean latency on the session recorded
        # meanDSlatBin <- graphmin/binw+(as.numeric(as.character((toplot2$CSplusLat)))*1000)/binw
        # 
        # sapply(seq(1, length(meanDSlatBin)), function(s){
        #         rect(xleft = meanDSlatBin[s], xright = meanDSlatBin[s]+1,
        #              ybottom=s, ytop=s+1, pch=19)
        # })
        
        #Cue exc or not
        CUEexcIdx <- (1:nrow(toplot2))[as.logical(toplot2$DSExc)] #Index of cue excited units
        
        sapply(seq(1, length(CUEexcIdx)), function(u){
                rect(xleft=xat[xlabels==0]-2, xright=xat[xlabels==0], 
                     ybottom=CUEexcIdx[u], ytop=CUEexcIdx[u]+1, col="black")
        })
        
        #Cue inh or not
        CUEinhIdx <- (1:nrow(toplot2))[as.logical(toplot2$DSInh)] #Index of cue excited units
        
        sapply(seq(1, length(CUEinhIdx)), function(u){
                rect(xleft=xat[xlabels==0]-4, xright=xat[xlabels==0]-2, 
                     ybottom=CUEinhIdx[u], ytop=CUEinhIdx[u]+1, col="white")
        })
        
       
        
        #Legend
        legendBoxLeft <- box1left+subboxWidth+10
        legendBoxRight <- box1left+subboxWidth+15
        legendsubboxHt <- nrow(toplot2)/length(colorScale)
        
        sapply(seq(1, length(colorScale)), function(c){
                rect(xleft=legendBoxLeft, xright=legendBoxRight, 
                     ybottom=legendsubboxHt*c, ytop=legendsubboxHt*(c+1),
                     col=colorScale[c], border=colorScale[c])
        })
        
        legLabels <- seq(minFR, maxFR, by=1)
        legLabelsFixed <- c(paste("<= ", minFR, sep=""), legLabels[-c(1, length(legLabels))], paste(">= ", maxFR, sep=""))
        stepHt <- nrow(toplot2)/(length(legLabels)-1)
        
        axis(side=4, line=-2, at=seq(1, nrow(toplot2)+1, by=stepHt), labels=legLabelsFixed, las=2)
        mtext(side=4, text="Firing rate (Z sc.)", font=2, cex=1.5)
        
        legend("bottomright", pch=22, pt.bg=c("white", "black"), legend=c("Cue inhibited", "Cue excited"))
        
        
        dev.off()
        
        return(toplotTC)
        
}




save(megaplotTC, file="E:/Dropbox/NMDA/R Functions/megaplotTC.R")
