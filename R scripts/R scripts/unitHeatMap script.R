
UnitHeatMap <- function(data=masterDF, sessFromCP=0, FRparameters=allNeuronsDS$parameters, folder=BySessFolder, winmin=0, winmax=400, BLmin=-2000, BLmax=0){
        

        # Install and call necessary packages
        if(!require(matlab)){install.packages("matlab")}
        library(matlab)
        
        ## Here's where I select the columns of interest for calculating mean FR in the period defined by winmin winmax
        CueBin <- unique(data$CueBin) #Bin right after cue onset
        binw <- FRparameters$binw #bin width for FR in ms
        
        minBin <- CueBin + winmin/binw #Bin that corresponds to the beginning of the period of interest (defined by winmin)
        maxBin <- CueBin + winmax/binw #Bin that corresponds to the end of the period of interest (defined by winmax)
        WOI <- minBin:maxBin #Window of interest
        
        #BL window
        minBLBin <- CueBin + BLmin/binw #Start of baseline period (in bins)
        maxBLBin <- CueBin + BLmax/binw #End of baseline period (in bins)
        BL_win <- minBLBin:maxBLBin
        
        #Window of interest firing rate columns
        WOI_FRcols <- match(WOI, colnames(data))
        BL_win_FRcols <- match(BL_win, colnames(data))
        
        daySel <- data[data$sessfromCPsess==sessFromCP,]

        dayUnits <- unique(daySel$allUnitIdx) #All of the units recorded on the day defined by sessFromCP
        
        # colValues <- seq(-4, 10, by=0.01)
        # colorpalette = jet.colors(length(valsRange))
        # 
        minFR <- -5 #Anything below this will be categorized as this value
        maxFR <- 6 #Anything above this will be categorized as this value
        colValues <- seq(minFR, maxFR, by=0.1)
        colValues <- round(colValues, 1)
        colorpalette <- colorRampPalette(colors=c("darkblue", "orangered"))(length(colValues))
        
        
        graphname <- paste(folder, "Cue evoked FR per unit on session ", sessFromCP, " from CP.pdf", sep="")
        
        pdf(file = graphname)
        
        plot.new()
        
        xmax <- 40*(sessFromCP+1)
        xmin <- xmax-80
        
        plot.window(xlim=c(xmin, xmax), ylim=c(0, length(dayUnits)+1))

        excRanking <- sapply(seq(1, length(dayUnits)), function(x){
                
                selUnit <- daySel[daySel$allUnitIdx==dayUnits[x],]
                
                FR_WOI <- selUnit[,WOI_FRcols] #Raw firing rate on the bins of interest (Window Of Interest)
                FR_WBl <- selUnit[,BL_win_FRcols] #Raw firing rate on the bins defined as baseline
                
                
                #For some reason FR_WOI is interpreted by R as a list of lists, I need to convert it into a matrix in this roundabout way.
                FR_WOI <- apply(FR_WOI, MARGIN=2, as.numeric)
                FRavg <- rowMeans(FR_WOI)
                
                
                FR_WBl <- apply(FR_WBl, MARGIN=2, as.numeric)
                BLavg <- as.numeric(format(selUnit$BLavg, digits=2))
                BLsd <- as.numeric(format(selUnit$BLsd, digits=2))
                
                # BLavg <- mean(rowMeans(FR_WBl)) #Baseline firing of that unit on that trial
                # BLsd <- sd(rowMeans(FR_WBl))
                # 
                
                FR_WOI_Zsc <- ZscoreCalc(x=FRavg, avg=BLavg, sd=BLsd)
                trialfromCPidx <- selUnit$trialfromCP
                
                #This calculates the criterion I'm going to use to sort neurons in the y axis on my graph
                if(sessFromCP==0){
                        mean(FR_WOI_Zsc[trialfromCPidx<0]) #Mean activity before CP
                } else {
                        mean(FR_WOI_Zsc)
                }
                
        })
        
        #Ranking of the units based 
        excRankingIdx <- (1:length(dayUnits))[order(excRanking)]
        
       datForRegr <- do.call("rbind", lapply(seq(1, length(excRankingIdx)), function(x){
                
                unitIdx <- dayUnits[excRankingIdx][x]
                selUnit <- daySel[daySel$allUnitIdx==unitIdx, ]
                FR_WOI <- selUnit[,WOI_FRcols] #Raw firing rate on the bins of interest (Window Of Interest)
                FR_WBl <- selUnit[,BL_win_FRcols] #Raw firing rate on the bins defined as baseline
                
                
                #For some reason FR_WOI is interpreted by R as a list of lists, I need to convert it into a matrix in this roundabout way.
                FR_WOI <- apply(FR_WOI, MARGIN=2, as.numeric)
                FRavg <- rowMeans(FR_WOI)
                
                
                FR_WBl <- apply(FR_WBl, MARGIN=2, as.numeric)
                BLavg <- as.numeric(format(selUnit$BLavg, digits=2))
                BLsd <- as.numeric(format(selUnit$BLsd, digits=2))
                
                # BLavg <- mean(rowMeans(FR_WBl)) #Baseline firing of that unit on that trial
                # BLsd <- sd(rowMeans(FR_WBl))
                # 
                
                FR_WOI_Zsc <- ZscoreCalc(x=FRavg, avg=BLavg, sd=BLsd)
                trialfromCPidx <- selUnit$trialfromCP
                
                sapply(seq(1, nrow(selUnit)), function(y){
                        valPos <- findInterval(FR_WOI_Zsc[y], colValues)
                        colPick <- colorpalette[valPos]
                        rect(ybottom=x, ytop=x+1, xleft = trialfromCPidx[y], xright=trialfromCPidx[y]+1, col=colPick, border=NULL)
                })
                
                cbind(rat=selUnit$rat, session=selUnit$session, group=selUnit$Session,
                      allUnitIdx=dayUnits[excRankingIdx[x]], TrialFromCP=trialfromCPidx, 
                      Fr_Zsc=FR_WOI_Zsc, Fr_Raw=FRavg, BLavg=BLavg, BLsd=BLsd, 
                      CSplusResponse=selUnit$CSplusResponse, CSplusLat=selUnit$CSplusLat,
                      CSplusSpecif=selUnit$CSplusSpecif, ITIlat=selUnit$ITIlatency
                      )
        }))
       
        datForRegr <- as.data.frame(datForRegr)
                
                abline(v=0, col="yellow")
                axis(side=1, cex.axis=1.4)
                mtext(side=1, line=2.5, text="Trial from CP", font=2, cex=1.5) 
                axis(side=2, las=2, at=seq(1.5, length(dayUnits)+0.5), labels=1:length(dayUnits), cex.axis=0.8)
                mtext(side=2, line=2.5, text="Unit Idx", cex=1.5, font=2)
                title(main=paste("Session ", sessFromCP, " from CP", sep=""))
                
                dev.off()
                
                
                # Scatterplots with lm
                
                PerfIndexes <- c("CSplusLat", "CSplusResponse", "CSplusSpecif", "ITIlat")
        
                if(sessFromCP==0){
                        
                        sapply(seq(1, length(PerfIndexes)), function(i){
                                
                                graphname <- paste(folder, PerfIndexes[i], " vs FR on session ", sessFromCP, " from CP.pdf", sep="")
                                PerfColIdx <- match(PerfIndexes[i], colnames(datForRegr)) #Select column with the beh index of interest
                                
                                pdf(file = graphname)
                                
                                plot.new()
                                
                                perfData <- datForRegr[,PerfColIdx]
                                
                                if(PerfIndexes[i]=="CSplusResponse"){
                                        perfData[!is.na(perfData)] <- 1
                                        perfData[is.na(perfData)] <-0
                                }
                                
                                minx <- floor(min(perfData)); maxx <- ceiling(max(perfData))
                                
                                plot.window(xlim=c(minx, maxx), ylim=c(-4, 30))
                                points(x=perfData[datForRegr$TrialFromCP<0], y=datForRegr$Fr_Zsc[datForRegr$TrialFromCP<0], pch=16, cex=0.8)
                                points(x=perfData[datForRegr$TrialFromCP>=0], y=datForRegr$Fr_Zsc[datForRegr$TrialFromCP>=0], pch=16, cex=0.8, col="blue")
                                axis(side=1, at=seq(0, 10), cex.axis=1.4)
                                mtext(side=1, line=2.5, text=PerfIndexes[i], cex=1.5, font=2)
                                axis(side=2, at=seq(-4, 30, by=2), cex.axis=1.4, las=2)
                                mtext(side=2, line=2.5, text="Firing rate (Z sc.)", cex=1.5, font=2)
                                
                                preCPlm <- lm(perfData[datForRegr$TrialFromCP<0] ~ datForRegr$Fr_Zsc[datForRegr$TrialFromCP<0])
                                postCPlm <- lm(perfData[datForRegr$TrialFromCP>=0] ~ datForRegr$Fr_Zsc[datForRegr$TrialFromCP>=0])
                                
                                abline(preCPlm)
                                abline(postCPlm, col="blue")
                                
                                legend("topleft", paste("Pre CP. Adj. R sq. =", round(summary(preCPlm)[9]$`adj.r.squared`, 5)), lty=1, lwd=2, col="black")
                                legend("topright", paste("Post CP. Adj. R sq. =", round(summary(postCPlm)[9]$`adj.r.squared`, 5)), lty=1, lwd=2, col="blue")
                                
                                abline(h=0, lty=3)
                                
                                title(main=paste("Session ", sessFromCP, " from CP", sep=""))
                                
                                dev.off()
                        })
                        
                
                } else {
                        
                        sapply(seq(1, length(PerfIndexes)), function(i){
                        
                                graphname <- paste(folder, PerfIndexes[i], " vs FR on session ", sessFromCP, " from CP.pdf", sep="")
                                PerfColIdx <- match(PerfIndexes[i], colnames(datForRegr)) #Select column with the beh index of interest
                        
                                pdf(file = graphname)
                        
                                plot.new()
                        
                                perfData <- datForRegr[,PerfColIdx]
                                
                                if(PerfIndexes[i]=="CSplusResponse"){
                                        perfData[!is.na(perfData)] <- 1
                                        perfData[is.na(perfData)] <-0
                                }
                        
                                minx <- floor(min(perfData)); maxx <- ceiling(max(perfData))
                               
                                plot.window(xlim=c(minx, maxx), ylim=c(-4, 30))
                                points(x=perfData, y=datForRegr$Fr_Zsc, pch=16, cex=0.8)
                                axis(side=1, at=seq(0, 10), cex.axis=1.4)
                                mtext(side=1, line=2.5, text=PerfIndexes[i], cex=1.5, font=2)
                                axis(side=2, at=seq(-4, 30, by=2), cex.axis=1.4, las=2)
                                mtext(side=2, line=2.5, text="Firing rate (Z sc.)", cex=1.5, font=2)
                                
                                sessLm <- lm(perfData ~ datForRegr$Fr_Zsc)
                                
                                abline(sessLm)
                                
                                legend("topright", paste("Adj. R sq. =", round(summary(sessLm)[9]$`adj.r.squared`, 5)), lty=1, lwd=2, col="black")
                                
                                abline(h=0, lty=3)
                                
                                title(main=paste("Session ", sessFromCP, " from CP", sep=""))
                                
                                dev.off()
               })
                
        }
        
}


save(UnitHeatMap, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/UnitHeatMap.R")
save(UnitHeatMap, file="E:/Dropbox/NMDA/R Functions/UnitHeatMap.R")
