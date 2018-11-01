CPdayFR <- function(experiment="Exp 4", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                    comp=c("VEH", "AP5"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                    correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, 
                    yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=0, WdwEnd=400, 
                    removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T){
        
        
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        
        if(imgFormat=="pdf"){pdf(file=paste(graphFolder, experiment, "FR Before vs After CP Boxplot", "Sessions", sessFromCP, dataProcess, trialSel, ".pdf", sep="_"))}
        if(imgFormat=="png"){png(filename=paste(graphFolder, experiment, "FR Before vs After CP Boxplot", Sessions, sessFromCP, dataProcess, trialSel, ".png", sep="_"))}
        
        binw <- neudata$parameters$binw
        minBin <- WdwStart/binw
        maxBin <- WdwEnd/binw
        
        selBins <- minBin:maxBin
        
        
        plot.new()
        
        # This function has 2 functions: 
        # a) Calculate and spit the mean FR per bin around the time of the event for each session (w respect to change point) for each group of units (VEH vs AP5)
        # b) Plot that info
        
        FRbyUnitBoth <- lapply(seq(1, length(masterDF)), function(c){
               
                CPdaySel <- filter(masterDF[[c]], sessfromCPsess==0)
                
                CPdaySel$BeforeCP <- ((CPdaySel$trialfromCP)>0)*(1) #0 is trials before CP, 1 is trials after CP
                PrePostCPidx <- unique(CPdaySel$BeforeCP)
                
                meanFRWOI <- sapply(seq(1, length(PrePostCPidx)), function(i){
                        
                        dataSel <- filter(CPdaySel, BeforeCP==PrePostCPidx[i]) #Trials before or after the CP
                        
                        FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                        
                        subFRcolNames <- (unique(masterDF[[c]]$CueBin)+minBin):(unique(masterDF[[c]]$CueBin)+maxBin) #Select bins to be plotted
                        subFRcols <- colnames(masterDF[[c]]) %in% as.character(subFRcolNames)
                        
                        ZscoreCalc <- function(x, avg, sd){(x-avg)/sd}
                        
                        if(correctOnly==TRUE){dataSel <- dataSel[!is.na(dataSel$CueResponse), ]}
                        
                        if(cueExcOnly==TRUE){dataSel <- dataSel[dataSel$CueExcited==T, ]}
                        
                        if(sum(is.na(dataSel[,1]))!=nrow(dataSel)){ #If no rows are left after the filters I just applied, then ignore the following code. Only apply if there are units to apply it to
                                
                                #All the units recorded on that session
                                uniqUnits <- unique(dataSel$allUnitIdx)
                                
                                byUnit <- do.call("rbind", lapply(seq(1, length(uniqUnits)), function(u){
                                        unitSel <- filter(dataSel, allUnitIdx==uniqUnits[u])
                                        numericDF <- apply(unitSel[, subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        BLaverage <- as.numeric(format(unique(unitSel$BLavg), digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                        BLsd <- as.numeric(format(unique(unitSel$BLsd), digits=2))
                                        if(is.null(nrow(numericDF))){
                                                MeanByUnit <- mean(numericDF, na.rm=T); 
                                                MeanByUnitZsc <-  ZscoreCalc(x=MeanByUnit, avg=BLaverage, sd=BLsd)
                                                
                                                } else {
                                                
                                                MeanByBin <- colMeans(numericDF, na.rm=T)
                                                MeanByUnit <- mean(MeanByBin, na.rm=T)
                                                MeanByUnitZsc <- ZscoreCalc(x=MeanByUnit, avg=BLaverage, sd=BLsd)
                                        }
                                        MeanByUnitZsc <- ZscoreCalc(x=MeanByUnit, avg=BLaverage, sd=BLsd)
                                        CueExcited <- unitSel$CueExcited[1]
                                        m <- data.frame(Unit=uniqUnits[u], FRbyUnit=MeanByUnit, FRZsc=MeanByUnitZsc, CueExcited=CueExcited)
                                        m
                                        return(m)
                                })
                                )
                                
                                if(cueExcOnly==T){
                                        byUnit <- filter(byUnit, CueExcited==T)
                                }
                                
                                
                                if(dataProcess=="Zscores"){
                                        
                                        MeanByUnit <- byUnit$FRZsc
                                        
                                        yAxMax=yAxMaxZ
                                        labelLeg="(Z sc.)"
                                        
                                } else {
                                        MeanByUnit <- byUnit$FRbyUnit
                                        
                                        yAxMax=yAxMaxRaw
                                        labelLeg="(Hz)"
                                }
                                
                                plot.window(xlim=c(0, length(PrePostCPidx)+1), ylim=c(yAxMin, yAxMax+3))
                                
                                MeanByUnit <- MeanByUnit[!is.nan(MeanByUnit)]
                                
                                barSide <- (i-2)+(i-1) #This will put PRE CP side to the left and POST CP side to the right
                                
                                Q1 <- summary(MeanByUnit)[2]
                                Q3 <- summary(MeanByUnit)[5]
                                IQR <- IQR(MeanByUnit)
                                Median <- summary(MeanByUnit)[3]
                                
                                #IQR rectangle
                                rect(xleft=c+(barSide)*0.3, xright=c, ybottom=Q1, ytop = Q3, col = colindx[c], border="white")
                                
                                #Median line
                                segments(x0=c+(barSide)*0.3, x1=c, y0=Median, y1=Median, lwd=2)
                                segments(x0=c+(barSide)*0.3, x1=c, y0=mean(MeanByUnit), y1=mean(MeanByUnit), lwd=2, col = "white")
                                
                                if(morethanIQR==T){
                                        #Whiskers: maximum value still within Q3+1.5*IQR (whatever is smaller) or minimum value Q1-1.5*IQR
                                        overTop <- MeanByUnit>(Q3+1.5*IQR); top <- max(MeanByUnit[overTop==F])
                                        underBottom <- MeanByUnit<(Q1-1.5*IQR); bottom <- min(MeanByUnit[underBottom==F])
                                        topWhisker <- min(max(MeanByUnit), top)
                                        bottomwhisker <- max(min(MeanByUnit), bottom)
                                        
                                        segments(x0=c+barSide*0.15, x1=c+barSide*0.15, y0=Q3, y1=topWhisker)
                                        segments(x0=c+barSide*0.15, x1=c+barSide*0.15, y0=Q1, y1=bottomwhisker)
                                        
                                        overWhisker <- MeanByUnit[overTop]
                                        underWhisker <- MeanByUnit[underBottom]
                                        
                                        #Outliers
                                        points(x=rep(c+((barSide)*0.15), length(overWhisker)), y=overWhisker, cex=0.2, pch=19)
                                        points(x=rep(c+((barSide)*0.15), length(underWhisker)), y=underWhisker, cex=0.2, pch=19)
                                        
                                }
                                
                                if(removeOutliers==T){
                                        outlierIdx <- (1:length(MeanByUnit))[(overTop==T | underBottom==T)]
                                        if(length(outlierIdx)>0){MeanByUnit <- MeanByUnit[-outlierIdx]}
                                }
                                
                        }
                        
                        return(MeanByUnit)
                        })
                })
        
                sapply(seq(1, length(FRbyUnitBoth)), function(x){
                        xpos <- c(x-0.1, x+0.1)
                        sapply(seq(1, nrow(FRbyUnitBoth[[x]])), function(u){
                                lines(x=xpos, y=FRbyUnitBoth[[x]][u, ])
                        })
                        
                        wilcox.test(FRbyUnitBoth[[x]][,1], FRbyUnitBoth[[x]][,2], paired=T)
                })
                
                axis(side=1, at=seq(1, length(PrePostCPidx)), labels=comp, cex.axis=1.4, tick = F)
                
                #Add axis, labels and legend
                if(dataProcess=="Zscores"){yAxMax=yAxMaxZ; yAxMin=yAxMinZ}
                if(dataProcess=="raw"){yAxMax=yAxMaxRaw; yAxMin=yAxMinRaw}
                
                axis(side=2, at=seq(yAxMin, yAxMax, by=2), las=2, cex.axis=1.4, pos=0.6)
                
                mtext(side=2, line=2.5, cex=1.5, font=2, text=paste("Firing rate", labelLeg, sep=" "))
                        
                     
                        
                        
                }
                