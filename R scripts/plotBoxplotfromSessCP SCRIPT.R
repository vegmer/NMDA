

plotBoxplotfromSessCP <- function(experiment="Exp 4", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                               comp=c("VEH", "AP5"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                               correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                               yAxMinZ = -2, yAxMaxZ = 12, yAxMaxRaw = 10, WdwStart=0, WdwEnd=400, 
                               factors=list(drugs=c("VEH", "AP5")), lines=FALSE, analysis=FALSE,
                               removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T){
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        if(!require(matrixStats)){install.packages("matrixStats")}
        library(matrixStats)
        if(!require(ez)){install.packages("ez")}
        library(ez)
        
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        
        if(imgFormat=="pdf"){pdf(file=paste(graphFolder, experiment, "FR Before vs After CP Boxplot", "Sessions", sessFromCP, dataProcess, trialSel, ".pdf", sep="_"))}
        if(imgFormat=="png"){png(filename=paste(graphFolder, experiment, "FR Before vs After CP Boxplot", Sessions, sessFromCP, dataProcess, trialSel, ".png", sep="_"))}
        
        binw <- neudata$parameters$binw
        minBin <- WdwStart/binw
        maxBin <- WdwEnd/binw
        
        selBins <- minBin:maxBin
        
        sessIdx <- sessFromCP[1]:sessFromCP[2]
        
        plot.new()
        
        # This function has 2 functions: 
        # a) Calculate and spit the mean FR per bin around the time of the event for each session (w respect to change point) for each group of units. Also, the second eleemnt of each list is the unit Index (for later)
        # b) Plot that info
        
        FRbyUnitBoth <- lapply(seq(1, length(masterDF)), function(c){
                
                
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                
                subFRcolNames <- (unique(masterDF[[c]]$CueBin)+minBin):(unique(masterDF[[c]]$CueBin)+maxBin) #Select bins to be plotted
                subFRcols <- colnames(masterDF[[c]]) %in% as.character(subFRcolNames)
                
                ZscoreCalc <- function(x, avg, sd){(x-avg)/sd}
                
                
                meanFRWOI <- lapply(seq(1, length(sessIdx)), function(i){
                        if(sum(masterDF[[c]]$sessfromCPsess==sessIdx[i], na.rm=T)>0){
                                
                                dataSel <- filter(masterDF[[c]], sessfromCPsess==sessIdx[i])
                                
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
                                                MeanByBin <- colMeans(numericDF, na.rm=T)
                                                MeanByUnit <- mean(MeanByBin, na.rm=T)
                                                MeanByUnitZsc <- ZscoreCalc(x=MeanByUnit, avg=BLaverage, sd=BLsd)
                                                CueExcited <- unitSel$CueExcited[1]
                                                return(data.frame(FRbyUnit=MeanByUnit, FRZsc=MeanByUnitZsc, CueExcited=CueExcited, UnitIdx=uniqUnits[u]))
                                        })
                                        )
                                        
                                        if(cueExcOnly==T){
                                                byUnit <- filter(byUnit, CueExcited==T)
                                        }
                                        
                                        
                                        if(dataProcess=="Zscores"){
                                                
                                                MeanByUnit <- byUnit$FRZsc
                                                        
                                                yAxMax=yAxMaxZ
                                                yAxMin=yAxMinZ
                                                labelLeg="(Z sc.)"
                                                
                                        } else {
                                                MeanByUnit <- byUnit$FRbyUnit
                                                
                                                yAxMax=yAxMaxRaw
                                                yAxMin=0
                                                labelLeg="(Hz)"
                                        }
                                        
                                        plot.window(xlim=c(0, length(sessIdx)+1), ylim=c(yAxMin, yAxMax+3))
                                        
                                        MeanByUnit <- MeanByUnit[!is.nan(MeanByUnit)]
                                        UnitIdx <- byUnit$UnitIdx[!is.nan(MeanByUnit)]
                                        
                                        #barSide <- (c-2)+(c-1) #This will put VEH side to the left and AP5 side to the right
                                        sessBarsWidth <- 0.25
                                        compBarBreaks <- seq(-sessBarsWidth, sessBarsWidth, length.out = length(masterDF)+1)
                                        leftRect <- compBarBreaks[c]
                                        rightRect <- compBarBreaks[c+1]
                                        
                                        Q1 <- summary(MeanByUnit)[2]
                                        Q3 <- summary(MeanByUnit)[5]
                                        IQR <- IQR(MeanByUnit)
                                        Median <- summary(MeanByUnit)[3]
                                        
                                        #IQR rectangle
                                        rect(xleft=i+leftRect, xright=i+rightRect, ybottom=Q1, ytop = Q3, col = color[c], border="white")
                                        
                                        #Mean (white) and Median (black) line
                                        segments(x0=i+leftRect, x1=i+rightRect, y0=mean(MeanByUnit), y1=mean(MeanByUnit), lwd=2, col = "black")
                                        segments(x0=i+leftRect, x1=i+rightRect, y0=Median, y1=Median, lwd=2, col="white")
                                        
                                        #Number of Units
                                        halfBar <- (compBarBreaks[c+1]-compBarBreaks[c])/2
                                        text(x=i+leftRect+halfBar, y=yAxMin, labels = length(uniqUnits), col=color[c], cex=0.7)
                                        
                                        if(morethanIQR==T){ 
                                                #Whiskers: maximum value still within Q3+1.5*IQR (whatever is smaller) or minimum value Q1-1.5*IQR
                                                overTop <- MeanByUnit>(Q3+1.5*IQR); top <- max(MeanByUnit[overTop==F])
                                                underBottom <- MeanByUnit<(Q1-1.5*IQR); bottom <- min(MeanByUnit[underBottom==F])
                                                topWhisker <- min(max(MeanByUnit), top)
                                                bottomwhisker <- max(min(MeanByUnit), bottom)
                                                
                                                segments(x0=i+leftRect+halfBar, x1=i+leftRect+halfBar, y0=Q3, y1=topWhisker)
                                                segments(x0=i+leftRect+halfBar, x1=i+leftRect+halfBar, y0=Q1, y1=bottomwhisker)
                                                
                                                overWhisker <- MeanByUnit[overTop]
                                                underWhisker <- MeanByUnit[underBottom]
                                                
                                                #Outliers
                                                points(x=rep(i+leftRect+halfBar, length(overWhisker)), y=overWhisker, cex=0.2, pch=19)
                                                points(x=rep(i+leftRect+halfBar, length(underWhisker)), y=underWhisker, cex=0.2, pch=19)
                                                
                                        }
                                        
                                        if(removeOutliers==T){
                                                outlierIdx <- (1:length(MeanByUnit))[(overTop==T | underBottom==T)]
                                                if(length(outlierIdx)>0){MeanByUnit <- MeanByUnit[-outlierIdx]}
                                        }
                                        
                                }
                                
                                return(list(MeanByUnit, UnitIdx)) 
                        }
                })
                
                if(dataProcess=="Zscores"){labelLeg="(Z sc.)"} else {labelLeg="(Hz)"}
                
                mtext(side=2, line=2.5, cex=1.5, font=2, text=paste("Firing rate", labelLeg, sep=" "))
                
                return(meanFRWOI)
                
        })
        
        axis(side=1, at=seq(1, length(sessIdx)), labels=(sessFromCP[1]:sessFromCP[2])+1, cex.axis=1.4, tick = F)
        mtext(side=1, line=2.5, cex=1.5, font=2, text="Session")
        
        #Add axis, labels and legend
        if(dataProcess=="Zscores"){yAxMax=yAxMaxZ; yAxMin=yAxMinZ}
        if(dataProcess=="raw"){yAxMax=yAxMaxRaw; yAxMin=yAxMinRaw}
        
        axis(side=2, at=seq(yAxMin, yAxMax, by=2), las=2, cex.axis=1.4, pos=0.6)

        #############################################################################
        ## Lines for individual neurons on S+ and S- trials
        
        if(lines==TRUE & sum(names(factors)=="Cue")>=1){
                
                #Identify what dataframes in the masterDF parameter are S+ or S- data
                DS_DF <- (1:length(comp))[comp %in% "S+"]
                NS_DF <- (1:length(comp))[comp %in% "S-"]
                
                minBoxArea <- 0.3
                maxBoxArea <- 0.3
                SubBoxArea <- (maxBoxArea+minBoxArea)/length(masterDF)
                halfBoxArea <- SubBoxArea/2
                xleftSeq <- seq(-minBoxArea, maxBoxArea-SubBoxArea, length.out = length(masterDF))
                
                sapply(seq(1, length(masterDF)/2), function(c){
                        
                        DS_DFidx <- DS_DF[c] 
                        NS_DFidx <- NS_DF[c]
                        
                        DStrials <- FRbyUnitBoth[[DS_DFidx]]
                        NStrials <- FRbyUnitBoth[[NS_DFidx]]
                        
                        xleft=xleftSeq[DS_DFidx]+halfBoxArea
                        xright=xleftSeq[NS_DFidx]+halfBoxArea
                        
                        nSess <- length(DStrials)
                        
                        #For each session
                        lapply(seq(1, nSess), function(i){
                                
                                DStrialsPerBin <- DStrials[[i]][[1]]
                                NStrialsPerBin <- NStrials[[i]][[1]]
                                
                                #For each unit
                                lapply(seq(1, length(DStrialsPerBin)), function(l){ 
                                        
                                        lines(x=c(i+xleft, i+xright), y=c(DStrialsPerBin[l], NStrialsPerBin[l]), col="black")
                                        
                                })
                        })
                        
                })
        }
                
        
        
      
        
        #############################################################################

        #Add info about analysis. 
        if(analysis==TRUE){
                missingData <- unlist(lapply(seq(1, length(FRbyUnitBoth)), function(d){
                        lapply(seq(1, length(FRbyUnitBoth[[d]])), function(s){
                                is.null(FRbyUnitBoth[[d]][[s]])})}))
                
                #If I have no data for some of the cells, then just print the graph with no analyses. Else continue w the analyses
                if(sum(missingData)>0){dev.off()} else {
                        longDF <- do.call("rbind", lapply(seq(1, length(FRbyUnitBoth)), function(d){ #For veh and drug or S+ and S-, whatever comparison I have
                                factorIdx <- expand.grid(factors)
                                do.call("rbind",lapply(seq(1, length(sessIdx)), function(s){ #For each session
                                        
                                        if(length(factors)==2){
                                                
                                                out <- data.frame(sess=sessIdx[s]+1, Factor1=factorIdx[d, 1], Factor2=factorIdx[d, 2], FR=FRbyUnitBoth[[d]][[s]][[1]], UnitIdx=FRbyUnitBoth[[d]][[s]][[2]])
                                                
                                                if(sum(names(factors)=="Modality")>=1 & sum(names(factors)=="Cue") >=1){
                                                        if(factorIdx[d, 1]=="Tone S+" & factorIdx[d, 2]=="S+"){k <- "Tone"}
                                                        if(factorIdx[d, 1]=="Tone S+" & factorIdx[d, 2]=="S-"){k <- "Light"}
                                                        if(factorIdx[d, 1]=="Light S+" & factorIdx[d, 2]=="S+"){k <- "Light"}
                                                        if(factorIdx[d, 1]=="Light S+" & factorIdx[d, 2]=="S-"){k <- "Tone"}
                                                        
                                                        out <- data.frame(sess=sessIdx[s]+1, Factor1=factorIdx[d, 1], Factor2=factorIdx[d, 2], sense=k, FR=FRbyUnitBoth[[d]][[s]][[1]], UnitIdx=FRbyUnitBoth[[d]][[s]][[2]])
                                                        
                                                }
                                                
                                        }
                                        
                                        if(length(factors)==1){
                                                out <- data.frame(sess=sessIdx[s]+1, Factor1=comp[d], FR=FRbyUnitBoth[[d]][[s]][[1]], UnitIdx=FRbyUnitBoth[[d]][[s]][[2]])
                                        }
                                        
                                        out
                                        
                                })
                                )
                        })
                        )
                }
                
                
                if(length(factors)==1){
                        
                        anovaTest <- ezANOVA(data=longDF, dv=FR, wid=UnitIdx, between=sess, within=Factor1, type=3)
                        
                        wilcTest <- do.call("rbind", lapply(seq(1, length(sessIdx)), function(s){
                                comp1 <- FRbyUnitBoth[[1]][[s]][[1]]
                                comp2 <- FRbyUnitBoth[[2]][[s]][[1]]
                                if(mean(comp1, na.rm = T) > mean(comp2, na.rm = T)){altpick="greater"}
                                if(mean(comp1, na.rm = T) < mean(comp2, na.rm = T)){altpick="less"}
                                if(mean(comp1, na.rm = T) == mean(comp2, na.rm = T)){altpick="two.sided"}
                                
                                TEST <- wilcox.test(comp1, comp2, paired=F, alt=altpick)
                                
                                data.frame(W=TEST$statistic, p.val=TEST$p.value)
                        })
                        )
                        
                        wilcTest$p.adj <- p.adjust(wilcTest$p.val, method="holm")
                        wilcTest$sig <- giveStars(wilcTest$p.adj)
                        
                        #Plot stars where diff is significant
                        sapply(seq(1, length(sessIdx)), function(s){
                                text(x=s, y=yAxMax+2, labels=wilcTest$sig[s], font=2, cex=1.5)
                        })
                }
                
                if(length(factors)==2){
                        
                        longDF$Factor1 <- factor(longDF$Factor1, levels = factors[[1]],labels =  factors[[1]])
                        longDF$Factor2 <- factor(longDF$Factor2, levels = factors[[2]],labels = factors[[2]])
                        
                        anovaTest <- ezANOVA(data=longDF, dv=FR, between=c(sess, Factor1), within=Factor2, wid=UnitIdx, type=3)
                        
                        
                        # wilcTest <- do.call("rbind", lapply(seq(1, length(sessIdx)), function(s){
                        #         comp1 <- FRbyUnitBoth[[1]][[s]]
                        #         comp2 <- FRbyUnitBoth[[2]][[s]]
                        #         if(mean(comp1, na.rm = T) > mean(comp2, na.rm = T)){altpick="greater"}
                        #         if(mean(comp1, na.rm = T) < mean(comp2, na.rm = T)){altpick="less"}
                        #         if(mean(comp1, na.rm = T) == mean(comp2, na.rm = T)){altpick="two.sided"}
                        #         
                        #         TEST <- wilcox.test(FRbyUnitBoth[[1]][[s]], FRbyUnitBoth[[2]][[s]], paired=F, alt=altpick)
                        #         
                        #         data.frame(W=TEST$statistic, p.val=TEST$p.value)
                        # })
                        # )
                        # 
                        # wilcTest$p.adj <- p.adjust(wilcTest$p.val, method="holm")
                        # wilcTest$sig <- giveStars(wilcTest$p.adj)
                        # 
                        # #Plot stars where diff is significant
                        # sapply(seq(1, length(sessIdx)), function(s){
                        #         text(x=s, y=yAxMax+2, labels=wilcTest$sig[s], font=2, cex=1.5)
                        # })
                }
                
                dev.off()
                
                save(longDF, file = paste(dataForRdir, "longDF.rdat", sep=""))
                if(length(factors)==1){testresults <- list(anovaTest, wilcTest)}
                if(length(factors)==2){testresults <- list(anovaTest)}
                
                save(testresults, file = paste(dataForRdir, "testresults.R", sep=""))  
                return(list(longDF, testresults))
                
                
        }
     
                dev.off()
               
        }
        
        


save(plotBoxplotfromSessCP, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/plotBoxplotfromSessCP.R")
save(plotBoxplotfromSessCP, file="E:/Dropbox/NMDA/R Functions/plotBoxplotfromSessCP.R")
save(plotBoxplotfromSessCP, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R Functions/plotBoxplotfromSessCP.R")
