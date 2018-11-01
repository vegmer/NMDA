

plotPSTHfromSessCP <- function(experiment="Exp 1", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                               graphFolder=MixedGraphFolder, dataProcess="Zscores", comp=c("Tone S+", "Light S+"),
                               correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                               yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                               imgFormat="pdf", neudata=allNeuronsDS){
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        if(!require(matrixStats)){install.packages("matrixStats")}
        library(matrixStats)
        
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        
        if(imgFormat=="pdf"){pdf(file=paste(graphFolder, experiment, "FR by sess from CP PSTH", "Sessions", sessFromCP, dataProcess, trialSel, ".pdf", sep="_"))}
        if(imgFormat=="png"){png(filename=paste(graphFolder, experiment, "FR by sess from CP PSTH", Sessions, sessFromCP, dataProcess, trialSel, ".png", sep="_"))}
        
        binw <- neudata$parameters$binw
        psthWins <- ((psthmin*1000)/binw):((psthmax*1000)/binw)
        psthWinMin <- (psthmin*1000)/binw
        psthWinMax <- (psthmax*1000)/binw
        
        sessIdx <- sessFromCP[1]:sessFromCP[2]
        
        plot.new()
        
        # This function has 2 functions: 
        # a) Calculate and spit the mean FR per bin around the time of the event for each session (w respect to change point) for each group of units (VEH vs AP5)
        # b) Plot that info
        
        FRbyBinBoth <- lapply(seq(1, length(masterDF)), function(c){
                
                
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                
                subFRcolNames <- (unique(masterDF[[c]]$CueBin)+psthWinMin):(unique(masterDF[[c]]$CueBin)+psthWinMax) #Select bins to be plotted
                subFRcols <- colnames(masterDF[[c]]) %in% as.character(subFRcolNames)
               
                ZscoreCalc <- function(x, avg, sd){(x-avg)/sd}
                
                FRbyBin <- lapply(seq(1, length(sessIdx)), function(i){
                        if(sum(masterDF[[c]]$sessfromCPsess==sessIdx[i], na.rm=T)>0){
                                
                                dataSel <- masterDF[[c]][masterDF[[c]]$sessfromCPsess==sessIdx[i], ]
                                
                                if(correctOnly==TRUE){dataSel <- dataSel[!is.na(dataSel$CueResponse), ]}
                                
                                if(cueExcOnly==TRUE){dataSel <- dataSel[dataSel$CueExcited==T, ]}
                                
                                if(sum(is.na(dataSel[,1]))!=nrow(dataSel)){ #If no rows are left after the filters I just applied, then ignore the following code. Only apply if there are units to apply it to
                                        numericDF <- apply(dataSel[, subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                                        BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                                        BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                                        
                                        if(dataProcess=="Zscores"){
                                                numericDFZsc <- sapply(seq(1, nrow(numericDF)), function(w){
                                                        sapply(seq(1, length(numericDF[w,])), function(v){
                                                                ZscoreCalc(x=numericDF[w,][v], avg=BLaverage[w], sd=BLsd[w])})
                                                })
                                                numericDFZsc <- t(numericDFZsc)
                                                
                                                MeanByBin <- colMeans(numericDFZsc, na.rm=T)
                                                SEMByBin <- colSds(numericDFZsc, na.rm=T)/sqrt(nrow(numericDFZsc))
                                                yAxMax=yAxMaxZ
                                                labelLeg=paste(yAxMax, "(Zsc.)")
                                                
                                        } else {
                                                MeanByBin <- colMeans(numericDF, na.rm=T)
                                                SEMByBin <- colSds(numericDFna.rm=T)/sqrt(nrow(numericDF))
                                                yAxMax=yAxMaxRaw
                                                labelLeg=paste(yAxMax, "(Hz)")
                                        }
                                        
                                        plot.window(xlim=c(1, length(subFRcolNames)), ylim=c(0, 3*(length(sessIdx)+1)*yAxMax))
                                        
                                        
                                        abline(h=(3*i*yAxMax)+0, col="gray60", lty=3)
                                        errCloud(col=color[c], x=1:ncol(numericDF), y=(3*i*yAxMax)+MeanByBin, err=SEMByBin)
                                        lines(x=1:ncol(numericDF), y=(3*i*yAxMax)+MeanByBin, col=color[c], lwd=1.5, type='l')
                                        nUnits <- length(unique(dataSel$allUnitIdx))
                                        text(x=10, y=(3*i*yAxMax)+c+(c-1)*2, col=color[c], labels = paste(comp[c], "n=", nUnits, sep=" "))
                                }
                                
                                return(MeanByBin) 
                        }
                })
                
                return(FRbyBin)
                
                
        })

        
        #Add axis, labels and legend
        if(dataProcess=="Zscores"){yAxMax=yAxMaxZ}
        if(dataProcess=="raw"){yAxMax=yAxMaxRaw}
        
        axis(side=2, las=2, at=seq(3*1*yAxMax, 3*length(sessIdx)*yAxMax, length.out = length(sessIdx)), labels=sessIdx, cex.axis=1.4)
        mtext(side=2, line=2.5, text = "Session from CP session", font=2, cex=1.5)
        
        #Redefine subFRcolNames bc it was inside the function
        subFRcolNames <- (unique(masterDF[[1]]$CueBin)+psthWinMin):(unique(masterDF[[1]]$CueBin)+psthWinMax) #Select bins to be plotted
        
        labs <- seq(psthmin, psthmax, by=binw/100)
        steps <- seq(1, length(subFRcolNames), by=10)
        axis(side=1, at=steps, labels=labs, cex.axis=1.4)
        mtext(side=1, "Time from S+ onset (s)", cex=1.5, font=2, line=2.5)
        abline(v=seq(1, length(subFRcolNames))[seq(psthmin, psthmax, by=binw/1000)==0], lwd=1.5)
        #abline(h=40, lty=2)
        
        #FR legend
        if(dataProcess=="Zscores"){yAxMax=yAxMaxZ; labelLeg=paste(yAxMax, "(Zsc.)")} 
        if(dataProcess=="raw"){yAxMax=yAxMaxRaw; labelLeg=paste(yAxMax, "(Zsc.)")}
        
        segments(x0=length(subFRcolNames)-2, x1=length(subFRcolNames), y0=length(sessIdx)*3*yAxMax, y1=length(sessIdx)*3*yAxMax, lwd=1.5)
        segments(x0=length(subFRcolNames), x1=length(subFRcolNames), y0=length(sessIdx)*3*yAxMax, y1=(length(sessIdx)*3*yAxMax)+3*yAxMax, lwd=1.5)
        text(labelLeg, x=length(subFRcolNames)+0.5, y=(length(sessIdx)*3*yAxMax)+1.5*yAxMax, srt=90, cex=0.8)
        
        dev.off()
        
        #############################################################################
        
        #Add info about analysis. 
        #I'm going to compare, for each session from CP, FR in both kinds of units using Wilcoxon for unpaired samples (Wilcoxon rank sum test)
        
        # #Window of interest for my analysis
        # WOIstart=000 #in ms
        # WOIend=400 #in ms
        # 
        # #Index of the values that contain that info in FRbyBin
        # WOIcolNames <- (unique(masterDF[[c]]$CueBin)+WOIstart/binw):(unique(masterDF[[c]]$CueBin)+WOIend/binw) #Select bins to be plotted
        # WOIcols <- subFRcolNames %in% as.character(WOIcolNames)
        # 
        # #
        # Group1data <- sapply(seq(1, length(FRbyBinBoth[[1]])), function(i){FRbyBinBoth[[1]][[i]][WOIcols]})
        # Group2data <- sapply(seq(1, length(FRbyBinBoth[[2]])), function(i){FRbyBinBoth[[2]][[i]][WOIcols]})
        # 
        # wilcoxResults <- do.call("rbind", lapply(seq(1, length(FRbyBinBoth[[1]])), function(i){
        #         meanDiff <- mean(Group1data[,i], na.rm = T)-mean(Group2data[, i], na.rm = T)
        #         if(meanDiff>0){altpick="greater"}
        #         if(meanDiff<0){altpick="less"}
        #         if(meanDiff==0){altpick="two.sided"}
        #         
        #         test <- wilcox.test(Group1data[,i], Group2data[,i], paired = F, alternative = altpick)
        #         
        #         return(data.frame(SessFromCP=sessIdx[i], mean1=mean(Group1data[,i], na.rm=T),
        #                           mean2=mean(Group2data[,i], na.rm=T), altpick=altpick,
        #                           W=test$statistic, p.val=test$p.value))
        #         
        # })
        # )
        # 
        # wilcoxResults$p.adj <- p.adjust(wilcoxResults$p.val, method="holm")
        # 
        # ### PLOT A GRAPH WITH THE NEURONAL DATA PER SESSION IN BIG BINS
        

}


save(plotPSTHfromSessCP, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/plotPSTHfromSessCP.R")
save(plotPSTHfromSessCP, file="E:/Dropbox/NMDA/R Functions/plotPSTHfromSessCP.R")
save(plotPSTHfromSessCP, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R Functions/plotPSTHfromSessCP.R")
