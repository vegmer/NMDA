

plotBoxplotPrePostCP <- function(experiment="Exp 4", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                  comp=c("VEH", "AP5"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                                  correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                                  factors=list(drug=c("VEH", "AP5")), ANOVA=FALSE,
                                  yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, WdwStart=0, WdwEnd=400, removeCPsess=F,
                                  removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T){
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        if(!require(matrixStats)){install.packages("matrixStats")}
        library(matrixStats)
        if(!require(ez)){install.packages("ez")}
        library(ez)
        
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        
        if(imgFormat=="pdf"){pdf(file=paste(graphFolder, experiment, "FR Pre vs. Post CP Boxplot", "Sessions", sessFromCP[1], dataProcess, trialSel, ".pdf", sep="_"))}
        if(imgFormat=="png"){png(filename=paste(graphFolder, experiment, "FR Pre vs. Post CP Boxplot", Sessions, sessFromCP[1], dataProcess, trialSel, ".png", sep="_"))}
        
        binw <- neudata$parameters$binw
        minBin <- WdwStart/binw
        maxBin <- WdwEnd/binw
        
        selBins <- minBin:maxBin
        
        sessIdx <- sessFromCP[1]:sessFromCP[2]
        
        PrePostCP <- function(dat){
                new <- sapply(seq(1, length(dat)), function(d){
                        if(is.na(dat[d])){a <- NA} else {
                                if(dat[d]<0){a <- -1}
                                if(dat[d]==0){a <- 0}
                                if(dat[d]>0){a <- 1}
                        }
                        return(a)
                })
                return(new)
        }
        
        PrePostCPidx <- PrePostCP(sessIdx)
        
        ZscoreCalc <- function(x, avg, sd){(x-avg)/sd}
        
        
        plot.new()
        
        # This function has 2 functions: 
        # a) Calculate and spit out the mean FR per bin around the time of the event for each session (w respect to change point) for each group of units (VEH vs AP5)
        # b) Plot that info
        
        FRbyUnitBoth <- lapply(seq(1, length(masterDF)), function(c){
                
                
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                
                subFRcolNames <- (unique(masterDF[[c]]$CueBin)+minBin):(unique(masterDF[[c]]$CueBin)+maxBin) #Select bins to be plotted
                subFRcols <- colnames(masterDF[[c]]) %in% as.character(subFRcolNames)
                
                PrePostCPidx <- PrePostCP(masterDF[[c]]$sessfromCPsess)
                
                masterDF[[c]]$BeforeCP <- PrePostCPidx
                
                divs <- unique(PrePostCPidx)[!is.na(unique(PrePostCPidx))]
                
                meanFRWOI <- lapply(seq(1, length(divs)), function(i){
                        
                        
                        if(sum(masterDF[[c]]$BeforeCP==divs[i], na.rm=T)>0){
                                
                                dataSel <- filter(masterDF[[c]], BeforeCP==divs[i])
                                
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
                                                return(data.frame(FRbyUnit=MeanByUnit, FRZsc=MeanByUnitZsc, CueExcited=CueExcited))
                                        })
                                        )
                                        
                                        if(cueExcOnly==T){
                                                byUnit <- filter(byUnit, CueExcited==T)
                                        }
                                        
                                        
                                        if(dataProcess=="Zscores"){
                                                
                                                MeanByUnit <- byUnit$FRZsc
                                                yAxMin=yAxMinZ
                                                yAxMax=yAxMaxZ
                                                labelLeg="(Z sc.)"
                                                
                                        } else {
                                                MeanByUnit <- byUnit$FRbyUnit
                                                yAxMin=0
                                                yAxMax=yAxMaxRaw
                                                labelLeg="(Hz)"
                                        }
                                        
                                        
                                        
                                        plot.window(xlim=c(0, length(divs)+1), ylim=c(yAxMin, yAxMax+3))
                                        
                                        MeanByUnit <- MeanByUnit[!is.nan(MeanByUnit)]
                                        
                                        sessBarsWidth <- 0.25
                                        compBarBreaks <- seq(-sessBarsWidth, sessBarsWidth, length.out = length(masterDF)+1)
                                        leftRect <- compBarBreaks[c]
                                        rightRect <- compBarBreaks[c+1]
                                        
                                        #barSide <- (c-2)+(c-1) #This will put VEH side to the left and AP5 side to the right
                                        
                                        Q1 <- summary(MeanByUnit)[2]
                                        Q3 <- summary(MeanByUnit)[5]
                                        IQR <- IQR(MeanByUnit)
                                        Median <- summary(MeanByUnit)[3]
                                        
                                        #IQR rectangle
                                        rect(xleft=i+leftRect, xright=i+rightRect, ybottom=Q1, ytop = Q3, col = color[c], border="white")
                                        
                                        #Mean (BLACK) and Median (WHITE) line
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
                                
                                return(MeanByUnit) 
                        }
                })
                
                mtext(side=2, line=2.5, cex=1.5, font=2, text=paste("Firing rate", labelLeg, sep=" "))
                
                
                return(meanFRWOI)
                
                
        })
        
        
        divs <- unique(PrePostCP(masterDF[[1]]$sessfromCPsess))
        divs <- divs[!is.na(divs)]
        
        if(length(factors)==1){names(FRbyUnitBoth) <- factors[[1]]}
                
      
        axis(side=1, at=seq(1, length(divs)), labels=c("Before CP", "CP session", "After CP"), cex.axis=1.4, tick = F)
        #mtext(side=1, line=2.5, cex=1.5, font=2, text="Session from CP session")
        
        #Add axis, labels and legend
        if(dataProcess=="Zscores"){yAxMax=yAxMaxZ; yAxMin=yAxMinZ}
        if(dataProcess=="raw"){yAxMax=yAxMaxRaw; yAxMin=yAxMinRaw}
        
        axis(side=2, at=seq(yAxMin, yAxMax, by=2), las=2, cex.axis=1.4, pos=0.6)
        
        
        #############################################################################
        
        #Add info about analysis. 
        
        longDF <- do.call("rbind", lapply(seq(1, length(FRbyUnitBoth)), function(d){ #For veh and drug
                factorIdx <- expand.grid(factors)
                do.call("rbind",lapply(seq(1, length(FRbyUnitBoth[[d]])), function(s){ #For each session
                        if(length(factors)==1){
                                pr <- data.frame(sess=divs[s], comp=comp[d], FR=FRbyUnitBoth[[d]][[s]])
                        }
                        
                        if(length(factors)==2){
                                pr <- data.frame(sess=divs[s], mod=factorIdx[d, 1], down=factorIdx[d, 2], FR=FRbyUnitBoth[[d]][[s]])
                        }   
                        
                        pr
                        
                })
                )
        })
        )
        
        
        if(removeCPsess==T){ #If I want to discard data from the CP session (good idea bc the behavior fluctuates a lot in this session and also I have many less units per condition)
                FRbyUnitBoth <- lapply(seq(1, length(FRbyUnitBoth)), function(x){
                        FRbyUnitBoth[[x]][c(divs!=0)]
                })
                
                longDF <- longDF[longDF$sess !=0, ]
        }
        
        longDF$UnitIdx<- 1:nrow(longDF)
        
       
        
        #Convert my factors to factors 
        if(length(factors)==1){
                longDF$sess <- factor(longDF$sess, levels = divs,labels = as.character(divs))
                longDF$comp <- factor(longDF$comp, levels = comp,labels = comp)
                longDF$UnitIdx <- factor(longDF$UnitIdx, levels=1:nrow(longDF), as.character(1:nrow(longDF)))
                
                if(ANOVA==TRUE){anovaTest <- ezANOVA(data=longDF, dv=FR, between=c(sess, comp), wid=UnitIdx, type=3)}
                
                
                #### POST-HOC COMP WITH WILCOXON RANK SUM TESTS
                sess_lab <- c("Pre CP", "Post CP")
                factor_lab <- c(factors[[1]])
                
                lab_matrix <- sapply(seq(1, length(sess_lab)), function(q){
                        paste(sess_lab[[q]], factor_lab)
                })
                
                #Comparison, within each session (or session mashup), of the effect of the factor in question
                wilcTest1 <- do.call("rbind", lapply(seq(1, length(FRbyUnitBoth)), function(s){
                        comp1 <- FRbyUnitBoth[[1]][[s]] #Down pre (s=1) or post (s=2) vs. Not Down pre(s=1) or post(s=2)/ Or VEH pre (s=1) or post (s=2) vs. AP5 pre or post (s=1 or s=2)
                        comp2 <- FRbyUnitBoth[[2]][[s]] 
                        if(mean(comp1, na.rm = T) > mean(comp2, na.rm = T)){altpick="greater"}
                        if(mean(comp1, na.rm = T) < mean(comp2, na.rm = T)){altpick="less"}
                        if(mean(comp1, na.rm = T) == mean(comp2, na.rm = T)){altpick="two.sided"}
                        
                        TEST <- wilcox.test(FRbyUnitBoth[[1]][[s]], FRbyUnitBoth[[2]][[s]], paired=F, alt=altpick)
                        
                        compName <- paste(lab_matrix[1, s], "vs.", lab_matrix[2, s])
                        
                        data.frame(Comparison=compName, W=TEST$statistic, p.val=TEST$p.value)
                })
                )
                
                #Comparison, within each factor condition (e.g. Down or Not Down), of the effect of training (before CP vs. after CP)
                #Similar for VEH vs. AP5
                wilcTest2 <- do.call("rbind", lapply(seq(1, length(unique(longDF$sess))), function(s){
                        
                        comp1 <- FRbyUnitBoth[[s]][[1]] ##Down pre (s=1) or Not Down pre (s=2) vs. Down post (s=1) or Not Down post (s=2)
                        comp2 <- FRbyUnitBoth[[s]][[2]]
                        
                        if(mean(comp1, na.rm = T) > mean(comp2, na.rm = T)){altpick="greater"}
                        if(mean(comp1, na.rm = T) < mean(comp2, na.rm = T)){altpick="less"}
                        if(mean(comp1, na.rm = T) == mean(comp2, na.rm = T)){altpick="two.sided"}
                        
                        TEST <- wilcox.test(comp1, comp2, paired=F, alt=altpick)
                        
                        compName <- paste(lab_matrix[s, 1], "vs.", lab_matrix[s, 2])
                        
                        data.frame(Comparison=compName, W=TEST$statistic, p.val=TEST$p.value)
                })
                )
                
                wilcTest <- rbind(wilcTest1, wilcTest2)
                
                wilcTest$p.adj <- p.adjust(wilcTest$p.val, method="holm")
                wilcTest$sig <- giveStars(wilcTest$p.adj)
                
                #Plot stars where diff is significant
                sapply(seq(1, length(sess_lab)), function(s){
                        if(s==1){xpos=1}
                        if(s==2){xpos=3} #Bc I excluded CP session from the analysis
                        text(x=xpos, y=yAxMax+2, labels=wilcTest$sig[s], font=2, cex=1.5)
                })
                
                 
        }
        
        if(length(factors)==2){
                longDF$sess <- factor(longDF$sess, levels = divs,labels = as.character(divs))
                longDF$mod <- factor(longDF$mod, levels = factors[[1]],labels =  factors[[1]])
                longDF$down <- factor(longDF$down, levels = factors[[2]],labels = factors[[2]])
                longDF$UnitIdx <- factor(longDF$UnitIdx, levels=1:nrow(longDF), as.character(1:nrow(longDF)))
                
                if(ANOVA==TRUE){anovaTest <- ezANOVA(data=longDF, dv=FR, between=c(sess, comp), wid=UnitIdx, type=3)}
                
                
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
        
        legend("topleft", col = color, legend=comp, pch=15)
        
        dev.off()
        
        if(ANOVA==TRUE){return(list(anovaTest, wilcTest))} else {
                return(wilcTest)
        }
        
        
}


save(plotBoxplotPrePostCP, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/plotBoxplotPrePostCP.R")
save(plotBoxplotPrePostCP, file="E:/Dropbox/NMDA/R Functions/plotBoxplotPrePostCP.R")
save(plotBoxplotPrePostCP, file="E:/Dropbox/NMDA/EXP4_Unilateral AP5/R Functions/plotBoxplotPrePostCPP.R")
