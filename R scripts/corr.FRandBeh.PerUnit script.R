corr.FRandBeh.PerUnit <- function(experiment= "Exp 3 pre and post CP", 
                                  masterDF=masterDF_DS, dataProcess="Zscores",
                                  WdwStart=100, WdwEnd=400, neudata=allNeuronsDS, 
                                  fromCP=TRUE, bySess=FALSE, sessRange = c(-5, 5), collapsePrePost=FALSE, 
                                  byTrial=TRUE, capped=TRUE, capValue=c(-90, 90), trialBinSize=15,
                                  cueExcOnly=FALSE, correctOnly=FALSE, 
                                  graphFolder=MixedGraphFolder){
        
        ############################################################################
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        
        
        
        #############################################################################
        ## For creating the PDF's name
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        if(cueExcOnly==TRUE){cueExcSel="Cue Exc Units"} else {cueExcSel="all units"}
        
        filename=paste(graphFolder, experiment, "Hist of corr coef FR x Beh", WdwEnd, cueExcSel, trialSel, ".pdf", sep=" ")
        
        pdf(file=filename)
        

        ############################################################################
        ### CREATE AN INDEX TO BREAK DOWN THE TRIALS/SESSIONS IN THE WAY DEFINED BY THE PARAMETERS AND CREATE A DATA OBJECT WITH THE TRIALS WITHIN THAT RANGE
        
        if(length(sessRange)==2){sessIdx <- sessRange[1]:sessRange[2]}
        if(length(sessRange)==1){sessIdx <- sessRange}
        
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
        
        masterDF$PrePostCPsess <- PrePostCP(masterDF$sessfromCPsess)
        
        assignBINindx <- function(a){
                name <- 1:length(unique(a))
                sapply(seq(1, length(a)), function(x){name[a[x]==unique(a)]})
        }
        
        #If I want to break down data by session
        if(bySess==TRUE){
                
                #If the sessions are indexed with respect to the CP session
                if(fromCP==TRUE){
                        data <- masterDF[masterDF$sessfromCPsess %in% sessIdx, ]
                        data$BIN <- assignBINindx(data$sessfromCPsess)
                        
                        #If I want to collapse all the sessions before and all teh sessions after CP
                        if(collapsePrePost==TRUE){
                                data$BIN <- assignBINindx(data$PrePostCPsess)
                        }
                }
                
                #If the sessions are indexed from the first session in order
                if(fromCP==FALSE){
                        data <- masterDF[masterDF$session %in% sessIdx, ]
                        data$BIN <- assignBINindx(data$session)
                        
                }
              
        }
        
        #If I want to break down the data by trial
        if(byTrial==TRUE){
        
             #By trial counting from change point
             if(fromCP==TRUE){
                     
                     #Within a specified range (specified by the parameter "capValue")
                     if(capped==T){
                             nDivisions <- round((capValue[2]-capValue[1])/trialBinSize, 0)
                             trialBins <- seq(capValue[1], capValue[2], by=trialBinSize)
                             correspBins <- findInterval(masterDF$trialfromCP, c(min(trialBins)-trialBinSize, trialBins, max(trialBins)+trialBinSize)) #I add one bin on each side bc I don't want the bins in the extremes to include trials that were outside the trialbinsize range under or over that bin
                             
                             #Assign NA to the bins that fall on the "minimum bin-trialBinSize" or "maximum bin+trialBinSize". I added them to the calculation to get rid of trials that fall outside that range. Now I have to assign NA to those bins to make sure they don't get selected later on.
                             correspBins[which(correspBins==0)] <- NA
                             correspBins[which(correspBins==max(correspBins, na.rm=T))] <- NA
                             
                             masterDF$correspBins <- correspBins
                             
                             data <- masterDF[!is.na(masterDF$correspBins), ]
                             data$BIN <- data$correspBins
                             
                     }
                     
                     #All trials from CP
                     if(capped==F){
                             nDivisions <- round(((max(masterDF$trialfromCP, na.rm=T)-min(masterDF$trialfromCP, na.rm=T))/trialBinSize), 0)
                             trialBins <- seq(min(masterDF$trialfromCP, na.rm=T), max(masterDF$trialfromCP, na.rm=T), by=trialBinSize)
                             correspBins <- findInterval(masterDF$trialfromCP, trialBins)
                             
                             masterDF$correspBins <- correspBins
                             
                             data <- masterDF[!is.na(masterDF$correspBins), ]
                             data$BIN <- data$correspBins
                     } 
             }
           
           #If I just want all the data of every trial by trial         
           if(fromCP==FALSE){
                   nDivisions <- round(((max(masterDF$trialCum, na.rm=T)-min(masterDF$trialCum, na.rm=T))/trialBinSize), 0)
                   trialBins <-  seq(min(masterDF$trialCum, na.rm=T), max(masterDF$trialCum, na.rm=T), by=trialBinSize)
                   correspBins <- findInterval(masterDF$trialCum, trialBins)
                   masterDF$correspBins <- correspBins
                   data <- masterDF
                   data$BIN <- data$correspBins
           }  
                
        }
        
        
        
        
        ############################################################
        ### POST-CUE BINS OF INTEREST
        
        #Columns with the Firing Rate information
        FRcols <- (1:ncol(data))[is.element(colnames(data), 1:ncol(data))]
        
        #Bins within the window of interest
        binw <- neudata$parameters$binw
        CueBin <- data$CueBin[1]
        TargetStartBin <- CueBin + (WdwStart/binw)
        TargetEndBin <- CueBin + (WdwEnd/binw)
        selBins <- TargetStartBin:TargetEndBin
        
        
        ###########################################################
        ### BASELINE BINS
        BLStartBin <- CueBin - (BLWdw/binw)
        BLEndBin <- CueBin-1
        
        ############################################################
        ### FILTER UNITS AND TRIALS OF INTEREST
        if(cueExcOnly==TRUE){
                data <- data[data$CueExcited==TRUE, ]
        }
        
        
        if(correctOnly==TRUE){
                data <- data[!is.na(data$CueResponse), ]
        }
        
        
        ###############################
        #Create a data frame with FR and BEH info per neuron in the division of interest (group of trials or session)
        
        correlData <- lapply(unique(data$BIN), function(i){
                
                dataSel <- data[data$BIN==i, ]
                uniqUnits <- unique(dataSel$allUnitIdx)
                
                do.call("rbind", lapply(seq(1, length(uniqUnits)), function(u){
                        
                        unitSel <- dataSel[dataSel$allUnitIdx==uniqUnits[u], ]
                        
                        FR <- apply(unitSel[, FRcols], MARGIN=2, as.numeric)
                        
                        # BL_FR <- mean(FR[ ,BLStartBin:BLEndBin], na.rm=T)
                        # BL_FR_SD <- sd(rowMeans(FR[ ,BLStartBin:BLEndBin], na.rm = T)) 
                        # 
                        # 
                         BL_FR <- as.numeric(as.character(unitSel$BLavg[1]))
                         BL_FR_SD <- as.numeric(as.character(unitSel$BLsd[1]))
                        # 
                        CUE_FR <- rowMeans(FR[ ,TargetStartBin:TargetEndBin], na.rm=T)
                        CUE_FR_TotalSpikes <- (rowSums(FR[ ,TargetStartBin:TargetEndBin], na.rm=T))*50/1000 #Convert back to spikes
                        FR_Zsc <- (CUE_FR-BL_FR)/BL_FR_SD #Firing rate per trial in Z scores
                        
                        if(dataProcess=="Zscores"){firing=FR_Zsc}
                        if(dataProcess=="Freq"){firing=CUE_FR}
                        if(dataProcess=="Spikes"){firing=CUE_FR_TotalSpikes}
                        
                        BEH <- unitSel$CueLat #Behavior by trial
                        
                        CORREL <- cor.test(x=BEH, y=firing, method="spearman")
                        
                        rho <- CORREL$estimate
                        pval <- CORREL$p.value
                        
                        return(data.frame(rho=rho, pval=pval, cueExcited=unitSel$CueExcited[1],
                                          unit=uniqUnits[u], BIN=i, session=unitSel$session[1], 
                                          sessfromCPsess=unitSel$sessfromCPsess[1], 
                                          beforeCP=unitSel$PrePostCPsess[1]
                                          ))
                        
                })
                )
                
        })
        
        
        ######################################################################
        ### PLOT HISTOGRAMS WITH THIS
        lapply(seq(1, length(correlData)), function(i){
               
                CueExc.rho <- correlData[[i]][!is.na(correlData[[i]]$rho) & correlData[[i]]$cueExcited==T, ]$rho
                CueExc.rho.sig <- correlData[[i]][!is.na(correlData[[i]]$rho) & correlData[[i]]$cueExcited==T & correlData[[i]]$pval<=0.05, ]$rho
                
                #CueExc.rho <- correlData[[i]][!is.na(correlData[[i]]$rho), ]$rho
                #CueExc.rho.sig <- correlData[[i]][!is.na(correlData[[i]]$rho) & correlData[[i]]$pval<=0.05, ]$rho

                
                hist(CueExc.rho, breaks = seq(-1, 1, by=0.1), col="gray80", 
                     xlab="Correlation coefficient (Spearman)", ylab="Number of Units", 
                     main="Distribution of corr. coeff. Post-cue FR x Latency across units")
                hist(CueExc.rho.sig, breaks = seq(-1, 1, by=0.1), col="gray30", add=T)
                
                arrows(x0=median(CueExc.rho), x1=median(CueExc.rho), y0=7, y1=5, lwd=2, col="black")
                
                legend("topright", legend=c("Significant", "Not significant"), col=c("gray30", "gray80"), lty=1, lwd=2)
                legend("topleft", legend=paste("BIN=", correlData[[i]]$BIN[1], "; n units=", length(CueExc.rho)))
                
        })
        
        dev.off()
        
}
