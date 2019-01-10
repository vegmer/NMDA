CPextract <- function(GallCrit=1.3, minSlope=0.5, adjDrop=-0, idx, plot=T, CPfuncFolder, CPGraphFolder, dataForRdir, dataForRCumulative){
        
        # Install and call necessary packages
        if(!require(R.utils)){install.packages("R.utils")}
        library(R.utils)
        
        # Define folder with scripts and test data
        Rscripts <- CPfuncFolder
        
        # Compile all functions in the folder Rscripts into the environment
        sourceDirectory(Rscripts, modifiedOnly=FALSE) #Now all the necessary functions are loaded into our environment
        
        #### Set important folders: 
        cumsumDataFolder <- dataForRCumulative
        graphFolder <- CPGraphFolder
        
        # All behavioral data objects for that group
        files <- paste(dataForRdir, list.files(dataForRdir), sep="")
        filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep="")
        for(i in 1:length(files)){load(files[[i]])}
        for(i in 1:length(filesCum)){load(filesCum[[i]])}
        
        #Apply CP wrapper
        Lindex <- list.files(cumsumDataFolder)
        index_names <- sapply(seq(1, length(Lindex)), function(m){
                dotpos <- gregexpr(".rdat", Lindex[m])[[1]][1] #Find position of '.rdat'
                substring(Lindex[m], 1, dotpos-1)})
        
        #This has 2 functions: make graphs with each rat's cumulative performance and CP analysis and save those graphs in the CP folder and create an object summarizing that ([[i]]-->rat; [[j]]--> index name)
        CPrawdata <- lapply(seq(1, length(rats)), function(i){
                lapply(seq(1, length(index_names)), function(j){
                        rat=rats[i]
                        index_name=index_names[j]
                        
                        if(index_name=="DSrespAll"){indexSel <- DSrespAll; idxlab='S+ response ratio'}
                        if(index_name=="DStaskAcc"){indexSel <- DStaskAcc; idxlab='S+ specificity'}
                        if(index_name=="DStimeToSpare"){indexSel <- DStimeToSpare; idxlab='S+ "time to spare"'}
                        if(index_name=="NSrespAll"){indexSel <- NSrespAll; idxlab='S- response ratio'}
                        if(index_name=="DSlatency"){indexSel <- DSlatency; idxlab='S+ latency'}
                        if(index_name=="NSlatency"){indexSel <- DSlatency; idxlab='S- latency'}
                        if(index_name=="ITIlatency"){indexSel <- DSlatency; idxlab='ITI latency'}
                        if(index_name=="NStaskAcc"){indexSel <- NStaskAcc; idxlab='S- specificity'}
                        if(index_name=="NStimeToSpare"){indexSel <- NStimeToSpare; idxlab='S- "time to spare"'}
                        if(index_name=="ITIrespRatio"){indexSel <- ITIrespRatio; idxlab="ITI response ratio"}
                        
                        dataset <- as.matrix(indexSel[[i]])
                        
                        if(is.na(match(index_name, list.files(CPGraphFolder)))){dir.create(paste(CPGraphFolder, index_name, sep=""))}
                        
                        #This line creates and saves PDF graphs and generates object with change point data
                        cp_wrapper2(dataset, isDiscrete=1, test=2, Crit=GallCrit, ratname=rat, index=c(index_name, idxlab), graphFolder=graphFolder)
                        
                })
        })
        
        
        save(CPrawdata, file=paste(dataForRdir, "CPrawdata.rdat"))
        
        ###########################################################################
        # Based on the results of Gallistel, I'm going to select the earliest change point that is followed by only "positive" change points (more than 0 in the cue specificity measure)
        #
        
        DSspecif <- lapply(seq(1, length(rats)), function(x){
                dat <- CPrawdata[[x]][[3]]
                upward <- diff(dat$slope)>0 #Is the change positive? 1st value applies to 2nd row of CP data, so add a FALSE to the beg. Also, switch the last one by a TRUE (bc the last one is the last trial it's always going to say slope=0)
                fixupward <- c(FALSE, upward[-length(upward)], TRUE)
                dat$upward <- fixupward
                
                #Are CPs followed by a slope of at least 1s (because I realized that sometimes I have change points followed by 0.5s slope or whatever, it doesn't look like there's learning going on but the algorithm marked it)
                positive <- dat$slope>minSlope
                fixpositive <- c(positive[-length(positive)], TRUE)
                dat$positive <- fixpositive
                
                #There's a few cases where the animal was doing really well and then for a few trials the slope goes slightly under 0. I don't want to miss a CP because of that, so I'll set a criterion: if the previous slope was more than 3 and the current slope is more than -0.5, count that as positive
                positiveAdj <- dat$positive
                negativeIdx <- (1:nrow(dat))[dat$positive==FALSE]
                #The first negativeIdx value is always 1, and it doesn't have a 'previous slope' value. So when I select [negativeIdx-1], it disappears. Assign a 2, so that when I select [negativeIdx-1], I get the slope of 1, which is never going to be over 3
                negativeIdx[(negativeIdx-1)==0] <- 2
                AdjIdx <- dat$slope[negativeIdx-1] > 3 & dat$slope[negativeIdx]>adjDrop
                if(sum(AdjIdx)>0){positiveAdj[negativeIdx[AdjIdx]] <- TRUE}
                dat$positiveAdj <- positiveAdj
                
                #dat$bothcrit <- dat$upward==T & dat$positiveAdj==T
                dat$realCP <- FALSE
                
                dat$nSessions <- length(idx[[x]])
                
                CPs <- dat$CP_trial
                
                #Previously responded to S+ trials
                CPs <- CPs+1 #The algorithm calls trial "1" trial "0". So every trial is really 1 trial more in my other objects
                #prevRespCSplus <- PreRespondedCSplus[[x]][CPs]
                #dat$prevRespCSplus <- prevRespCSplus
                
                #CTratio: sum up until that trial
                #CTratioSum <- sapply(seq(1, length(CPs)), function(n){sum(CTratio[[x]][0:CPs[n]])})
                #dat$CTratioSum <- CTratioSum
                
                trends <- rle(dat$positiveAdj)
                last <- function(w){w[length(w)]}
                
                if(last(trends$values)==TRUE & last(trends$lengths)>1){ #If there is a string of CPs
                        CPpos <- length(dat$realCP)-(last(trends$lengths)-1)
                        dat$realCP[CPpos] <- TRUE
                }
                
                #There are a few cases where there's only 2 values (1st and last trial). It's possible that the CP is trial 1 if the animal starts off by doing really well. But otherwise it's not a real CP. How to distinguish?
                #Only accept this as a CP if the slope is at least 2.
                if(nrow(dat)==2 & dat$slope[1]<2){dat$realCP=FALSE}
                
                return(dat)
        })
        
        save(DSspecif, file=paste(dataForRdir, 'DSspecif.rdat', sep=""))
        
        
        ##########################################################################
        ### PLOT
        
        if(plot==T){
                sapply(seq(1, length(DSspecif)), function(x){
                        
                        #Info for saving graph
                        ratname=rats[x]
                        lastTrial <- length(DStaskAcc[[x]])
                        filetitle = paste(graphFolder, "/", paste(ratname, "DStaskAcc"), ".pdf", sep="")
                        pdf(file=filetitle)
                        par(mfrow=c(3, 1))
                        
                        #Plot DS specificity per trial/per bin
                        plot(DStaskAcc[[x]], type="l", lwd=1.5, xlab='Trial', ylab="S+ specificity")
                        abline(h=seq(-7.5, 7.5, by=2.5), lwd=0.5, col="gray70")
                        title(paste(ratname, "S+ specificity"))
                        
                        abline(h=0, lwd=1.5)
                        
                        
                        #Draw blue line with avg S+ specificity by bin
                        binsize <- 5
                        nbins <- length(DStaskAcc[[x]])/binsize
                        allbins <- 1:nbins
                        trialidx <- 1:lastTrial
                        binidx <- c(); for(i in 1:nbins){binidx <- c(binidx, rep(allbins[i], binsize))}
                        meanPerBin <- sapply(seq(1, nbins), function(y){
                                sel <- binidx==allbins[y]
                                mean(DStaskAcc[[x]][sel], na.rm=T)
                        })
                        meanToPlot <- c(); for(i in 1:nbins){meanToPlot <- c(meanToPlot, rep(meanPerBin[i], binsize))}
                        lines(meanToPlot, lwd=2, col="blue")
                        
                        #Draw lines demarcating the sessions
                        nTrialSess <- sapply(seq(1, length(idx[[x]])), function(k){
                                sessIdx <- idx[[x]][k]
                                length(alldata[[sessIdx]]$CSpluscue)
                        })
                        
                        abline(v=cumsum(nTrialSess), lwd=0.5, col="gray20")
                        
                        #Plot cumulative performance
                        plot(cumsum(DStaskAcc[[x]]), type="l", lwd=1.5, xlab='Trial', ylab="Cumulative S+ spec.")
                        abline(h=0, lwd=1.5)
                        abline(v=c(0, cumsum(nTrialSess)), lwd=0.5, col="gray50")
                        CP <- DSspecif[[x]]$CP_trial[DSspecif[[x]]$realCP]
                        points(x=DSspecif[[x]]$CP_trial[-1], y=cumsum(DStaskAcc[[x]])[DSspecif[[x]]$CP_trial])
                        if(sum(DSspecif[[x]]$realCP)>0){
                                points(x=CP, y=DSspecif[[x]]$cum_DStaskAcc[DSspecif[[x]]$realCP], col="red")
                                abline(v=CP, col="red", lwd=2)}
                        
                        
                        #Plot mean pre and post change point accroding to Gallistel's method
                        plot.new()
                        plot.window(xlim=c(1, lastTrial), ylim=c(-7, 7))
                        
                        lines(DStaskAcc[[x]])
                        
                        #MEAN PRE AND POST CP by Gallistel's method
                        if(sum(DSspecif[[x]]$realCP)==0){
                                mean <- mean(DStaskAcc[[x]], na.rm=T)
                                sd <- sd(DStaskAcc[[x]], na.rm=T)
                                
                                segments(x0=1, x1=lastTrial, y0=mean, y1=mean, lwd=2)
                                
                                abline(h=0, lwd=1.5)
                                
                                axis(side = 1, at=seq(0, lastTrial, by=10), labels = seq(0, lastTrial, by=10))
                                axis(side = 2, at=seq(-7, 7, 1))
                                
                                mtext("Trial", side=1, line=2.5, cex=1.5)
                                mtext("S+ specificity", side=2, line=2.5, cex=1.5)
                                
                        } else {
                                ###Calculate mean and sd pre and post change point
                                meanPre <- mean(DStaskAcc[[x]][1:CP], na.rm=T)
                                sdPre <- sd(DStaskAcc[[x]][1:CP], na.rm=T)
                                meanPost <- mean(DStaskAcc[[x]][(CP+1):lastTrial], na.rm=T)
                                sdPost <- sd(DStaskAcc[[x]][(CP+1):lastTrial], na.rm=T)
                                
                                segments(x0=1, x1=CP, y0=meanPre, y1=meanPre, lwd=2, col="red")
                                segments(x0=CP, x1=lastTrial, y0=meanPost, y1=meanPost, lwd=2, col="red")
                                segments(x0=CP, x1=CP, y0=meanPre, y1=meanPost, lwd=2, col="red")
                                
                                
                                abline(h=0, lty=2)
                                
                                axis(side = 1, at=seq(0, lastTrial, by=10), labels = seq(0, lastTrial, by=10))
                                axis(side = 2, at=seq(-7, 7, 1))
                                
                                mtext("Trial", side=1, line=2.8, cex=1.5)
                                mtext("S+ specificity", side=2, line=2.5, cex=1.5)
                                
                        }
                        
                        dev.off()
                }) 
        }
        
        #############################################3
        
        CPsummary <- sapply(seq(1, length(DSspecif)), function(x){
                CP <- DSspecif[[x]]$CP_trial[DSspecif[[x]]$realCP]
                if(length(CP)!=0){
                        slopePre <- mean(DStaskAcc[[x]][1:CP], na.rm=T)
                        slopePost <- mean(DStaskAcc[[x]][(CP+1):length(DStaskAcc[[x]])], na.rm=T)
                        
                        #In what session did the CP take place?
                        nSess <- length(idx[[x]])
                        trialsPerSess <- sapply(seq(1, nSess), function(y){length(alldata[[idx[[x]][y]]]$CSpluscue)})
                        trialDistr <- c(0, cumsum(trialsPerSess))
                        CPsess <- findInterval(CP, trialDistr)
                        PreCPRewTrials <- DSspecif[[x]]$prevRespCSplus[DSspecif[[x]]$realCP] #Number of rewarded trials before CP
                        PercCPRewTrials <- PreCPRewTrials/CP #Percentage of rewarded trials before CP
                } else {
                        nSess <- length(idx[[x]])
                        CP <- NA
                        CPsess <- NA
                        slopePre <- mean(DStaskAcc[[x]], na.rm=T)
                        slopePost <- slopePre
                        PreCPRewTrials <- NA
                        PercCPRewTrials <- NA
                }
                
                perc80Asympt <- 0.8*slopePost
                #First CP with a slope over 80% of the asymptote (mean slope post CP if there actually was a CP)
                if(!is.na(CP)){
                        CPOver80pc <- DSspecif[[x]]$CP_trial[DSspecif[[x]]$slope>=perc80Asympt]
                        FirstOver80pc <- min(CPOver80pc[CPOver80pc>=DSspecif[[x]]$CP_trial[DSspecif[[x]]$realCP]])
                        DynamicInterv <- FirstOver80pc-DSspecif[[x]]$CP_trial[DSspecif[[x]]$realCP]
                } else {FirstOver80pc=NA; DynamicInterv <- NA} 
                
                
                
                results <- c(CP, CPsess, slopePre, slopePost, FirstOver80pc, DynamicInterv, nSess)
                return(results)
        })
        
        CPdata <- data.frame(rat=rats, CP=CPsummary[1,], CPsess=CPsummary[2,], slopePre=CPsummary[3,], slopePost=CPsummary[4,], FirstCPOver80pc=CPsummary[5,], DynamicInterv=CPsummary[6,], nSess=CPsummary[7,])
        
        NS.PreCP <- sapply(seq(1, length(rats)), function(k){
                allcues <- unlist(lapply(seq(1, length(idx[[k]])), function(l){
                      alldata[[idx[[k]][l]]]$orderCues 
                }))
                
                CPrat <- CPdata$CP[k]
                globalPosCP <- min((1:length(allcues))[cumsum(allcues==1)==CPrat]) #Position of DS CP in the global order of cues
                numNSPreCP <- max(cumsum(allcues[1:globalPosCP]==2))
                
        })
        
        CPdata$NS.PreCP <- NS.PreCP
        
        save(CPdata, file=paste(dataForRdir, "CPdata.rdat", sep=""))
        
        return(CPdata)
}

save(CPextract, file="E:/Dropbox/NMDA/R functions/CPextract.R")
