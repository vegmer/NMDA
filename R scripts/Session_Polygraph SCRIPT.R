#################################################################################################################################
# FUNCTION TO CREATE A RECORD, TRIAL BY TRIAL, OF EACH NEURON'S CUE-EVOKED FIRING RATE ON THE FIRST SESSION (OR OTHER SESSIONS) #
#################################################################################################################################


Session_Polygraph <- function(data=Day1_TrialUnit, cues=c("S+", "S-"), day=1, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder){
        
        
        ### BEHAVIOR INFO
        ratsWithUnits <- as.character(unique(data[[1]]$rat))
        sessIdx <- idx[c(rats %in% ratsWithUnits)]
        sess_Day_idx <- sapply(sessIdx, function(r){r[day]}) #Select the first day of each rat and find the index in the global list of files
        
        
        cuesorder <- lapply(1:length(sessIdx), function(r){
                alldata[[sess_Day_idx[[r]]]]$orderCues
        })
        
        if(("S+" %in% cues)==TRUE){
                DSresponse <- lapply(1:length(sess_Day_idx), function(r){
                        !is.na(alldata[[sess_Day_idx[[r]]]]$CSplusresponse)})
                if(length(cues)==1){
                        cuesorder <- lapply(1:length(sessIdx), function(r){
                                rep(1, length(unique(data[[1]]$bin)))})
                        
                        response <-  list(DSresp=DSresponse)
                }
                
        }
        
        if(("S-" %in% cues)==TRUE){
                NSresponse <- lapply(1:length(sessIdx), function(r){
                        !is.na(alldata[[sess_Day_idx[[r]]]]$CSminusresponse)})
                if(length(cues)==1){
                        cuesorder <- lapply(1:length(sessIdx), function(r){
                                rep(2, length(unique(data[[1]]$bin)))})
                        
                        response <-  list(NSresp=NSresponse)
                }
        }
        
        if(length(cues)==2){
                response <- list(DSresp=DSresponse, NSresp=NSresponse)
        }
        
        filename <- paste(plotFolder, "By unit by trial Day ", day, " on ", paste(cues, collapse = ""), " trials.pdf", sep="")
        
        pdf(file = filename, onefile = TRUE)
        
        lapply(seq(1, length(ratsWithUnits)), function(r){
                
                plot.new() 
                
                dataSel <- data[[1]][data[[1]]$rat==ratsWithUnits[r], ]
                
                nTrials <- length(cuesorder[[r]])
                nUnits <- length(unique(dataSel$unitIdx))
                
                plot.window(xlim=c(1, nTrials), ylim=c(0, (0.5+nUnits)*(FRmax-FRmin)))
                
                if(("S+" %in% cues)==TRUE){ #If I have S+ trials in my selection
                        pick <- (1:length(cues))[cues %in% "S+"]
                        DSdataSel <- data[[pick]][data[[pick]]$rat==ratsWithUnits[r], ]
                        
                        DStrialPos <- (1:length(cuesorder[[r]]))[cuesorder[[r]]==1]
                        
                        sapply(1:length(DStrialPos), function(t){
                                rect(xleft = DStrialPos[t]-0.5, xright=DStrialPos[t]+0.5,
                                     ybottom = 0, ytop=(0.5+nUnits)*(FRmax-FRmin), col="gray90", border = "white")
                                
                                entered <- response$DSresp[[r]][t]
                                if(entered==TRUE){
                                        symbols(x=DStrialPos[t], y=((0.5+nUnits)*(FRmax-FRmin))+2, 
                                                inches=F, circles = 0.4, bg="black", add=T)
                                }
                        })
                        
                }
                
                if(("S-" %in% cues)==TRUE){ #If I have S- trials in my selection
                        
                        pick <- (1:length(cues))[cues %in% "S-"]
                        NSdataSel <- data[[pick]][data[[pick]]$rat==ratsWithUnits[r], ]
                        NStrialPos <- (1:length(cuesorder[[r]]))[cuesorder[[r]]==2]
                        
                        sapply(1:length(NStrialPos), function(t){
                                
                                entered <- response$NSresp[[r]][t]
                                if(entered==TRUE){
                                        symbols(x=NStrialPos[t], y=((0.5+nUnits)*(FRmax-FRmin))+2, inches=F, circles = 0.4, col="red", add=T)
                                }
                        })
                }
                
                
                if(length(cues)==1){
                        FRdat <- lapply(1:nUnits, function(u){
                                dataSel[dataSel$unitIdx==unique(dataSel$unitIdx)[u], ]$byUnitFR
                        })    
                }
                
                if(length(cues)==2){
                        
                        #Get the Firing rate data vector per unit on that rat
                        FRdat <- lapply(1:nUnits, function(u){
                                DS_FRdat <- DSdataSel[DSdataSel$unitIdx==unique(dataSel$unitIdx)[u], ]$byUnitFR
                                NS_FRdat <- NSdataSel[NSdataSel$unitIdx==unique(dataSel$unitIdx)[u], ]$byUnitFR
                                v <- list(DS_FRdat, NS_FRdat)
                                
                                #This interleaves the FR of the unit on S+ trials and S- trials based on the order of the trials during teh session (given by "cuesorder")
                                sapply(seq(1, length(cuesorder[[r]])), function(t){
                                        trialKind <- cuesorder[[r]][t]
                                        trialKindCum <- sum(cuesorder[[r]][1:t]==trialKind)
                                        FR <- v[[trialKind]][trialKindCum]
                                })
                                
                        })
                }
                
                
                lapply(seq(1, nUnits), function(u){
                        
                        ypos <- ((FRmax-FRmin)/2)+(u-1)*(FRmax-FRmin)
                        
                        #Line indicating 0 Z scores
                        abline(h=ypos, lty=2)
                        
                        # +2, 4, and 6 Zscores
                        abline(h=ypos-2); abline(h=ypos+2); abline(h=ypos+4); abline(h=ypos+6)
                        
                        #Actual FR in Z scores
                        lines(x=1:length(cuesorder[[r]]), y=ypos+FRdat[[u]], lwd=2)
                        
                        axis(side=2, at=c(ypos-2, ypos, ypos+2, ypos+4, ypos+6), 
                             labels=c(-2, 0, 2, 4, 6), las=2, cex.axis=0.5)
                        
                })
                
                axis(side=1)
                mtext(side=1, line=2.5, text = "Trial", font=2, cex=1.4)
                mtext(side=2, line=2.5, text = "Firing rate (Zsc)", font=2, cex=1.4)
                legend(x=0, y=1, fill=c("gray90", "white"), legend=c("S+", "S-"), ncol=2, cex=0.5)
                legend(x=nTrials-10, y=1, pt.bg=c("black", "white"), pch=21, legend=c("S+ entered", "S- entered"), ncol=2, cex=0.5)
                
                title(main=paste(ratsWithUnits[r], "day ", day, "FR per unit per trial", sep=" "))
                
        })
        
        
        dev.off()  
        
}


save(Session_Polygraph, file=paste(funcdirect, "Session_Polygraph.R", sep=""))
