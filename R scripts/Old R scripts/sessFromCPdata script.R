
sessFromCPdata <- function(data=masterDF_DS, sessFromCP=0, FRparameters=allNeuronsDS$parameters, 
                           folder=BySessFolder, winmin=0, winmax=400, BLmin=-2000, BLmax=0){
        
       
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
        
        lapply(seq(1, length(sessFromCP)), function(d){
                
                sess <- sessFromCP[d]
                
                daySel <- data[data$sessfromCPsess==sess, ]
                
                dayUnits <- unique(daySel$allUnitIdx) #All of the units recorded on the day defined by sessFromCP
                
        
                datForRegr <- do.call("rbind", lapply(seq(1, length(dayUnits)), function(x){
                        
                        unitIdx <- dayUnits[x]
                        selUnit <- daySel[daySel$allUnitIdx==unitIdx, ]
                        FR_WOI <- selUnit[ ,WOI_FRcols] #Raw firing rate on the bins of interest (Window Of Interest)
                        FR_WBl <- selUnit[ ,BL_win_FRcols] #Raw firing rate on the bins defined as baseline
                        
                        
                        #For some reason FR_WOI is interpreted by R as a list of lists, I need to convert it into a matrix in this roundabout way.
                        FR_WOI <- apply(FR_WOI, MARGIN=2, as.numeric)
                        FRavg <- rowMeans(FR_WOI)
                        
                        
                        FR_WBl <- apply(FR_WBl, MARGIN=2, as.numeric)
                        BLavg <- as.numeric(format(selUnit$BLavg, digits=2))
                        BLsd <- as.numeric(format(selUnit$BLsd, digits=2))
                        
                        
                        FR_WOI_Zsc <- ZscoreCalc(x=FRavg, avg=BLavg, sd=BLsd)
                        trialfromCPidx <- selUnit$trialfromCP
                        
                        data.frame(rat=selUnit$rat, session=selUnit$session, sessFromCP=sessFromCP[d], 
                              allUnitIdx=dayUnits[x], TrialFromCP=trialfromCPidx, TrialNumber=selUnit$trialIdx, 
                              Fr_Zsc=FR_WOI_Zsc, Fr_Raw=FRavg, BLavg=BLavg, BLsd=BLsd, CueExcited=selUnit$CueExcited,
                              CueResponse=selUnit$CueResponse, CueLatency=selUnit$CueLatency,
                              CueSpecif=selUnit$CueSpecif, ITIlat=selUnit$ITIlatency
                        )
                }))
                
                #datForRegr <- as.data.frame(datForRegr)
                
                #cbind converted my logical vector to a factor vector of 1 (FALSE) and 2s (TRUE). So fix that:
                #datForRegr$CueExcited <- datForRegr$CueExcited==2
                
                datForRegr
                
        })
        
}


save(sessFromCPdata, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/sessFromCPdata.R")
save(sessFromCPdata, file="E:/Dropbox/NMDA/R Functions/sessFromCPdata.R")
