#Possible values for behIndex:
# - "S+ specificity"
# - "S+ response ratio"
# - "S+ latency"
# - "ITI latency"

#Values for "datapoint":
# - datapoint="unitsMean". In the plot, one data point per trial, aggregating units from same and different animals.
# - datapoint="individual units". In the plot, one data point per unit recorded per trial


bySessionPlots <- function(data=masterDF, neudataParameters=allNeuronsDS$parameters, 
                          winmin=100, winmax=300, behIndex="S+ specificity", datapoint="unitsMean",
                          setsOfSess <- list(c(-5:-1), c(1:4))){
        
        sessIdx <- sort(unique(data$sessfromCPsess))
        
        if(!is.na(setsOfSess){sessIdx}
        
        ## Here's where I select the columns of interest for calculating mean FR in the period defined by winmin winmax
        CueBin <- unique(selSess$CueBin) #Bin right after cue onset
        binw <- neudataParameters$binw #bin width for FR in ms
        
        minBin <- CueBin + winmin/binw #Bin that corresponds to the beginning of the period of interest (defined by winmin)
        maxBin <- CueBin + winmax/binw #Bin that corresponds to the end of the period of interest (defined by winmax)
        WOI <- minBin:maxBin
        
        #Window of interest firing rate columns
        WOI_FRcols <- colnames(masterDF)==WOI
        
        #Find the index of the column that contains the performance index defined in the function parameter "behIndex".
        if(behIndex=="S+ specificity"){PerfCol <- colnames(data)=="CSplusSpecif"}
        if(behIndex=="S+ response ratio"){PerfCol <- colnames(data)=="CSplusResponse"}
        if(behIndex=="S+ latency"){PerfCol <- colnames(data)=="CSplusLat"}
        if(behIndex=="ITI latency"){PerfCol <- colnames(data)=="ITIlatency"}
        
        #Now, I'm going to make a plot with performance on one axis and firing rate in the period of interest (winmin:winmax) on the other axis
        
        sapply(seq(1, length(sessIdx)), function(x){
                
                selSess <- data[data$sessfromCPsess==sessIdx[x], ]
                nUnit <- length(unique(selSess$allUnitIdx))
                nrats <- length(unique(selSess$rat))
                
                FRraw <- selSess[,WOI_FRcols]
                
                #For some reason, R treats FRraw as a list of ncol vectors and it won't let me treat it as a data frame. Make numeric.
                #It also treats BLavg in a weird way (as an integer instead of numeric, I lose the decimals if I try to convert to numeric, so I have to do this:)
                FRraw <- apply(FRraw, MARGIN=2, as.numeric)
                BLavg <- as.numeric(format(selSess$BLavg, digits=2))
                BLsd <- as.numeric(format(selSess$BLsd, digits=2))
                
                FR <- ZscoreCalc(x=rowMeans(FRraw, na.rm=T), avg=BLavg, sd=BLsd)
                Perf <- selSess[, PerfCol]
                if(behIndex=="S+ response ratio"){
                        Perf[!is.na(Perf)] <- 1
                        Perf[is.na(Perf)] <- 0
                        }
                
                #Determine axes range
                xmin <- floor(min(Perf))
                xmax <- ceiling(max(Perf))
                ymin <- floor(min(FR))
                ymax <- ceiling(max(FR))
                
                #Make a graph with 2 columns and the necessary number of rows to fit them
                graphRows <- ceiling(length(sessIdx/2))
                
                if(datapoint=="unitsMean"){
                        uniqTrial <- unique(selSess$uniqTrial)
                        plot.new()
                        plot.window(xlim=c(xmin, xmax), ylim=c(ymin, 15))
                        sapply(seq(1, length(uniqTrial)), function(y){
                                rowSel <- selSess$uniqTrial==uniqTrial[y]  #Rows that contribute to data from one single trial
                                points(x=mean(Perf[rowSel]) , y=mean(FR[rowSel]), pch=19, cex=1)
                                axis(side=1, cex.axis=1.4)
                                mtext(side=1, line=2.5, text=behIndex, cex=1.5, font=2)
                                axis(side=2, cex.axis=1.4, las=2)
                                mtext(side=2, line=2.5, text="Firing rate (Z sc.)", cex=1.5, font=2)
                                abline(lm(mean( mean(FR)[rowSel] ~ Perf[rowSel])), col="red")
                                
                        })
                
                }
                
                if(datapoint=="individual units"){
                        abline(lm(FR ~ Perf), col="red")
                        plot(x=Perf, y=FR)  
                }
                
        
                
        })
}