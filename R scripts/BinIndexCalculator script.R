#Function to identify bincuts for each event (DS, NS and all cues together). It creates one object for each and inside each object, one vector per rat

BinIndexCalculator <- function(data=alldata, binsize=300, sessLength=9000){
        binCuts <- seq(0, sessLength, by=binsize)
        DSbinIdx <- lapply(seq(1, length(alldata)), function(x){
                DSstart <- alldata[[x]]$CSpluscue
                DSbinIdx <- findInterval(DSstart, binCuts)})
        
        NSbinIdx <- lapply(seq(1, length(alldata)), function(x){
                NSstart <- alldata[[x]]$CSminuscue
                AllCueStart <- alldata[[x]]$allCues
                NSbinIdx <- findInterval(NSstart, binCuts)
        })
        
        AllCueBinIdx <- lapply(seq(1, length(alldata)), function(x){
                AllCueStart <- alldata[[x]]$allCues
                AllCueBinIdx <- findInterval(AllCueStart, binCuts)
        })
        
        save(DSbinIdx, file=paste(dataForRCumulative, "DSbinIdx.rdat", sep=""))
        save(NSbinIdx, file=paste(dataForRCumulative, "NSbinIdx.rdat", sep=""))
        save(AllCueBinIdx, file=paste(dataForRCumulative, "AllCueBinIdx.rdat", sep=""))
}


funcdirect <- "E:/Dropbox/NMDA/R functions/"
save(BinIndexCalculator, file=paste(funcdirect, "BinIndexCalculator.R", sep=""))
