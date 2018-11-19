cumulative_CP <- function(Exp="Exp3", CPdata=CPdata, numSess=6, byAnimal=TRUE, byNeuron=FALSE, 
                          nexdata=allNeuronsDS, graphFolder=MixedGraphFolder){
        
        
        if(byAnimal==TRUE){
                cumCP <- sapply(seq(1, numSess), function(x){
                        sum(CPdata$CPsess<=x)/numSess
                })
                
                byWhat <- "By Animal"
        }
        
        if(byNeuron==TRUE){
                
                unitsPerSess <- do.call("rbind", lapply(seq(1, length(nexdata)), function(x){
                        ratname <- nexdata[[x]]$ratname
                        expt <- nexdata[[x]]$expt
                        nUnits <- sum(grepl("sig", names(nexdata[[x]])))
                        
                        #Some rats have a 0 before their number in nexdata but not in CPdata, fix that. Also the ones with 3 numbers, drop the first one
                        if(ratname=="MV06"){ratname <- "MV6"}
                        if(ratname=="MV05"){ratname <- "MV5"}
                        if(ratname=="MV08"){ratname <- "MV8"}
                        if(ratname=="MV86"){ratname <- "MV186"}
                        if(ratname=="MV90"){ratname <- "MV190"}
                        
                        
                        CP <- CPdata[CPdata$rat==ratname, ]$CPsess
                        
                        data.frame(ratname, expt, CPsess=CP, nUnits)
                })
                )
                
                totalUnits <- sum(unitsPerSess[unitsPerSess$CPsess <= nSess, ]$nUnits, na.rm=T)
                
                cumCP <- sapply(seq(1, numSess), function(x){
                        
                        sum(unitsPerSess[unitsPerSess$CPsess <= x, ]$nUnits, na.rm=T)/totalUnits
                        
                })
                
                byWhat <- "By Unit"
        }
        
        
        
        filename <- paste(graphFolder, Exp, " Cumulative CP by session ", byWhat, ".pdf", sep="")
        
        pdf(file = filename)
        
        par(oma=c(2,2,2,2))
        plot.new()
        plot.window(xlim=c(1, numSess), ylim=c(0, 1))
        
        lines(x=seq(1, numSess), y=cumCP, lwd=2)
        #points(x=seq(1, numSess), y=cumCP, pch=19)
        
        axis(side=1, cex.axis=1.4, font=2)
        axis(side=4, cex.axis=1.4, font=2, las=2)
        
        dev.off()
}


save(cumulative_CP, file="E:/Dropbox/NMDA/R Functions/cumulative_CP.R")
