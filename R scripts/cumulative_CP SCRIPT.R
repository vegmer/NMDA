cumulative_CP <- function(Exp="Exp3", CPdata=CPdata, numSess=6, byAnimal=TRUE, byNeuron=FALSE, 
                          neudata=allNeuronsDS, graphFolder=MixedGraphFolder){
        
        
        if(byAnimal==TRUE){
                cumCP <- sapply(seq(1, numSess), function(x){
                        sum(CPdata$CPsess<=x)/numSess
                })
                
                byWhat <- "By Animal"
        }
        
        if(byNeuron==TRUE){
                
                nexdata <- neudata$nexdata
                
                unitsPerSess <- do.call("rbind", lapply(seq(1, length(nexdata)), function(x){
                        ratname <- nexdata[[x]]$ratname
                        expt <- nexdata[[x]]$expt
                        nUnits <- sum(grepl("sig", names(nexdata[[x]])))
                        CP <- CPdata[CPdata$rat==ratname, ]$CPsess
                        
                        data.frame(ratname, expt, CPsess=CP, nUnits)
                })
                )
                
                totalUnits <- sum(unitsPerSess[unitsPerSess$CPsess <= nSess, ]$nUnits)
                
                cumCP <- sapply(seq(1, numSess), function(x){
                        
                        sum(unitsPerSess[unitsPerSess$CPsess <= x, ]$nUnits)/totalUnits
                        
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
