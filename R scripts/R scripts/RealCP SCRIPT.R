RealCP <- function(data=CPdataAllCrit, method="intersect", index_names=c("DSlatency", "DSrespAll", "DStaskAcc", "DStimeToSpare", "ITIlatency", "NSlatency", "NSrespAll", "NStaskAcc", "NStimeToSpare"), rats=rats){
        #In the object "CPdataAllCrit", these are the levels:
        # [[x]] --> Rat
        # [[y]] --> Index
        # [[z]] --> Signifiance criterion in Gallistel's algorithm
        
        CPbyIndex <- sapply(seq(1, length(data)), function(x){
                sapply(seq(1, length(data[[x]])), function(y){
                        allCrit <- data[[x]][[y]]
                        CPbyCrit <-  sapply(seq(1, length(allCrit)), function(z){
                                allCrit[[z]]$CP_trial
                        })
                        
                        if(method=="mostconservative"){
                                mostconsIdx <- length(CPbyCrit)
                                candidates <- CPbyCrit[[mostconsIdx]]
                                
                                #For some reason, the first and last trials are always considered CPs, get rid of them.
                                candidates <- candidates[-c(1, length(candidates))]
                                
                                if(length(candidates)!=0){CP <- min(candidates)} else {CP <- NA}
                        }
                        
                        if(method=="intersect"){
                                #This methods finds the FIRST trial that appears as a CP in every level of the criterion
                                commonCPs <- Reduce(intersect, CPbyCrit)
                                
                                #For some reason, the first and last trials are always considered CPs, get rid of them.
                                commonCPs <- commonCPs[-c(1, length(commonCPs))] 
                                
                                if(length(commonCPs)!=0){CP <- min(commonCPs)} else {CP <- NA}
                        }
                        
                       
                        
                })
        })
        
        colnames(CPbyIndex) <- rats
        rownames(CPbyIndex) <- index_names
        
        return(CPbyIndex)
        
}   




save(RealCP, file="E:/Dropbox/NMDA/R functions/RealCP.R")
save(RealCP, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R functions/RealCP.R")
