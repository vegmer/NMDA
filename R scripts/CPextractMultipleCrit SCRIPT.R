CPextractMultipleCrit <- function(GallCrit=c(1.3, 2, 4, 5, 6), idx=idx, CPfuncFolder=CPfuncFolder, CPGraphFolder=CPGraphFolder, dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative){
        
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
        CPrawdata <- lapply(seq(1, length(rats)), function(i){                  #Per rat
                lapply(seq(1, length(index_names)), function(j){                #Per index
        
                        rat=rats[i]
                        index_name=index_names[j]
                        
                        if(index_name=="DSrespAll"){indexSel <- DSrespAll; idxlab='S+ response ratio'}
                        if(index_name=="DStaskAcc"){indexSel <- DStaskAcc; idxlab='S+ specificity'}
                        if(index_name=="DStimeToSpare"){indexSel <- DStimeToSpare; idxlab='S+ "time to spare"'}
                        if(index_name=="NSrespAll"){indexSel <- NSrespAll; idxlab='S- response ratio'}
                        if(index_name=="DSlatency"){indexSel <- DSlatency; idxlab='S+ latency'}
                        if(index_name=="NSlatency"){indexSel <- DSlatency; idxlab='S- latency'}
                        if(index_name=="ITIlatency"){indexSel <- ITIlatency; idxlab='ITI latency'}
                        if(index_name=="NStaskAcc"){indexSel <- NStaskAcc; idxlab='S- specificity'}
                        if(index_name=="NStimeToSpare"){indexSel <- NStimeToSpare; idxlab='S- "time to spare"'}
                        if(index_name=="ITIrespRatio"){indexSel <- ITIrespRatio; idxlab="ITI response ratio"}
                        
                        dataset <- as.matrix(indexSel[[i]]) #Data has to be in a nx1 matrix format
                        
                        if(is.na(match(index_name, list.files(CPGraphFolder)))){dir.create(paste(CPGraphFolder, index_name, sep=""))}
                        
                        lapply(seq(1, length(GallCrit)), function(k){               #Per value of Gallistel criterion
                                if(is.na(match(GallCrit[k], list.files(paste(CPGraphFolder, index_name, sep=""))))){dir.create(paste(CPGraphFolder, index_name, "/", GallCrit[k], sep=""))}
                        
                                #This line creates and saves PDF graphs and generates object with change point data. I made this version for when I input more than one GallCrit value
                                cp_wrapper3(dataset, isDiscrete=1, test=2, Crit=GallCrit[k], ratname=rat, index=c(index_name, idxlab), graphFolder=graphFolder)
                                
                                })
                        
                
                       
                })
        })
        
}





save(CPextractMultipleCrit, file="E:/Dropbox/NMDA/R functions/CPextractMultipleCrit.R")