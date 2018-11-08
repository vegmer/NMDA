neuralhist <- function(funcdirect=funcdirect, path, startt=0, endt, binw, psthmin, psthmax, event, cueexonly=F, side="both", allResults=F) { 
        
        # Load necessary functions
        functionFolder = funcdirect
        load(file = paste(functionFolder, "eventfinder.r", sep=""))
        load(file = paste(functionFolder, "neuraldataextractMV.r", sep=""))
        
        #specify path of nex files
        nexfiles = list.files(path)
        
        #extracts all NeuroExplorer data
        nexdata = lapply(nexfiles, function(x) {
                neuraldataextractMV(inpath=paste(path, x, sep = ""), funcdirect=funcdirect)})
        
        
        # The first one will select neurons in either drug-treated side, saline-treated side or both sides (based on "side" parameter)
        neurons = lapply(seq(1, length(nexdata)), function(x, SIDE=side) {
                
                drugSide = nexdata[[x]]$drug #This gives me the side of the brain in which the drug was infused (given by the name of the NEX file)
                ratname = nexdata[[x]]$ratname
                nexnames = names(nexdata[[x]])
                neuronidx = grep("sig", nexnames)
                
                #### SELECTION OF HEMISPHERE BASED ON "SIDE" PARAMETER ("drug": AP5-treated side; "vehicle":veh treated side; "both":both sides)
                
                neuronumbers <- as.numeric(substr(nexnames[neuronidx], 4, 6)) #Channel of each recorded unit
                
                #"code" will tell me on what side of the brain the recorded units were
                code <- sapply(seq(1, length(neuronumbers)), function(i){
                        if(neuronumbers[i]<=8){code="L"} else {code="R"} #Channels 1-8 are LEFT hemisphere and 9-16 are RIGHT hemisphere in rats with left connector in the front
                        
                        if( ratname=="MV90" | ratname=="MV05"){
                                if(neuronumbers[i]<=8){code="R"} else {code="L"} #Rats 90 and 05 had the RIGHT array in front, so the sides are actually opposite of what I just defined for the rest of the rats
                        }
                        
                        return(code)
                })
                
                if(SIDE=="drug"){
                        neuronidx <- neuronidx[code==drugSide]
                }
                
                if(SIDE=="vehicle"){
                        neuronidx <- neuronidx[code!=drugSide]
                }
                
                if(SIDE=="both"){
                        neuronidx <- neuronidx
                }
                
                
                ####
                if(length(neuronidx)==0){return(NA)} #If no neurons on that side were recorded, return NA
                if(length(neuronidx)>0) {return(nexdata[[x]][neuronidx])}
                        
        })
        
        
        
        events = eventfinder(nexdata, event)
        
        #######################################################################
        #specify the interval of timestamps you want to examine  
        sstart = startt
        send = endt
        
        #these two functions simply restrict the timestamp values
        #to the window specified above
        #######################################################
        
        neurons = lapply(seq(1, length(neurons)), function(x) {
                lapply(seq(1, length(neurons[[x]])), function(y) {
                        
                        neurons[[x]][[y]][which(neurons[[x]][[y]] >= sstart & neurons[[x]][[y]] <= send)]    
                        
                })
        })
        
        events = lapply(seq(1, length(events)), function(x) {
                events[[x]][which(events[[x]] >= sstart & events[[x]] <= send)] })
        
        
        
        #########################################################
        #Now that the relevant data is extracted, we can begin to 
        #analyze it by constructing rasters, histograms, etc.
        #########################################################
        binsize = binw #in ms
        winmin = psthmin #in s
        winmax = psthmax #in s
        
        
        
        #####################################################################
        #This function returns a list of lists of lists:
        #Level 1 is the experiment; Level 2 is each neuron in the experiment;
        #Level 3 is the time relative to the specified event for each trial
        #The object masterlist can then be used to construct PSTHs or rasters
        ########################################################################
        masterlist = lapply(seq(1, length(neurons)), function(x) {
                lapply(seq(1, length(neurons[[x]])), function(y) {
                        lapply(seq(1, length(events[[x]])), function(z) {
                                stampsidx = which(neurons[[x]][[y]] >= events[[x]][z] - winmin & neurons[[x]][[y]] <= events[[x]][z] + winmax)
                                relstamps = neurons[[x]][[y]][stampsidx] - events[[x]][z]
                        })
                })    
        })
        
        #########################################################################
        #This function will flag cue-excited neurons if that is the only population   
        #you wish to examine. Criterion: 2 consecutive 50 ms bins in the 500ms window after the cue
        #in which the FR exceeds the 99.9% conf. interval
        #for a Poisson distribution using 2 s pre-cue as baseline.
        #########################################################################
        cueexidx = lapply(seq(1, length(masterlist)), function(x) {
                lapply(seq(1, length(masterlist[[x]])), function(y) {
                        pbin = .05
                        threshold = 3
                        BLwdw=5
                        
                        allvals = unlist(masterlist[[x]][[y]])
                        hcounts = hist(allvals[which(allvals >= -BLwdw & allvals <= .5)], breaks = seq(-BLwdw, .5, pbin), plot = F)$counts 
                        
                        freq = (hcounts/pbin)/length(masterlist[[x]][[y]])  #just for comparison's sake, to see how counts compare to actual frequency
                        baseline = mean(freq[1:(BLwdw/pbin)]) 
                        
                        baselinecounts=hcounts[1:(BLwdw/pbin)]
                        
                        critwin = hcounts[(BLwdw/pbin + 1):((BLwdw+.5)/pbin)]         #this is the window (500 ms) in which we look for excitation
                        critval = poisson.test(x = sum(baselinecounts), T = length(baselinecounts), conf.level = 0.999)$conf.int[2]   #computes the upper limit of the confindence interval based on the baseline
                        
                        diffs = diff(which(critwin > critval))       #computes the differences in indices for bins exceeding the critical value (to check whether excited bins are consecutive or not)
                        cueex = F
                        if(length(which(rle(diffs)$values == 1 & rle(diffs)$lengths >= threshold))>0) (cueex = T)   #looks for consecutive bins (diff equal to 1) that are at least 2 bins long -or whatever #of bins specified by the argument "threshold" (rle length of at least 2)
                        
                        return(cueex)  
                })
        })
        
        
        propcueex <- lapply(seq(1, length(cueexidx)), function(x){
                length(which(cueexidx[[x]]==T))/length(cueexidx[[x]])
        })
        
        
        #########################################################################
        #This function will flag cue-inhibited neurons if that is the only population   
        #you wish to examine. Criterion (Kevin's paper, 2018): 2 consec 50 ms bins in the 
        #500ms post-cue window in which the FR exceeds the lowerlimit of the 99.9% 
        #confidence interval for a Poisson distribution using 2 s pre-cue as baseline.
        #########################################################################
        cueinhidx = lapply(seq(1, length(masterlist)), function(x) {
                lapply(seq(1, length(masterlist[[x]])), function(y) {
                        pbin = .05
                        threshold = 2
                        BLwdw=5
                        
                        allvals = unlist(masterlist[[x]][[y]])
                        hcounts = hist(allvals[which(allvals >= -BLwdw & allvals <= .5)], breaks = seq(-BLwdw, .5, pbin), plot = F)$counts 
                        
                        freq = (hcounts/pbin)/length(masterlist[[x]][[y]])  #just for comparison's sake, to see how counts compare to actual frequency
                        baseline = mean(freq[1:(BLwdw/pbin)]) 
                        baselinecounts=hcounts[1:(BLwdw/pbin)]
                        
                        critwin = hcounts[(BLwdw/pbin + 1):((BLwdw+.5)/pbin)]         #this is the window (500 ms) in which we look for excitation
                        critval = poisson.test(x = sum(baselinecounts), T = length(baselinecounts), conf.level = 0.999)$conf.int[1]   #computes the lower limit of the confindence interval
                        
                        diffs = diff(which(critwin < critval))       #computes the differences in indices for bins lower than the critical value (to check whether excited bins are consecutive or not)
                        cueinh = F
                        if(length(which(rle(diffs)$values == 1 & rle(diffs)$lengths >= threshold))>0) (cueinh = T)   #looks for consecutive bins (diff equal to 1) that are at least 2 bins long (rle length of at least 2)
                        
                        return(cueinh)  
                })
        })
        
        
        propcueinh <- lapply(seq(1, length(cueinhidx)), function(x){
                length(which(cueinhidx[[x]]==T))/length(cueinhidx[[x]])
        })
        
        
        
        
        #######################################################################
        #This function returns a list of lists:
        #Level 1 is the experiment
        #Level 2 contains the the frequency histogram (PSTH)
        #for each neuron for the event specified above.
        #
        #The object neurohist can then be used to plot histograms individually or
        #averaged together as in the next function.
        #########################################################################
        neurohist = lapply(seq(1, length(masterlist)), function(x) {
                lapply(seq(1, length(masterlist[[x]])), function(y) {
                        allvals = unlist(masterlist[[x]][[y]])
                        hcounts = hist(allvals, breaks = seq(-winmin, winmax, binsize/1000), plot = F)$counts 
                        freq = (hcounts/(binsize/1000))/length(masterlist[[x]][[y]])
                })
        })
        
        firingmat = do.call("cbind", lapply(seq(1,length(neurohist)), function(x) {
                do.call("cbind", neurohist[[x]]) }))
        
        #When there are no neurons in the selected side, firingmat comes out with columns with only zeroes
        #NonZeroColumns = colSums(firingmat)!=0
        #firingmat = firingmat[,NonZeroColumns]
        
        
        # Average firing of all the neurons in one session per bin and number of neurons per session
        meanFreqPerBin = lapply(seq(1, length(neurohist)), function(x){
                perBin <- as.data.frame(neurohist[[x]], col.names = c(1:length(neurohist[[x]])))
                mean <- rowMeans(perBin)
                return(mean)
        })
        
        SEM <- function(x){sd(x, na.rm=T)/sqrt(length(x))}
        semFreqPerBin = lapply(seq(1, length(neurohist)), function(x){
                perBin <- as.data.frame(neurohist[[x]], col.names = c(1:length(neurohist[[x]])))
                SdErrMean <- apply(perBin, 1, SEM)
                return(SdErrMean)
        })
        
        nNeurons = lapply(seq(1, length(neurohist)), function(x){
                n <- length(neurohist[[x]])
                return(n)
        })
        
        #### What data to use/show
        
        # Cue excited neurons only?
        if(cueexonly == T) ({
                cueexneurons = which(unlist(cueexidx) == T)
                firingmat = firingmat[,cueexneurons]
        })
        
        ### Parameters used 
        parameters = data.frame(path=path, startt=startt, endt=endt, binw=binw, psthmin=psthmin, psthmax=psthmax, event=event, cueexonly=cueexonly, side=side)
        
        
        if(allResults==T){
                return(list(nexdata=nexdata, neurons=neurons, events=events, 
                            masterlist=masterlist, cueexidx=cueexidx, propcueex=propcueex, 
                            cueinhidx=cueinhidx, propcueinh=propcueinh,
                            neurohist=neurohist, firingmat=firingmat, meanFreqPerBin=meanFreqPerBin, 
                            nNeurons=nNeurons, semFreqPerBin=semFreqPerBin, 
                            parameters=parameters))
                
        } else {return(firingmat)}
        
}

funcdirect <- paste(getwd(), "/R functions/", sep="")
save(neuralhist, file=paste(funcdirect, "neuralhist.R", sep=""))
