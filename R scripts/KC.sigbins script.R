KC.sigbins = function(path, startt, endt, event, BLwdw=10, PostEvent_wdw=3, pbin=0.05, funcdirect=funcdirect) {
        
        # Load necessary functions
        functionFolder = funcdirect
        load(file = paste(functionFolder, "eventfinder.r", sep=""))
        load(file = paste(functionFolder, "neuraldataextractMV.r", sep=""))
        
        #specify path of nex files
        neuralpath = path
        nexfiles = list.files(neuralpath)
        
        #extracts all NeuroExplorer data
        nexdata = lapply(nexfiles, function(x) {
                neuraldataextractMV(inpath=paste(path, x, sep = ""), funcdirect=funcdirect)})
        
        
        #these two functions split up the neurons and the behavioral events for each experiment
        neurons = lapply(seq(1, length(nexdata)), function(x) {
                neuronidx = grep("sig", names(nexdata[[x]]))
                nexdata[[x]][neuronidx] })
        
        events = eventfinder(nexdata, event)
        
        #if some nex files have no neurons, remove them here
        torem = which(unlist(lapply(neurons, function(x) length(x))) == 0)
        if(length(torem) != 0) ({
                neurons = neurons[-torem]
                events =  events[-torem]
                nexdata = nexdata[-torem]
        }) 
        
        
        #specify the interval of timestamps you want to examine  
        sstart = startt
        send = endt
        
        
        
        #########################################################
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
        
        
        spikelist = list()  #List in which each element is the vector of timestamps of each neuron
        for(i in 1:(length(neurons))) ({
                spikelist = c(spikelist, neurons[[i]])
        })  
        
        
        
        allevents = lapply(seq(1, length(events)), function(x) {
                rep(events[x], length(neurons[[x]]))
        })
        
        eventslist = list()  
        for(i in 1:(length(allevents))) ({
                eventslist = c(eventslist, allevents[[i]])
        })     
        
        
        
        #binsize = binw #in ms
        winmin = 10 #in s
        winmax = 10 #in s
        
        
        
        #####################################################################
        #This function returns a list of lists of lists:
        #Level 1 is the experiment; Level 2 is each neuron in the experiment;
        #Level 3 is the time relative to the specified event for each trial
        #The object masterlist can then be used to construct PSTHs or rasters
        ########################################################################
        masterlist = lapply(seq(1, length(spikelist)), function(x) {
                lapply(seq(1, length(eventslist[[x]])), function(y) {
                        stampsidx = which(spikelist[[x]] >= eventslist[[x]][y] - winmin & spikelist[[x]] <= eventslist[[x]][y] + winmax)
                        relstamps = spikelist[[x]][stampsidx] - eventslist[[x]][y]
                })
        }) 
        
        
        
        sigbins = lapply(seq(1, length(masterlist)), function(x) {
                
                
                intmin = -BLwdw #in s. Will be used as BL
                intmax = PostEvent_wdw #in s. Post event period of interest
                pbin = pbin #in s (50ms bin)
                
                allvals = unlist(masterlist[[x]])
                hcounts = hist(allvals[which(allvals >= intmin & allvals <= intmax)], breaks = seq(intmin, intmax, pbin), plot = F)$counts
                
                baseline = hcounts[1:abs(intmin/pbin)] #Counts per bin on the baseline period only (defined by intmin in bins of width defined by pbin)
                
                
                #Applying Vince's code to my data
                premean = mean( (baseline / length(masterlist[[x]])) )   # find the mean fr per bin within the pre period    
                probvect = ppois(hcounts, premean*length(masterlist[[x]])) # the lambda argument is the mean fr per bin times the number of trials
                
                critwin = probvect[(abs(intmin/pbin)+1):(abs(intmin/pbin)+intmax/pbin)]
                
                
                sigbins = rep(F, length(critwin))
                sigbins[which(critwin >= .999)] = T #If the probability of observing the count value exceeds a 99.9% confidence interval of our Poisson distribution (defined in ppois using the 10s pre-cue window), then flag that neuron as excited on that bin
                
                
                
                return(sigbins)
                
                
        })
        
        sigbindf = t(do.call("rbind", sigbins)) 
        
        return(sigbindf)
        
}


save(KC.sigbins, file="E:/Dropbox/DISSERTATION/R functions/KC.sigbins.R")
save(KC.sigbins, file="E:/Dropbox/NMDA/EXP1_Performance/R Functions/KC.sigbins.R")
save(KC.sigbins, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/KC.sigbins.R")
save(KC.sigbins, file="E:/Dropbox/NMDA/R Functions/KC.sigbins.R")


