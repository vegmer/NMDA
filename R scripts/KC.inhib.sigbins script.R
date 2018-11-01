KC.inhib.sigbins = function(path, funcdirect, startt, endt, event, BLwdw=10, PostEvent_wdw=3, pbin=0.05) {
        
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
        
        
        #specify the interval of timestamps you want to examine  
        sstart = startt
        send = endt
        
        #specify whether to only return cue excited neurons
        
        
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
        
        
        
        #########################################################
        #Now that the relevant data is extracted, we can begin to 
        #analyze it by constructing rasters, histograms, etc.
        #########################################################
        #binsize = binw #in ms
        winmin = 10 #in s
        winmax = 10 #in s
        
        
        
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
        
        
        
        sigbins = lapply(seq(1, length(masterlist)), function(x) {
                lapply(seq(1, length(masterlist[[x]])), function(y) {    
                        
                        intmin = -BLwdw #in s. Will be used as BL
                        intmax = PostEvent_wdw #in s. Post event period of interest
                        pbin = pbin #in s (50ms bin)
                        
                        allvals = unlist(masterlist[[x]][[y]])
                        hcounts = hist(allvals[which(allvals >= intmin & allvals <= intmax)], breaks = seq(intmin, intmax, pbin), plot = F)$counts
                        
                        baseline = hcounts[1:abs(intmin/pbin)] #Counts on the baseline period only
                        
                        
                        #Applying Vince's code to my data
                        premean = mean( (baseline / length(masterlist[[x]][[y]])) )   # find the mean fr per bin within the pre period  (bc "baseline" adds the counts on all of the trials throughout the session)  
                        probvect = ppois(hcounts, premean*length(masterlist[[x]][[y]])) # the lambda argument is the mean fr per bin times the number of trials. That gives us the count of the number of spikes in the typical baseline bin. The ppois function gives us the probability of the observed counts given a Poisson distribution with the lambda we entered.
                        
                        critwin = probvect[(abs(intmin/pbin)+1):(abs(intmin/pbin)+intmax/pbin)]
                        
                        
                        sigbins = rep(F, length(critwin))
                        sigbins[which(critwin <= .01)] = T #If the probability of observing the count value is lower than the 0.01% confidence interval of our Poisson distribution (defined in ppois using the 10s pre-cue window), then flag that neuron as inhibited on that bin
                        
                        
                        
                        return(sigbins)
                        
                })
        })
        
        sigbindf = t(do.call("rbind", lapply(seq(1, length(sigbins)), function(x) do.call("rbind", sigbins[[x]])))) 
        
        return(sigbindf)
        
}


save(KC.inhib.sigbins, file="E:/Dropbox/DISSERTATION/R functions/KC.inhib.sigbins.R")
save(KC.inhib.sigbins, file="E:/Dropbox/NMDA/EXP1_Performance/R Functions/KC.inhib.sigbins.R")
save(KC.inhib.sigbins, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/KC.inhib.sigbins.R")
save(KC.inhib.sigbins, file="E:/Dropbox/NMDA/R Functions/KC.inhib.sigbins.R")

