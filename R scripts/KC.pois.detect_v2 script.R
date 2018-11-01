KC.pois.detect_v2 = function(spikes, events) {
        
        #load("C:/Users/Kevin Caref/Google Drive/RScripts/Functions/eventfinder.rFunc")
        #
        ##these two functions split up the neurons and the behavioral events for each experiment
        #neurons = lapply(seq(1, length(nexdata)), function(x) {
        #  neuronidx = grep("sig", names(nexdata[[x]]))
        #  nexdata[[x]][neuronidx] })
        #  
        #events = eventfinder.rFunc(nexdata, 1) #for this function, we will find cue-excited neurons using behavioral event 1, all CS+ cues
        #
        # 
        #
        #
        ##specify the interval of timestamps you want to examine  
        #sstart = startt
        #send = endt
        #
        ##specify whether to only return cue excited neurons
        #    
        #
        ##########################################################
        ##these two functions simply restrict the timestamp values
        ##to the window specified above
        ########################################################
        #neurons = lapply(seq(1, length(neurons)), function(x) {
        #  lapply(seq(1, length(neurons[[x]])), function(y) {
        #    neurons[[x]][[y]][which(neurons[[x]][[y]] >= sstart & neurons[[x]][[y]] <= send)]
        #   })
        #  })
        #
        #events = lapply(seq(1, length(events)), function(x) {
        #   events[[x]][which(events[[x]] >= sstart & events[[x]] <= send)] })
        #  
        #spikelist = list()  
        #for(i in 1:(length(neurons))) ({
        #  spikelist = c(spikelist, neurons[[i]])
        #  })  
        #
        #
        #
        #allevents = lapply(seq(1, length(events)), function(x) {
        #  rep(events[x], length(neurons[[x]]))
        #  })
        #  
        #eventslist = list()  
        #for(i in 1:(length(allevents))) ({
        #  eventslist = c(eventslist, allevents[[i]])
        #  })     
        #  
        
        #########################################################
        #Now that the relevant data is extracted, we can begin to 
        #analyze it by constructing rasters, histograms, etc.
        #########################################################
        #binsize = binw #in ms
        winmin = 10 #in s
        winmax = .5 #in s
        
        
        
        #####################################################################
        #This function returns a list of lists of lists:
        #Level 1 is the experiment; Level 2 is each neuron in the experiment;
        #Level 3 is the time relative to the specified event for each trial
        #The object masterlist can then be used to construct PSTHs or rasters
        ########################################################################
        masterlist = lapply(seq(1, length(spikes)), function(x) {
                lapply(seq(1, length(events[[x]])), function(y) {
                        stampsidx = which(spikes[[x]] >= events[[x]][y] - winmin & spikes[[x]] <= events[[x]][y] + winmax)
                        relstamps = spikes[[x]][stampsidx] - events[[x]][y]
                })
        }) 
        
        
        
        
        cueex = lapply(seq(1, length(masterlist)), function(x) {
                
                
                intmin = -10
                intmax = .5
                pbin = .05
                threshold = 1
                
                allvals = unlist(masterlist[[x]])
                hcounts = hist(allvals[which(allvals >= intmin & allvals <= intmax)], breaks = seq(intmin, intmax, pbin), plot = F)$counts
                
                baseline = hcounts[1:abs(intmin/pbin)]
                
                
                #Applying Vince's code to my data
                premean = mean( (baseline / length(masterlist[[x]])) )   # find the mean fr per bin within the pre period    
                probvect = ppois(hcounts, premean*length(masterlist[[x]])) # the lambda argument is the mean fr per bin times the number of trials
                
                critwin = probvect[(abs(intmin/pbin)+1):(abs(intmin/pbin)+intmax/pbin)]
                
                
                
                #diffs = diff(which(critwin >= .999))       #computes the differences in indices for bins exceeding the critical value
                cueex = F
                #if(length(which(rle(diffs)$values == 1 & rle(diffs)$lengths >= threshold))>0) (cueex = T)   #looks for consecutive bins (diff equal to 1) that are at least 3 bins long (rle length of at least 3)
                if(length(which(critwin >= .999)) > 0) (cueex = T)
                
                return(cueex)
                
                
        })
        
        return(cueex)
        
}


save(KC.pois.detect_v2, file="E:/Dropbox/DISSERTATION/R functions/KC.pois.detect_v2.R")
save(KC.pois.detect_v2, file="E:/Dropbox/NMDA/EXP1_Performance/R Functions/KC.pois.detect_v2.R")
save(KC.pois.detect_v2, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R Functions/KC.pois.detect_v2.R")
save(KC.pois.detect_v2, file="E:/Dropbox/NMDA/R Functions/KC.pois.detect_v2.R")
