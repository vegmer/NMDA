neuraldataextractMV <- function (inpath, funcdirect) {
        
        load(paste(funcdirect, "/nexreaderwaves.f", sep=""))
        
        filename = inpath
        
        starttimeVec <- list()
        eventstampslength <- list()
        neuralstampslength <- list()
        
        for(i in 1:length(filename)){
                
                #filepath = paste(NEXfilesSource, files[i], sep="")
                fullname =  filename[i]
                file = basename(fullname)
                underscPos = gregexpr("_", file)[[1]] #Position of underscores in the filename
                dotPos = gregexpr(".nex", file)[[1]]
                RatPos <- underscPos[1]+1 #Position of rat's number
                SessPos <- underscPos[2]+8 #Position of session number
                KPulsePos <- underscPos[3]+1 #Position of info about task version (98,99)
                InfusionPos <- underscPos[4]+1 #Position of info about drug infusion: (L) left side; (R) right; (N) none; (B) both
                Drovepos <- underscPos[5]+1 #Position of info about driving down electrodes or not
                
                
                exptdate <- paste(substr(file, 1, 2), "/", substr(file, 3, 4), "/", substr(file, 5, 6), sep="")
                ratname <- substr(file, RatPos, underscPos[2]-1)
                expt <-  substr(file, SessPos, underscPos[3]-1)
                taskVersion <- substr(file, KPulsePos, underscPos[4]-1)
                drugSide <- substr(file, InfusionPos, underscPos[5]-1)
                droveDown <- substr(file, Drovepos,dotPos[1]-1)
                
                nexfile <- nexreaderwaves.f(paste(fullname, sep = ""))
                neuralindx = grep("sig", names(nexfile))
                eventsindx = names(nexfile)[-c(neuralindx, grep("AllFile", names(nexfile)), grep("frequency", names(nexfile)))] #I removed from the index all the neurons and the variables that were not recorded with timestamps
                eventdata = nexfile[eventsindx]
                neuraldata = nexfile[neuralindx]
                eventstamps = lapply(eventdata, function(x) x$tstamp)
                neuralstamps = lapply(neuraldata, function(x) x$tstamp)
                neuronnames = names(neuralstamps)
                eventnames = names(eventstamps)
                
                if(length(eventstamps$SessionStart) < 1) (eventstamps$SessionStart = 0)
                
                #Substract SESSION START timestamp from all the other timestamps (so that data is "aligned" to session start)
                #There's one file (MV184, expt 3) in which, for some reason, I have more than one value for starttime. I have 28 values and the first 27 are BS (they are regular pulses at 3.77s. I checked and the last one is the right one, so choose that one)
                starttime = eventstamps$SessionStart
                if(length(starttime)>1){starttime <- last(starttime)}
                eventstamps = lapply(seq(1, length(eventstamps)), function(x) eventstamps[[x]] - starttime)
                names(eventstamps) <- eventnames
                correctedstamps = lapply(seq(1, length(neuralstamps)), function(x) neuralstamps[[x]] - starttime)
                neuralstamps = lapply(seq(1, length(neuralstamps)), function(x) correctedstamps[[x]][which(correctedstamps[[x]] > 0)])
                names(neuralstamps) <- neuronnames
                
                starttimeVec[[i]] <- list(starttime)
                eventstampslength[[i]] <- list(eventstamps)
                neuralstampslength[[i]] <- list(neuralstamps)
                
                
                #I just want to include the previous events in my "eventstamps" object (eliminate irrelevant events)
                #eventstamps = eventstamps[-which(is.na(names(eventstamps)) == T)]
                
                eventstamps$ratname = ratname
                eventstamps$date = exptdate
                eventstamps$filename = filename
                eventstamps$expt = expt
                eventstamps$drug = drugSide
                eventstamps$down = droveDown
                eventstamps$version= taskVersion
                
                return(c(eventstamps, neuralstamps))
        }
}



save(neuraldataextractMV, file="E:/Dropbox/NMDA/EXP1_Performance/R Functions/neuraldataextractMV.R")
save(neuraldataextractMV, file="E:/Dropbox/NMDA/R Functions/neuraldataextractMV.R")
