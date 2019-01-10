### FUNCTION THAT EXTRACTS AND SAVES THE FOLLOWING BEHAVIOR-RELATED DATA OBJECTS:

# a) In the "datafolder" folder:
#       csacqidx -> index of all sessions with summary indices
#       alldata -> list of objects with detailed about relevant session events (each object corresponds to each session, the number of each object corresponds to the row of the same number in 'csacqidx')
#       rats -> unique names of all rats
#       idx -> list where each object corresponds to one rat (the object rats will tell you which rat). Inside each object, there is a vector with the indices of each one of that animal's sessions in csacqidx and alldata

# b) In the "dataForRCumulative" folder:
# DSlatency -> on every S+ trial: latency to enter the compartment after S+ onset. If the animal didn't enter, a 10 (max. length of cue) was assigned. This is because the
# task accuracy index (or "performance index") is based on ITIlatency-DSlatency and, to be able to do that, 10 is assigned when the animal doesn't respond. But if you're 
# just interested in latency, then swith the values that equal 10 for NA.
# DSrespAll -> every DS trial of every session: 1=responded to the cue; 0=ignored the cue
# DStaskAcc -> How much faster was the animal to enter the rec. after the cue vs. during the last 10s of the ITI on that trial (ITIlatency or ITIpseudolatency)
# DStimeToSpare -> opposite of DSlatency (10-latency to enter after DS onset)
# ITIlatency -> also known as "ITI pseudolatency". Latency to enter into the compartment in the 10 s ITI window that precedes the cue.Note that, in this index, all trials are included (ITI windows that precede DSs and NSs. Further analyses break this down by cue)
# NSlatency -> Same as DS latency but for NS (or S-) trials
# NSrespAll -> Same thing as DS response ratio but for NS
# NStaskAcc -> same as DS task accuracy but for NS (this is also known as "performance index")
# NStimeToSpare -> same as DStimeToSpare but for NS trials
# ITIrespRatio -> Whether the animal made an entry during the 10s ITI window that preceded the cue or not

# The parameters of the function are:
# "MovAvg" --> Select how to calculate the ITI pseudolatency. Options: "All", "No" or "Impinged only"
#    a) "All":To each trial, assign the average ITI pseudolatency in the 5 trial window around it (5 trials before, 5 trials after and the actual trial); 
#    b) "No": to each trial assign its own ITI pseudolatency, whatever it is. If the animal was already inside the receptacle when the ITI window started, it won't count that. It will only count the first "active" entry after the beginning of the ITI period.
#    c) "Impinged only": In trials in which the animal's head was already inside at the time the ITI pseudolatency period started, assign the moving average of the 5 trial window around that trial. If the head of the animal was not inside at that time, just assign the ITI pseudolatency value of that trial.
# "funcdirect" -->  folder where the functions that are called by the function are
# "datafolder" -->  folder with the MedPC files
# "dataForRdir" -->  folder where the function is going to store the objects it creates
# "dataForRCumulative" -->  folder where the function is going to store objects that will be used later for cumulative performance graphs (trial by trial performance by animal in different indices)
# "cuelength" --> length of the cue in case the animal doesn't make an entry (in seconds)

MedPCextract <- function(MovAvg="Impinged only", funcdirect=funcdirect, datafolder=datafolder, 
                         dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative, cuelength){
        
        load (paste(funcdirect,"CStask.Rfunc",sep=""))
        load (paste(funcdirect,"mpcextract_blockSingle.Rfunc",sep=""))
        load (paste(funcdirect,"mpcextractSingle.Rfunc",sep=""))
        
        mpcdatadir = datafolder
        allfiles = list.files(mpcdatadir)
        
        #Create list alldata in which each element is the data of each animal during a session
        alldata = list()
        for (i in(1:length(allfiles))) { 
                
                rawmpc = mpcextractSingle.Rfunc(paste(mpcdatadir, allfiles[i], sep = ""))
                toanalyze = CStask.Rfunc(rawmpc)
                
                EXPTindx=grep("Experiment",rawmpc[[1]])
                EXPT=substr (rawmpc[[1]][EXPTindx],13,nchar(rawmpc[[1]][EXPTindx]))
                testname=unlist (strsplit(rawmpc$info[1],"Subject"))
                ratname=substr (testname[2],2,nchar(testname[2]))
                idate=unlist (strsplit(rawmpc$info[1],"!"))
                testdate=substr(idate[2],1,10)
                
                group= substr(toanalyze$info[6], 8, 17)
                CSpluscue = toanalyze$CSpluscue
                CSminuscue = toanalyze$CSminuscue
                receptacleentries = toanalyze$receptacleentries
                receptacleexits = toanalyze$receptacleexits
                rewarddelivery = toanalyze$rewarddelivery 
                
                CSplusresponse = vector()
                for(j in 1:length(CSpluscue)) {
                        CSplusentries = receptacleentries[receptacleentries > CSpluscue[j] & receptacleentries < CSpluscue[j]+cuelength]
                        if (length(CSplusentries) > 0)(CSplusresponse[j] = min(CSplusentries))
                        else(CSplusresponse[j] = NA)
                }
                
                #CSplusresponse = CSplusresponse[-which(is.na(CSplusresponse))]
                
                correctTrials <- findInterval(CSplusresponse[-which(is.na(CSplusresponse))], toanalyze$CSpluscue)
                
                #In July 2017 I changed the task so that the cue extends 2s after correct entry, so calculating the CSplusEnd will be different after 07/01/2017
                CSplusEnd = vector()
                
                if(as.Date(substr(toanalyze$info[2], 13, 20), format='%m/%d/%y')<as.Date("07/01/17", format="%m/%d/%y")) #Before change
                {
                        for(j in 1:length(CSpluscue)){
                                if(sum(j==correctTrials)>0){
                                        selTrial <- as.numeric(correctTrials[j==correctTrials])
                                        CSplusEnd[j] = toanalyze$rewarddelivery[selTrial]
                                } else {
                                        CSplusEnd[j] = toanalyze$CSpluscue[j] + cuelength
                                }
                        }
                } else { #after change
                        for(j in 1:length(CSpluscue)){
                                if(sum(j==correctTrials)>0){
                                        selTrial <- as.numeric(correctTrials[j==correctTrials])
                                        CSplusEnd[j] = (CSplusresponse[selTrial])+2
                                } else {
                                        
                                        CSplusEnd[j] = toanalyze$CSpluscue[j] + cuelength
                                }
                        }
                }
                
                CSminusEnd = vector()
                for(j in 1:length(CSminuscue)){
                        CSminusEnd[j] = toanalyze$CSminuscue[j] + cuelength
                }
                
                
                CSminusresponse = vector()
                for(k in 1:length(CSminuscue)) {
                        CSminusentries = receptacleentries[receptacleentries > CSminuscue[k] & receptacleentries < CSminuscue[k]+cuelength]
                        if (length(CSminusentries) > 0)(CSminusresponse[k] = min(CSminusentries))
                        else(CSminusresponse[j] = NA)
                }
                
                CSminusresponseONLY = CSminusresponse[-which(is.na(CSminusresponse))]
                
                # Duration of entries
                if(length(receptacleentries) > length(receptacleexits)){receptacleentries <- receptacleentries[-(length(receptacleentries))]} #Drop last entry in case session ended mid entry
                entryDur <- receptacleexits - receptacleentries #Duration of all entries
                RewEnt <- receptacleentries %in% CSplusresponse #Index of CS Plus entries (rewarded)
                CSMinusEnt <- receptacleentries %in% CSminusresponse #Index of CS Minus entries
                ITIEnt <- RewEnt==F & CSMinusEnt==F #Index of ITI entries
                
                CSPlusEntryDur <- entryDur[RewEnt]
                CSMinusEntryDur <- entryDur[CSMinusEnt]
                ITIentryDur <- entryDur[ITIEnt]
                
                
                allCues = sort(c(toanalyze$CSpluscue, toanalyze$CSminuscue))
                allCueEnds = sort(c(CSplusEnd, CSminusEnd))
                
                ITIlength <- vector()
                for(l in 1:length(allCues)){
                        if(l==1){
                                ITIlength[l]<-allCues[l]
                        }else{
                                ITIlength[l]<-allCues[l]-allCueEnds[l-1]
                        }
                }
                
                lengthCue <- vector()
                for(m in 1:length(allCues)){
                        lengthCue[m]<- allCueEnds[m]-allCues[m]
                }
                
                CSplusLat <- CSplusresponse - CSpluscue
                CSminusLat <- CSminusresponse - CSminuscue
                
                trialStart <- c(0, allCueEnds[1:length(allCueEnds)-1])
                
                allCuesOrder <- allCues
                
                allCuesOrder[is.element(allCues, toanalyze$CSpluscue)] <- 1
                allCuesOrder[is.element(allCues, toanalyze$CSminuscue)] <- 2
                
                firstTr <- 1
                
                ITIlatency <- sapply(seq(firstTr, length(ITIlength)), function(x){
                        ITIwdwStarts <-allCues[x]-cuelength
                        ITIwdwEnds <- allCues[x]
                        
                        ITIlat <- sapply(seq(1, length(ITIwdwStarts)), function(y){
                                ITIwdwEntries <- receptacleentries[receptacleentries>=ITIwdwStarts[y] & receptacleentries<=ITIwdwEnds[y]]
                                if(length(ITIwdwEntries)==0){ITIlat <- cuelength}
                                else {
                                        firstEnt <- min(ITIwdwEntries)
                                        ITIlat <- firstEnt-ITIwdwStarts[y]}
                        })
                        
                        return(ITIlat)
                        
                })
                
                
                #Go to cued entry latency values and replace NA values with the length of the cue
                CSplusLat[is.na(CSplusLat)] <- cuelength
                CSminusLat[is.na(CSminusLat)] <- cuelength
                
                #ITI latency - cued latency. If 1st cue is long CS+ after 5s, get rid of that trial for this calculation (bc the ITI is too short, 5s)
                #if(ITIlength[1]<=6){allCuesOrder2 <- allCuesOrder[-1]; CSplusLat2 <- CSplusLat[-1]} else {
                #        allCuesOrder2 <- allCuesOrder; CSplusLat2 <- CSplusLat}
                
                allCuesOrder2 <- allCuesOrder; CSplusLat2 <- CSplusLat
                CSplusTA <- ITIlatency[which(allCuesOrder2==1)] - CSplusLat #TA stands for "Task Accuracy", which I also later call performance index
                CSminusTA <- ITIlatency[which(allCuesOrder2==2)] - CSminusLat
                
                alldata[[i]] = list(ratname = ratname, testdate = testdate, expt = EXPT, group=group, receptacleentries = receptacleentries, receptacleexits = toanalyze$receptacleexits,
                                    CSpluscue = CSpluscue, CSplusEnd=CSplusEnd, CSplusresponse = CSplusresponse,rewarddelivery = toanalyze$rewarddelivery, 
                                    CSplusLat=CSplusLat, CSminuscue = CSminuscue, CSminusresponse = CSminusresponse, CSminusEnd = CSminusEnd, CSminusLat=CSminusLat, 
                                    allCues=allCues, allCueEnds=allCueEnds, CueLength=lengthCue, ITIlength=ITIlength, trialStart=trialStart, orderCues = allCuesOrder, 
                                    CSPlusEntryDur=CSPlusEntryDur, CSMinusEntryDur=CSMinusEntryDur, ITIentryDur=ITIentryDur, ITIlatency=ITIlatency, CSplusTA=CSplusTA, 
                                    CSminusTA=CSminusTA)
                
                
        }
        
        allrats=unlist (lapply(alldata,FUN=function(x)({
                toget=x$ratname
                return(toget)
        })))
        
        allexpts=unlist (lapply(alldata,FUN=function(x)({
                toget2=x$expt
                return(toget2)
        })))
        
        alldates=unlist (lapply(alldata,FUN=function(x)({
                toget3=x$testdate
                return(toget3)
        })))     
        
        allgroups=unlist (lapply(alldata,FUN=function(x)({
                toget3=x$group
                return(toget3)
        })))     
        
        #Create empty objects for the idx and the data
        csacqidx = data.frame(subject = allrats, session = allexpts, group = allgroups, filename = allfiles, date = alldates)
        
        
        ##############################################################
        #ANALYSIS 1: BASIC PARAMETERS PER SESSION PER ANIMAL       ###
        ##############################################################
        
        DSperc <- c()
        NSperc <- c()
        ITItotal <- c()
        ITIperSec <- c()
        ITIlatency <- c()
        latencyCSplus <- c()
        latencyCSminus <- c()
        latencyCSplusSEM <- c()
        latencyCSminusSEM <- c()
        DSentryDur <- c()
        NSentryDur <- c()
        CSplusTaskAcc <- c()
        CSminusTaskAcc <- c()
        ITIPeriodentryDur <- c()
        PercShortDSent <- c()
        
        for (i in 1:length(alldata)){
                DSperc [[i]]<- length(alldata[[i]]$rewarddelivery)/length(alldata[[i]]$CSpluscue)
                NSperc [[i]]<- length(alldata[[i]]$CSminusresponse[which(alldata[[i]]$CSminusresponse>0)])/length(alldata[[i]]$CSminuscue)
                ITItotal[[i]] <- sum(alldata[[i]]$ITIlength)
                ITIperSec[[i]] <- round((length(alldata[[i]]$receptacleentries) - length(c(alldata[[i]]$CSplusresponse, alldata[[i]]$CSminusresponse)))/ITItotal[[i]], 2)
                #if(ITIperSec[[i]]<0){ITIperSec[[i]]=0} else {ITIperSec[[i]]=ITIperSec[[i]]}
                ITIlatency[[i]] <- mean(alldata[[i]]$ITIlatency, na.rm=T)
                latencyCSplus[[i]] <- mean(alldata[[i]]$CSplusLat, na.rm=T)
                latencyCSminus[[i]] <- mean(alldata[[i]]$CSminusLat, na.rm=T)
                latencyCSplusSEM[[i]] <-sd(alldata[[i]]$CSplusLat, na.rm=T)/sqrt(length(alldata[[i]]$CSplusLat[which(alldata[[i]]$CSplusLat >0)]))
                latencyCSminusSEM[[i]] <-sd(alldata[[i]]$CSminusLat, na.rm=T)/sqrt(length(alldata[[i]]$CSminusLat[which(alldata[[i]]$CSminusLat >0)]))
                
                DSentryDur[[i]] <- mean(alldata[[i]]$CSPlusEntryDur, na.rm=T)
                NSentryDur[[i]] <- mean(alldata[[i]]$CSMinusEntryDur, na.rm=T)
                
                CSplusTaskAcc[[i]] <- mean(alldata[[i]]$CSplusTA, na.rm=T)
                CSminusTaskAcc[[i]] <- mean(alldata[[i]]$CSminusTA, na.rm=T)
                
                ITIPeriodentryDur[[i]] <- mean(alldata[[i]]$ITIentryDur, na.rm=T)
                PercShortDSent[[i]] <- sum(alldata[[i]]$CSPlusEntryDur<=2)/length(alldata[[i]]$CSPlusEntryDur)
        }
        
        #Now I'm adding a column to csacqidx with the rate of HE of each animal each day
        csacqidx["DSperc"]<- DSperc
        csacqidx["NSperc"]<- NSperc
        csacqidx["ITItotal"]<- ITItotal
        csacqidx["ITIpersec"]<- ITIperSec
        csacqidx["ITIlatency"] <- ITIlatency
        csacqidx["LatencyCSplus"] <- latencyCSplus
        csacqidx["LatencyCSminus"] <- latencyCSminus
        csacqidx["LatencyCSplusSEM"] <- latencyCSplusSEM
        csacqidx["LatencyCSminusSEM"] <- latencyCSminusSEM
        csacqidx["CSplusTaskAcc"] <- CSplusTaskAcc
        csacqidx["CSminusTaskAcc"] <- CSminusTaskAcc
        csacqidx["DSentryDur"] <- DSentryDur
        csacqidx["NSentryDur"] <- NSentryDur
        csacqidx["ITIentryDur"] <- ITIPeriodentryDur
        csacqidx["PercShortDSent"]<- PercShortDSent
        
        
        #index to select sessions of the same animals
        rats <- unique(csacqidx$subject)
        v <- 1:length(csacqidx$subject)
        idx<-c()
        for(i in 1:length(rats)){
                idx[[i]] <- v[csacqidx$subject==rats[i]]
        }
        
        ### OBJECTS TO USE FOR CUMULATIVE DATA
        #Make a list with one object per rat. Each object is a vector with probability of response (1 or 0) to every single trial of every session of that rat (in order)
        DSrespAll <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        DSresp <- (!is.na(alldata[[sessIdx]]$CSplusresponse))*1
                }))
        })
        
        NSrespAll <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        NSresp <- (!is.na(alldata[[sessIdx]]$CSminusresponse))*1
                }))
        })
        
        #BC of an issue with MedPC, sometimes the last trial was not recorded or MedPC allowed one more trial to be presented. Different lengths of DS and NS give me errors, so make them the same length (length of the shortest one)
        DSlength <- sapply(DSrespAll, length)
        NSlength <- sapply(NSrespAll, length); lengthIdx <- NSlength #Bc NSlength is always shorter
        
        DSrespAll <- lapply(seq(1, length(DSrespAll)), function(k){DSrespAll[[k]][1:lengthIdx[k]]})
        
        #DSlatency
        DSlatency <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        DSresp <- (!is.na(alldata[[sessIdx]]$CSplusresponse))*1
                        DSlat <- alldata[[sessIdx]]$CSplusLat
                        DSlat[which(is.na(DSlat))] <- cuelength
                        DSlat
                }))
        })
        
        DSlatency <- lapply(seq(1, length(DSlatency)), function(k){DSlatency[[k]][1:lengthIdx[k]]})
        
        #DStimetoSpare (opposite of latency)
        DStimeToSpare <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        DSresp <- (!is.na(alldata[[sessIdx]]$CSplusresponse))*1
                        DSlat <- alldata[[sessIdx]]$CSplusLat
                        DSlat[which(is.na(DSlat))] <- cuelength
                        timeToSpare <- cuelength-DSlat
                }))
        })
        
        DStimeToSpare <- lapply(seq(1, length(DStimeToSpare)), function(k){DStimeToSpare[[k]][1:lengthIdx[k]]})
        
        
        #NStimeToSpare
        NStimeToSpare <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        NSresp <- (!is.na(alldata[[sessIdx]]$CSminusresponse))*1
                        NSlat <- alldata[[sessIdx]]$CSminusLat
                        NSlat[which(is.na(NSlat))] <- cuelength
                        timeToSpare <- cuelength-NSlat
                }))
        })
        
        #NSlatency
        NSlatency <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        NSresp <- (!is.na(alldata[[sessIdx]]$CSminusresponse))*1
                        NSlat <- alldata[[sessIdx]]$CSminusLat
                        NSlat[which(is.na(NSlat))] <- cuelength
                        NSlat
                }))
        })
        
        #DS task accuracy (latency before the cue vs latency after the cue). How much faster to enter the receptacle after the cue vs before the cue
        DStaskAcc <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        DStaskAcc <- alldata[[sessIdx]]$CSplusTA
                        DStaskAcc <- DStaskAcc[!is.na(DStaskAcc)]#Remove NAs bc it won't work otherwise (there shouldn't be any)
                }))
        })
        
        DStaskAcc <- lapply(seq(1, length(DStaskAcc)), function(k){DStaskAcc[[k]][1:lengthIdx[k]]})
        
        
        #NS task accuracy
        NStaskAcc <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        NStaskAcc <- alldata[[sessIdx]]$CSminusTA
                        NStaskAcc <- NStaskAcc[!is.na(NStaskAcc)] #Remove NAs bc it won't work otherwise (there shouldn't be any)
                }))
        })
        
        
        ITIlatency <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        ITIlatency <- alldata[[sessIdx]]$ITIlatency
                        ITIlatency
                }))
        })
        
        #### Now I'm going to correct the ITI latency of problematic trials as specified by the parameter "MovAvg"
        
        ### Calculate ITI latency using a moving average of +/- 5 trials
        cuesOrderAll <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        cuesOrder <- alldata[[sessIdx]]$orderCues
                }))
        })
        
        ITIlatMovingAvg <- lapply(seq(1, length(rats)), function(k){
                
                #Number of trials around the trial of interest that I want to use for the moving window. Change this parameter to change the number of trials in that window
                nTrialWdw <- 5 
                nTrials <- length(ITIlatency[[k]])
                
                #Let's say the moving average is the moving average of 5 trials. If the trial in question is number 3, let's say, 
                #then I can't choose the 5 trials before it. In that case, average the first 11 trials and assign that value to 
                #trial 3 (bc in trials in the middle of the session, I average 11 trials to get the value: the 5 before, 5 after 
                #and the trial in question). Same thing if the trial is one of the last ones.
                
                sapply(seq(1, nTrials), function(l){
                        if(l<(1+nTrialWdw)){
                                trlIdx <- 1:((nTrialWdw*2)+1)
                        } 
                        if(l>(nTrials-nTrialWdw)){
                                trlIdx <- (nTrials-(nTrialWdw*2)):nTrials
                        }
                        if(l>nTrialWdw & l<=(nTrials-nTrialWdw)){
                                trlIdx <- (l-nTrialWdw):(l+nTrialWdw)
                        }
                        
                        avg <- mean(ITIlatency[[k]][trlIdx], na.rm=T)
                })
        })
        
        #Task accuracy in S+ trials using the moving average
        DStaskAccMA <- lapply(seq(1, length(rats)), function(k){
                DS_ITIlatMA <- ITIlatMovingAvg[[k]][cuesOrderAll[[k]]==1]
                DS_ITIlatMA <- DS_ITIlatMA[1:lengthIdx[k]]
                
                DStaskAccMA <- DS_ITIlatMA-DSlatency[[k]]
        })
        
        #Task accuracy in S- trials using the moving average
        NStaskAccMA <- lapply(seq(1, length(rats)), function(k){
                NS_ITIlatMA <- ITIlatMovingAvg[[k]][cuesOrderAll[[k]]==2]
                NStaskAccMA <- NS_ITIlatMA-NSlatency[[k]]
        })
        
        #If no moving average is to be used to calculate ITI pseudolatency. This is not necessary but just to represent every possible value of the parameter "MovAvg".
        if(MovAvg=="No"){
                DStaskAcc <- DStaskAcc
                NStaskAcc <- NStaskAcc
                ITIlatency <- ITIlatency
        }
        
        #If moving average is to be used to calculate ITI pseudolatency in all trials.
        if(MovAvg=="All"){
                ITIlatency <- ITIlatMovingAvg
                DStaskAcc <- DStaskAccMA
                NStaskAcc <- NStaskAccMA
        }
        
        #If moving average is to be used to calculate ITI pseudolatency ONLY in trials in which the animal's head was already inside
        if(MovAvg=="Impinged only"){
                
                #Just for calculating how pervasive the problem of the rat being inside the compartment in the ITI pseudolat. window is. This gives me an index of the trials in which the animals were inside the compartment at the time the critical "ITI window" started
                Impinge <- lapply(seq(1, length(alldata)), function(x){
                        entries <- alldata[[x]]$receptacleentries
                        exits <- alldata[[x]]$receptacleexits
                        
                        #Pairs of entries-exits
                        pairs <- lapply(seq(1, length(entries)), function(y){
                                a <- c(entries[y], exits[y])
                        })
                        trialStarts <- alldata[[x]]$trialStart
                        ITIprecue <- alldata[[x]]$allCues-10
                        
                        #Index: ITI window that started in between what two entries. Substitute any 0 value (usually the first one) by 1, so that when we use this as an index we don't lose one trial
                        ITItoEntryIdx <- findInterval(ITIprecue, entries)
                        ITItoEntryIdx[ITItoEntryIdx==0]<- 1
                        
                        #Subset the pairs of entries-exits that might be problematic (an ITI window happened between two consecutive entries. Now we want to test if the rat exited before the start of the ITI window)
                        toTest <- pairs[ITItoEntryIdx]
                        
                        #Select trials based on what cue came on before (modify accordingly). 1=DS, 2=NS.
                        #preTrialCue <- c(0, alldata[[x]]$orderCues[-length(alldata[[x]]$orderCues)])
                        
                        ImpingeTrials <- sapply(seq(1, length(toTest)), function(z){
                                test <- ITIprecue[z]>toTest[[z]][1] & ITIprecue[z]<toTest[[z]][2] 
                                return(test)
                        })
                        return(ImpingeTrials)
                })
                
                #This will give me the index of trials in which the animals were inside the receptacle by rat
                ImpingeIdx <- lapply(seq(1, length(rats)), function(k){
                        selSess <- idx[[k]]
                        do.call("c", lapply(seq(1, length(selSess)), function(l){
                                Impinge[[selSess[l]]]
                        }))
                })
                
                #For each rat, use normal ITI latency (the one that is calculated with the precue window alone) and, in the trials in which the animal's head is still in at the beginning of the ITI window, use the average of the +/-5 trials around the trial
                ITIlatCorrected <- lapply(seq(1, length(rats)), function(k){
                        sapply(seq(1, length(ImpingeIdx[[k]])), function(l){
                                if(ImpingeIdx[[k]][l]==FALSE){ITIlat <- ITIlatency[[k]][l]} #If the head of the animal was NOT inside at the beginning of the window, just choose the normal ITI latency for that trial
                                if(ImpingeIdx[[k]][l]==TRUE){ITIlat <- ITIlatMovingAvg[[k]][l]} #If his head was inside at the time the ITI window started, then replace the value with the moving average for that trial (the +/- 5 trial average)
                                ITIlat
                        })
                })
                
                DStaskAccCorrected <- lapply(seq(1, length(rats)), function(k){
                        DS_ITIlatCorr <- ITIlatCorrected[[k]][cuesOrderAll[[k]]==1]
                        DS_ITIlatCorr <- DS_ITIlatCorr[1:lengthIdx[[k]]]
                        DStaskAccMA <- DS_ITIlatCorr-DSlatency[[k]]
                })
                
                NStaskAccCorrected <- lapply(seq(1, length(rats)), function(k){
                        NS_ITIlatCorr <- ITIlatCorrected[[k]][cuesOrderAll[[k]]==2]
                        NStaskAccMA <- NS_ITIlatCorr-NSlatency[[k]]
                })
                
                #Assign the new corrected values to ITI latency, DS task accuracy and NS task accuracy
                ITIlatency <- ITIlatCorrected
                DStaskAcc <- DStaskAccCorrected
                NStaskAcc <- NStaskAccCorrected
        }
        
        
        
        ###### SAVING OBJECTS #########
        save(csacqidx, file=paste(dataForRdir, "csacqidx.ridx", sep=""))
        save(alldata, file=paste(dataForRdir, "alldata.rdat", sep=""))
        save(rats, file=paste(dataForRdir, "rats.rdat", sep=""))
        save(idx, file=paste(dataForRdir, "idx.rdat", sep=""))
        save(DSrespAll, file=paste(dataForRCumulative, "DSrespAll.rdat", sep=""))
        save(NSrespAll, file=paste(dataForRCumulative, "NSrespAll.rdat", sep=""))
        save(DStimeToSpare, file=paste(dataForRCumulative, "DStimeToSpare.rdat", sep=""))
        save(NStimeToSpare, file=paste(dataForRCumulative, "NStimeToSpare.rdat", sep=""))
        save(DStaskAcc, file=paste(dataForRCumulative, "DStaskAcc.rdat", sep=""))
        save(NStaskAcc, file=paste(dataForRCumulative, "NStaskAcc.rdat", sep=""))
        save(DSlatency,  file=paste(dataForRCumulative, "DSlatency.rdat", sep=""))
        save(NSlatency,  file=paste(dataForRCumulative, "NSlatency.rdat", sep=""))
        save(ITIlatency,  file=paste(dataForRCumulative, "ITIlatency.rdat", sep=""))
}


funcdirect <- "E:/Dropbox/NMDA/R functions/"
save(MedPCextract, file=paste(funcdirect, "MedPCextract.r", sep=""))
