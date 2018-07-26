### FUNCTION THAT EXTRACTS AND SAVES THE FOLLOWING BEHAVIOR-RELATED DATA OBJECTS:
# csacqidx -> index of all sessions with summary indices
# alldata -> list of objects with detailed about relevant session events (each object corresponds to each session, the number of each object corresponds to the row of the same number in 'csacqidx')
# rats -> unique names of all rats
# idx -> list where each object corresponds to one rat (the object rats will tell you which rat). Inside each object, there is a vector with the indices of each one of that animal's sessions in csacqidx and alldata
# DSrespAll -> every DS trial of every session: 1=responded to the cue; 0=ignored the cue
# NSrespAll -> Same thing but of NSs
# DStaskAcc -> How much faster was the animal to enter the rec. after the cue vs. during the ITI (average of 100 random 10s windows
# NStaskAcc -> same but of NS
# DStimeToSpare -> 10-latency to enter after DS onset
# NStimeToSpare -> same but of NS


MedPCextract <- function(funcdirect, datafolder, dataForRdir, dataForRCumulative, cuelength, consumeRewWdw){
        
        load (paste(funcdirect,"CStask.Rfunc",sep=""))
        load (paste(funcdirect,"mpcextract_blockSingle.Rfunc",sep=""))
        load (paste(funcdirect,"mpcextractSingle.Rfunc",sep=""))
        
        mpcdatadir = datafolder
        allfiles = list.files(mpcdatadir)
        
        cuelength = cuelength #in sec
        
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
                
                correctTrials <- findInterval(toanalyze$rewarddelivery, toanalyze$CSpluscue)
                
                #In July 2017 I changed the task so that the cue extends 2s after correct entry, so calculating the CSplusEnd will be different after 07/01/2017
                CSplusEnd = vector()
                
                if(as.Date(substr(toanalyze$info[2], 13, 20), format='%m/%d/%y')<as.Date("07/01/17", format="%m/%d/%y")) #Before change
                {
                        for(j in 1:length(CSpluscue)){
                                if(sum(j==correctTrials)>0){
                                        CSplusEnd[j] = toanalyze$rewarddelivery[j==correctTrials]
                                } else {
                                        CSplusEnd[j] = toanalyze$CSpluscue[j] + cuelength
                                }
                        }
                } else { #after change
                        for(j in 1:length(CSpluscue)){
                                if(sum(j==correctTrials)>0){
                                        CSplusEnd[j] = toanalyze$rewarddelivery[j==correctTrials]+2
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
                
                
                # Cue-specificity: to what extent are post-cue entries under the control of the cue as opposed to just random? Compare latency of entry after the cue with latency of entry during a random 10s ITI interval
                
                set.seed(1111)
                consumeRewWdw <- consumeRewWdw #Exclude beginning of the ITI in which the animal might be consuming reward (Let's say it takes the animal 4.5s to consume reward)
                Iter <- 100
                
                #If the 1st trial happened after 5s (version of the task with first cue being a long CS+ after 5s), the first ITI lantecy value will generate a NA, but that's ok. Will deal with that later
                
                firstTr <- 1
                
                ITIlatency <- sapply(seq(firstTr, length(ITIlength)), function(x){
                        eligibleITI <- ITIlength[x]-consumeRewWdw-cuelength #Let's discard the first few seconds of the ITI bc the animal might be consuming reward and the last 10 seconds bc we need a window to fit before the end of the cue
                        if(eligibleITI<=0.1) {consumeRewWdw <- 1} #A handful of trials have an ITI that is just 15s or slightly less so it'd give me no eligible ITI. If that's the case, reduce the consumption rel. window to 1
                        ITIwdwStarts <- trialStart[x]+consumeRewWdw+runif(Iter, min=0, max=eligibleITI)
                        ITIwdwEnds <- ITIwdwStarts+cuelength
                        
                        ITIlat <- sapply(seq(1, length(ITIwdwStarts)), function(y){
                                ITIwdwEntries <- receptacleentries[receptacleentries>=ITIwdwStarts[y] & receptacleentries<=ITIwdwEnds[y]]
                                if(length(ITIwdwEntries)==0){ITIlat <- cuelength} #If no entries were made during that window, determine that latency=maximum latency on that ITI window (i.e. cue length)
                                else {
                                        firstEnt <- min(ITIwdwEntries) #Select first entry after the window in case there's more than one
                                        ITIlat <- firstEnt-ITIwdwStarts[y]}
                        })
                        
                        ITIlatency <- mean(ITIlat, na.rm=T)
                })
                
                #Go to cued entry latency values and replace NA values with the length of the cue
                CSplusLat[is.na(CSplusLat)] <- cuelength
                CSminusLat[is.na(CSminusLat)] <- cuelength
                
                #If 1st cue is long CS+ after 5s, get rid of that trial for this calculation (bc the ITI is too short, 5s)
                if(ITIlength[1]<=6){allCuesOrder2 <- c(NA, allCuesOrder[-1]); CSplusLat2 <- c(NA, CSplusLat[-1])} else {
                        allCuesOrder2 <- allCuesOrder; CSplusLat2 <- CSplusLat}
                
                #allCuesOrder2 <- allCuesOrder; CSplusLat2 <- CSplusLat
                CSplusTA <- ITIlatency[which(allCuesOrder==1)] - CSplusLat2
                CSminusTA <- ITIlatency[which(allCuesOrder==2)] - CSminusLat
                
                # This block of code is for sth I explored once that is not very interesting
                # #Number of responded to trials up to each trial
                # PreviouslyRespondedCSplusTrials <- sapply(seq(1, length(CSpluscue)), function(h){
                #         if(h==1){preRespCSplusTrials <- 0
                #         } else {preRespCSplusTrials <- sum(!is.na(CSplusresponse[1:(h-1)]))}
                # })
                # PreviouslyRespondedCSminusTrials <- sapply(seq(1, length(CSminuscue)), function(h){
                #         if(h==1){preRespCSminusTrials <- 0
                #         } else {preRespCSminusTrials <- sum(!is.na(CSminusresponse[1:(h-1)]))}
                # })
                # 
                # #For each reward, time from previous reward divided by cue-rew latency ('C/T ratio' as per Balsam, etc.)
                # CTratioPerTrial <- sapply(seq(1, length(CSpluscue)), function(h){
                #         NotNAtrialsIdx <- (1:length(CSplusresponse))[!is.na(CSplusresponse)]
                #         PreNotNAtrialsIdx <- NotNAtrialsIdx[NotNAtrialsIdx<h]
                #         # C
                #         if(length(PreNotNAtrialsIdx)==0){C <- NA} else {C <- CSplusresponse[h]-CSplusresponse[max(PreNotNAtrialsIdx)]}
                #         #T
                #         t <- CSplusresponse[h]-CSpluscue[h]
                #         
                #         CTratio <- C/t
                # })
                
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
        latencyCSplus <- c()
        latencyCSminus <- c()
        latencyCSplusSEM <- c()
        latencyCSminusSEM <- c()
        DSentryDur <- c()
        NSentryDur <- c()
        CSplusTaskAcc <- c()
        CSminusTaskAcc <- c()
        ITIPeriodentryDur <- c()
        ITIlatency <- c()
        PercShortDSent <- c()
        
        for (i in 1:length(alldata)){
                DSperc [[i]]<- length(alldata[[i]]$rewarddelivery)/length(alldata[[i]]$CSpluscue)
                NSperc [[i]]<- length(alldata[[i]]$CSminusresponse[which(alldata[[i]]$CSminusresponse>0)])/length(alldata[[i]]$CSminuscue)
                ITItotal[[i]] <- sum(alldata[[i]]$ITIlength)
                ITIperSec[[i]] <- round((length(alldata[[i]]$receptacleentries) - length(c(alldata[[i]]$CSplusresponse, alldata[[i]]$CSminusresponse)))/ITItotal[[i]], 2)
                if(ITIperSec[[i]]<0){ITIperSec[[i]]=0}
                latencyCSplus[[i]] <- mean(alldata[[i]]$CSplusLat, na.rm=T)
                latencyCSminus[[i]] <- mean(alldata[[i]]$CSminusLat, na.rm=T)
                latencyCSplusSEM[[i]] <-sd(alldata[[i]]$CSplusLat, na.rm=T)/sqrt(length(alldata[[i]]$CSplusLat[which(alldata[[i]]$CSplusLat >0)]))
                latencyCSminusSEM[[i]] <-sd(alldata[[i]]$CSminusLat, na.rm=T)/sqrt(length(alldata[[i]]$CSminusLat[which(alldata[[i]]$CSminusLat >0)]))
                ITIlatency[[i]]<- mean(alldata[[i]]$ITIlatency, na.rm=T)
                
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
        
        DSlatency <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        DSresp <- (!is.na(alldata[[sessIdx]]$CSplusresponse))*1
                        DSlat <- alldata[[sessIdx]]$CSplusLat
                        DSlat[which(is.na(DSlat))] <- cuelength
                        DSlat
                }))
        })
        
        DStimeToSpare <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        DSresp <- (!is.na(alldata[[sessIdx]]$CSplusresponse))*1
                        DSlat <- alldata[[sessIdx]]$CSplusLat
                        DSlat[which(is.na(DSlat))] <- cuelength
                        timeToSpare <- cuelength-DSlat
                }))
        })
        
        NStimeToSpare <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        NSresp <- (!is.na(alldata[[sessIdx]]$CSminusresponse))*1
                        NSlat <- alldata[[sessIdx]]$CSminusLat
                        NSlat[which(is.na(NSlat))] <- cuelength
                        timeToSpare <- cuelength-NSlat
                }))
        })
        
        NSlatency <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        NSresp <- (!is.na(alldata[[sessIdx]]$CSminusresponse))*1
                        NSlat <- alldata[[sessIdx]]$CSminusLat
                        NSlat[which(is.na(NSlat))] <- cuelength
                        NSlat
                }))
        })
       
        
        DStaskAcc <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        DStaskAcc <- alldata[[sessIdx]]$CSplusTA
                        DStaskAcc <- DStaskAcc[!is.na(DStaskAcc)]#Remove NAs bc it won't work otherwise (there shouldn't be any)
                }))
        })
        
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
        
        
        #Now make object of DStaskAcc without excluding the first trial of each session. Not for CP calculation
        DStaskAccAllTrials <- lapply(seq(1, length(rats)), function(k){
                do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
                        sessIdx <- idx[[k]][m]
                        DStaskAcc <- alldata[[sessIdx]]$CSplusTA
                        isnaIdx <- (1:length(DStaskAcc))[is.na(DStaskAcc)]  #Identify NA trials
                        replacement <- sapply(seq(1, length(isnaIdx)), function(l){
                                a <- DStaskAcc[isnaIdx[l]+1]
                        })
                        DStaskAcc[isnaIdx] <- replacement
                        DStaskAcc
                }))
        })
        
        
        # #CTratio <- lapply(seq(1, length(rats)), function(k){
        #         do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
        #                 sessIdx <- idx[[k]][m]
        #                 CTratio <- alldata[[sessIdx]]$CTratioPerTrial
        #                 CTratio[is.na(CTratio)] <- 0
        #                 CTratio
        #         }))
        # })
        # 
        # PreRespondedCSplus <- lapply(seq(1, length(rats)), function(k){
        #         
        #         cumsum(do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
        #                 sessIdx <- idx[[k]][m]
        #                 RespondedCSplusPerSess <- diff(alldata[[sessIdx]]$PreviouslyRespondedCSplusTrials)
        #         })
        #         )
        #         )
        # })
        # 
        # PreRespondedCSminus <- lapply(seq(1, length(rats)), function(k){
        #         
        #         cumsum(do.call("c", lapply(seq(1, length(idx[[k]])), function(m){
        #                 sessIdx <- idx[[k]][m]
        #                 RespondedCSminusPerSess <- diff(alldata[[sessIdx]]$PreviouslyRespondedCSminusTrials)
        #         })
        #         )
        #         )
        # })
        
        
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
        #save(CTratio, file=paste(dataForRdir, "CTratio.rdat", sep=""))
        #save(PreRespondedCSplus, file=paste(dataForRdir, "PreRespondedCSplus.rdat", sep=""))
        #save(PreRespondedCSminus, file=paste(dataForRdir, "PreRespondedCSminus.rdat", sep=""))
        #save(DStaskAccAllTrials, file=paste(dataForRdir, "DStaskAccAllTrials.rdat", sep=""))
}

funcdirect <- "E:/Dropbox/NMDA/R functions/"
save(MedPCextract, file=paste(funcdirect, "MedPCextract.r", sep=""))
