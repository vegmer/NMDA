compareCPs <- function(data, expNames, colindx=c("gray30", "gray45", "gray60", "gray75"), graphFolder, minSess=5, imgFormat="pdf", graphs=c("Change point", "Non learners", "Dynamic Interval", "Asymptote")){
        lapply(seq(1, length(graphs)), function(x){
               if(graphs[x]=="Change point"){
                       
                       if(imgFormat=="png"){filetitle=paste(graphFolder, graphs[x], expNames, ".png", sep=""); png(filename = filetitle)}
                       if(imgFormat=="pdf"){filetitle=paste(graphFolder, graphs[x], expNames, ".pdf", sep=""); pdf(file = filetitle)}
                       
                       plot.new()
                       CPcomp <- sapply(seq(1, length(data)), function(i){
                               minSessSel <- data[[i]][data[[i]]$nSess>=minSess,] #Select animals that got at least minSess sessions
                               minCPSessSel <- minSessSel[(is.na(minSessSel$CPsess) | minSessSel$CPsess<=minSess),] #Select animals whose CP either didn't happen at all or happened before minSess (if CP happened after minSess, it's not fair to say that animal had a CP)
                               dat <- minCPSessSel$CP
                               plot.window(xlim=c(0, length(data)+1), ylim=c(0, 250))
                               widthRect <- 0.5
                               widthNARect <- 0.4
                               rect(xleft=i-widthRect, xright=i+widthRect, ybottom=0, ytop=mean(dat, na.rm=T), col=colindx[i], border = "white")
                               jitter=rep(c(-0.08, 0.08, -0.06, 0.06, -0.04, 0.04),length(dat))
                               sapply(seq(1, length(dat)), function(n){
                                  points(x=i+jitter[n], y=dat[n], pch=19) #CP of animals that learned     
                               })
                               
                               axis(side=2, at=seq(0, 250, by=50), las=2)
                               mtext(side=1, at=i, expNames[i], cex=1.5, font=2)
                               mtext(side=2, text=graphs[x], line=2.5, cex=1.2, font=2)
                               mtext(side=3, text=paste("Change point (CP) trial in rats that had a CP before the end of session", minSess, sep=" "))
                               })
                       dev.off()
               } 
                
                if(graphs[x]=="Non learners"){
                        if(imgFormat=="png"){filetitle=paste(graphFolder, graphs[x], expNames, ".png", sep=""); png(filename = filetitle)}
                        if(imgFormat=="pdf"){filetitle=paste(graphFolder, graphs[x], expNames, ".pdf", sep=""); pdf(file = filetitle)}
                        
                        plot.new()
                        plot.window(xlim=c(0, length(data)+1), ylim=c(0, 100))
                        
                        NL <- sapply(seq(1, length(data)), function(i){
                        minSessSel <- data[[i]][data[[i]]$nSess>=minSess,] #Select animals that got at least minSess sessions
                        NumberNonLearners <- sum(is.na(minSessSel$CPsess)) + sum(minSessSel$CPsess>minSess, na.rm=T)
                        PercNonLearners <- (NumberNonLearners/nrow(minSessSel))*100
                        if(PercNonLearners==0){PercNonLearners=0.5} #If all animals learned, I want a tiny line to mark the bar, not just a white space
                        widthRect <- 0.5
                        
                        rect(xleft=i-widthRect, xright=i+widthRect, ybottom=0, ytop=PercNonLearners, col=colindx[i], border="white")
                        mtext(side=3, text=paste("% Animals with no change point before the end of session", minSess, sep=" "))
                        axis(side=2, las=2)
                        mtext(side=2, line=2.5, text="Percentage", cex=1.5, font=2)
                        axis(side=1, at=i, expNames[i], tick=F, cex=1.5, font=2)
                        axis(side=1, line=1.5, at=i, tick=F, labels=paste("(N=", nrow(minSessSel), ')'))
                        })
                  dev.off()
                }
                
                if(graphs[x]=="Dynamic Interval"){
                        
                        if(imgFormat=="png"){filetitle=paste(graphFolder, graphs[x], expNames, ".png", sep=""); png(filename = filetitle)}
                        if(imgFormat=="pdf"){filetitle=paste(graphFolder, graphs[x], expNames, ".pdf", sep=""); pdf(file = filetitle)}
                        
                        plot.new()
                        DIcomp <- sapply(seq(1, length(data)), function(i){
                                minSessSel <- data[[i]][data[[i]]$nSess>=minSess,] #Select animals that got at least minSess sessions
                                minCPSessSel <- minSessSel[(is.na(minSessSel$CPsess) | minSessSel$CPsess<=minSess),] #Select animals whose CP either didn't happen at all or happened before minSess (if CP happened after minSess, it's not fair to say that animal had a CP)
                                dat <- minCPSessSel$DynamicInterv
                                plot.window(xlim=c(0, length(data)+1), ylim=c(0, 100))
                                widthRect <- 0.5
                                jitter=rep(c(-0.08, 0.08, -0.06, 0.06, -0.04, 0.04),length(dat))
                                ytop=mean(dat, na.rm=T)
                                if(ytop==0){ytop=0.5}
                                rect(xleft=i-widthRect, xright=i+widthRect, ybottom=0, ytop=ytop, col=colindx[i], border = "white")
                                sapply(seq(1, length(dat)), function(n){
                                        points(x=i+jitter[n], y=dat[n], pch=19) #CP of animals that learned     
                                
                                axis(side=1, at=i, expNames[i], tick=F, cex=1.5, font=2)
                                axis(side=2, las=2)
                                mtext(side=2, line=2.5, text="Dynamic interval (#trials)")
                                mtext(side=3, text=paste("Abruptness of behavioral change in animals with CP before end of session", minSess))
                                })
                                
                        })
                        dev.off()
                } 
                
                if(graphs[x]=="Asymptote"){
                        
                        if(imgFormat=="png"){filetitle=paste(graphFolder, graphs[x], expNames, ".png", sep=""); png(filename = filetitle)}
                        if(imgFormat=="pdf"){filetitle=paste(graphFolder, graphs[x], expNames, ".pdf", sep=""); pdf(file = filetitle)}
                        
                        plot.new()
                        Asymptcomp <- sapply(seq(1, length(data)), function(i){
                                
                                minSessSel <- data[[i]][data[[i]]$nSess>=minSess,] #Select animals that got at least minSess sessions
                                minCPSessSel <- minSessSel[(is.na(minSessSel$CPsess) | minSessSel$CPsess<=minSess),] #Select animals whose CP either didn't happen at all or happened before minSess (if CP happened after minSess, it's not fair to say that animal had a CP)
                                onlyCPSel <- minCPSessSel[!is.na(minCPSessSel$CP),]
                                slopePre <- onlyCPSel$slopePre
                                slopePost <- onlyCPSel$slopePost
                                
                                plot.window(xlim=c(0, length(data)+1), ylim=c(-2, 5))
                                widthRect <- 0.4
                                rect(xleft=i-widthRect, xright=i, ybottom=0, ytop=mean(slopePre, na.rm=T), col=colindx[i], border = "white")
                                rect(xleft=i, xright=i+widthRect, ybottom=0, ytop=mean(slopePost, na.rm=T), col=colindx[i], border = "white")
                                
                                #points(x=rep(i-widthRect+widthRect/2, length(slopePre)), y=slopePre, pch=19, cex=0.7)
                                #points(x=rep(i+widthRect-widthRect/2, length(slopePost)), y=slopePost, pch=19, cex=0.7)
                            
                                sapply(seq(1, length(slopePre)), function(n){
                                        lines(x=c(i-widthRect+widthRect/2, i+widthRect-widthRect/2), y=c(slopePre[n], slopePost[n])) #CP of animals that learned     
                                
                                
                                axis(side=1, at=i, tick=F, labels=expNames[i], cex=1.5, font=2)
                                
                                text("Pre", x=i-widthRect/2, y=5, font=2)
                                text("Post", x=i+widthRect/2, y=5, font=2)
                                })
                                
                                axis(side=2, las=2)
                                mtext(side=2, text="ITI latency-Cued latency (s)", line=2.5, font=2)
                                
                                mtext(side=3, line=2, text="Performance pre vs. post CP in animals with CP before end of session 5")
                                
                        })
                        dev.off()
                } 
        })
}


funcdirect <- "E:/Dropbox/DISSERTATION/R functions/"
save(compareCPs, file=paste(funcdirect, "compareCPs.r", sep=""))
