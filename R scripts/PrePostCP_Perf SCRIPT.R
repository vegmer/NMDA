PrePostCP_Perf <- function(data, CP, y_axis_label, graphFolder=PerfRelToCPFolder, plot=T){
        
        #If I'm looking into ITI latency or ITI response ratio, then subset the ITIs that preceded S+s
        if(grepl("ITI", y_axis_label)){
                data <- lapply(seq(1, length(rats)), function(m){
                        sessIdx <- idx[[m]]
                        cueKind <- do.call("c", lapply(seq(1, length(sessIdx)), function(n){
                                alldata[[sessIdx[n]]]$orderCues
                        }))
                        
                        data[[m]][cueKind==1]
                })
                
                
        }
        
        perfAll <- t(sapply(seq(1, length(data)), function(x){
                dataSel <- data[[x]]
                if(is.na(CP[x])){PerfPre <- NA; PerfPost <- NA} 
                else {
                        PerfPre <- mean(dataSel[1:CP[x]], na.rm=T)
                        PerfPost <- mean(dataSel[CP[x]:length(dataSel)], na.rm=T)
                }
                data.frame(PerfPre=PerfPre, PerfPost=PerfPost)
        }))
        
        perfAll <- as.data.frame(perfAll)
        
        if(plot==T){
                
                imagename <- paste(graphFolder, "Pre vs Post ", y_axis_label, " Average.pdf",  sep="")
                
                pdf(file = imagename, width = 5, height=8)
                
                plot.new()
                ymin=floor(min(c(unlist(perfAll$PerfPre), unlist(perfAll$PerfPost)), na.rm=T))
                ymax=ceiling(max(c(unlist(perfAll$PerfPost), unlist(perfAll$PerfPre)), na.rm=T))
                
                plot.window(xlim=c(0, 0.5), ylim=c(ymin, ymax+0.2))
                
                rect(xleft=0, xright=0.25, ybottom=0, ytop=mean(unlist(perfAll$PerfPre), na.rm=T), border = "white", col="gray60")
                rect(xleft=0.25, xright=0.5, ybottom=0, ytop=mean(unlist(perfAll$PerfPost), na.rm=T), border = "white", col="gray60")
                
                sapply(seq(1, nrow(perfAll)), function(x){lines(x=c(0.125, 0.375), y=c(unlist(perfAll$PerfPre)[x], unlist(perfAll$PerfPost)[x]), lwd=2, col="gray20")})
                
                abline(h=0, lty=3)
                axis(side=1, at=c(0.125, 0.375), labels=c("Before CP", "After CP"), cex.axis=1.4, tick = FALSE)
                axis(side=2, las=2, cex.axis=1.4)
                
                mtext(side=2, line=3, text=y_axis_label, cex=1.5, font=2)  
                
                # Test significance of pre vs. post difference
                diff <- unlist(perfAll$PerfPre)-unlist(perfAll$PerfPost)
                normalityTest <- shapiro.test(diff)
                #IF the data (difference btwn pre and post) follows a normal distribution:
                if(normalityTest$p.value>=0.05){statTest <- t.test(x=unlist(perfAll$PerfPre), y=unlist(perfAll$PerfPost), paired = TRUE)}
                
                #If it doesn't:
                if(normalityTest$p.value<0.05){statTest <- wilcox.test(unlist(perfAll$PerfPre), unlist(perfAll$PerfPost), paired = TRUE)}
                
                p <- statTest$p.value
                
                segments(x0=0.125, x1=0.375, y0=ymax+0.1, y1=ymax+0.1, lwd=4)
                text(x=0.25, y=ymax+0.2, labels = giveStars(p), cex=1.5, font=2)
                
                
                dev.off()
                
                statTest
        }
        
        
}

save(PrePostCP_Perf, file="E:/Dropbox/NMDA/R functions/PrePostCP_Perf.R")
save(PrePostCP_Perf, file="E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R functions/PrePostCP_Perf.R")

