
FRandPerf_Scatterplot <- function(condition="All", data=toplotSumm, graphFolder=ScatterplotFolder, CPtoo=TRUE, CPmerged=FALSE,
                                  color=c(colindxC[1], colindxB[1], colindx[1]), behIndex="PerfIndex", chunk=NA){
        
        
        if(CPtoo==TRUE){CPincl="CPincl"} else {CPinc="No CP"}
        if(CPmerged==TRUE){CPmerg="CPmerged"} else {CPmerg="CP apart"}
        
        filename=paste(graphFolder, "FR vs Perf Scatterplot", condition, CPincl, CPmerg, behIndex, ".pdf", sep=" ")
        
        pdf(file=filename)
        
        plot.new()
        
        #For some reason these columns are in factor mode which is a headache, convert to double
        data$CSplusTA <- as.numeric(as.character(data$CSplusTA))
        data$CSplusLat <- as.numeric(as.character(data$CSplusLat))
        data$DSRR <- as.numeric(as.character(data$DSRR))
        
        if(!is.na(chunk)){
                data <- data[data$BeforeCP==chunk, ]
        }
        
        if(behIndex=="PerfIndex"){
                minX <- floor(min(data$CSplusTA, na.rm=T))
                maxX <- ceiling(max(data$CSplusTA, na.rm=T))
                labpos <- c(-1, 1)
                xAxTitle <- "Performance Index"
                perfAll <- data$CSplusTA
        }
        
        
        if(behIndex=="Latency"){
                minX <- 0; maxX <- 10
                labpos <- c(2, 5); xAxTitle <- "S+ latency (s)"
                perfAll <- data$CSplusLat
        }
        
        if(behIndex=="RespRatio"){
                minX <- 0; maxX <- 1
                labpos <- c(0.2, 0.7)
                xAxTitle <- "Response ratio"
                perfAll <- data$DSRR
                
        }
        
        minY <- floor(min(data$ZDS, na.rm=T))
        maxY <- 9
        
        
        plot.window(xlim=c(minX, maxX), ylim=c(minY, maxY))
        
        SessIdx <- seq(-1, 1, 1)
        
        if(CPtoo==FALSE){
                CPsessIndex <- data$BeforeCP==0
                data <- data[CPsessIndex==F, ]
                SessIdx <- c(-1, 1)
        }
        
        sapply(seq(1, length(SessIdx)), function(x){
                
                dataSel <- data[data$BeforeCP==SessIdx[x],]
                
                if(SessIdx[x]==1){pch=19; col=color[1]}
                
                if(CPmerged==TRUE){
                        if(SessIdx[x]<1){pch=19; col=color[2]}
                } else {
                        if(SessIdx[x]==0){pch=19; col=color[2]}
                        if(SessIdx[x]==-1){pch=19; col=color[3]}
                }
                
                if(behIndex=="PerfIndex"){perf <- dataSel$CSplusTA}
                if(behIndex=="Latency"){perf <- dataSel$CSplusLat}
                if(behIndex=="RespRatio"){perf <- dataSel$DSRR}
                
                FR <- as.numeric(as.character(dataSel$ZDS))
                
                points(x=perf, y=FR, pch=pch, cex=2, col=col)
                
                axis(side=1, cex.axis=1.4)
                axis(side=2, at=seq(-3, 9, by=1), cex.axis=1.4, las=2)
                
                abline(h=0, lty=2)
                abline(v=0)
                
                mtext(side=2, line=2.5, text="Firing rate (Z sc.)", cex=1.5, font=2)
                mtext(side=1, line=2.5, text=xAxTitle, cex=1.5, font=2)
                
        })
        
        
        FRAll <- as.numeric(as.character(data$ZDS))
        
        fit <- lm(FRAll ~ perfAll)
        Rsq <- summary(fit)$r.squared
        coeff <- summary(fit)$coefficients[2,1]
        p.val <- summary(fit)$coefficients[2,4]
        r <- cor.test(FRAll, perfAll)$estimate
        r.ttest <- round(cor.test(FRAll, perfAll)$statistic, 2)
        r.df <- cor.test(FRAll, perfAll)$parameter
        r.pval <- round(cor.test(FRAll, perfAll)$p.value, 5)
        
        
        text(x=labpos[1], y=6, labels = "With outliers", font=2, col="gray50")
        text(x=labpos[1], y=5.5, labels = paste("Coeff. = ", round(coeff, 2), sep=""))
        text(x=labpos[1], y=5, labels = paste("p = ", round(p.val, 5), sep=""))
        text(x=labpos[1], y=4.5, labels = paste("R.sq = ", round(Rsq, 2), sep=""))
        text(x=labpos[1], y=4, labels= paste("r = ", round(r, 2), sep=""))
        text(x=labpos[1], y=3.5, labels=paste(paste("t = ", r.ttest, sep=""),
                                              paste("df= ", r.df, sep=""),
                                              paste("r_pval = ", r.pval, sep=""), sep="; "))
        
        
        abline(a=summary(fit)$coefficients[1, 1], b=summary(fit)$coefficients[2,1], lwd=2, col="gray50")
        
        #Find outliers
        cooksd <- cooks.distance(fit)
        outliers <- as.numeric(names(cooksd)[(cooksd > 3*mean(cooksd, na.rm=T))])
        sapply(seq(1, length(outliers)), function(j){
                text(x=perfAll[outliers[j]], y=FRAll[outliers[j]]+0.2, labels=outliers[j])
        })
        
        
        fit <- lm(FRAll[-outliers] ~ perfAll[-outliers])
        Rsq <- summary(fit)$r.squared
        coeff <- summary(fit)$coefficients[2,1]
        p.val <- summary(fit)$coefficients[2,4]
        r <- cor.test(FRAll[-outliers], perfAll[-outliers])$estimate
        r.ttest <- round(cor.test(FRAll[-outliers], perfAll[-outliers])$statistic, 2)
        r.df <- cor.test(FRAll[-outliers], perfAll[-outliers])$parameter
        r.pval <- round(cor.test(FRAll[-outliers], perfAll[-outliers])$p.value, 5)
        
        
        text(x=labpos[2], y=6, labels = "Without outliers", font=2, col="black")
        text(x=labpos[2], y=5.5, labels = paste("Coeff. = ", round(coeff, 2), sep=""))
        text(x=labpos[2], y=5, labels = paste("p = ", round(p.val, 5), sep=""))
        text(x=labpos[2], y=4.5, labels = paste("R.sq = ", round(Rsq, 2), sep=""))
        text(x=labpos[2], y=4, labels= paste("r = ", round(r, 2), sep=""))
        text(x=labpos[2], y=3.5, labels=paste(paste("t = ", r.ttest, sep=""),
                                              paste("df= ", r.df, sep=""),
                                              paste("r_pval = ", r.pval, sep=""), sep="; "))
        
        abline(a=summary(fit)$coefficients[1, 1], b=summary(fit)$coefficients[2,1], lwd=2)
        
        forLegend <- c("After CP", "CP", "Before CP"); colorLeg <- color
        if(!is.na(chunk)){m <- c(1, 0, -1); forLegend <- forLegend[chunk==m]; colorLeg <- color[chunk==m]}
        
        legend("bottomright", col=colorLeg, pch=19, legend=forLegend)
        
        dev.off()
        
        
        #Plot outlier test
        
        filename=paste(graphFolder, "OUTLIER TEST", condition, CPincl, CPmerg, behIndex, ".pdf", sep=" ")
        
        pdf(file = filename)
        
        plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
        abline(h = 3*mean(cooksd, na.rm=T), col="red")  # add cutoff line. The cutoff point is 3 times the mean of cook´s distances
        text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>3*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
        
        dev.off()
        
}


save(FRandPerf_Scatterplot, file="E:/Dropbox/NMDA/R Functions/FRandPerf_Scatterplot.R")

