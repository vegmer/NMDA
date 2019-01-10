#### HEATMAPS FROM FIRINGMAT

#Load data
setwd("E:/Dropbox/DISSERTATION/")
setwd(paste(getwd(), "/Dropbox/DISSERTATION/", sep=""))
datadir <- paste(getwd(), "/Experiment 4/Data for R/", sep="")

load(file=paste(datadir, "DSentryALL.rdat", sep="")) #All data using REWARDED ENTRY as event with minpsth=2 and maxpsth=3 and bin=50ms. Object name might be allNeuronsDS
load(file=paste(datadir, "ITIentryALL.rdat", sep="")) #All data using ITI ENTRY as event with minpsth=2 and maxpsth=3 and bin=50ms. Object name might be allNeuronsNS
load(file=paste(datadir, "ZscDSentry.rdat", sep="")) #Firingmat in Zsc using REWARDED ENTRY as event with minpsth=2 and maxpsth=3 and bin=50ms. 2s pre event BL Object name might be ZDS
load(file=paste(datadir, "ZscITIentry.rdat", sep="")) #Firingmat in Zsc using ITI ENTRY as event with minpsth=2 and maxpsth=3 and bin=50ms. 2s pre event BL Object name might be ZNS


graphFolder <- paste(getwd(),"/Graphs/", sep="")
neurographFolder <- paste(graphFolder, "/Neuronal data/", sep="")
behgraphFolder <- paste(graphFolder,"/Behavioral data/", sep="")
mixedgraphFolder <- paste(graphFolder,"/Mixed graphs/", sep="")


#### HEATMAP WITH OVERALL EVOLUTION OF PERFORMANCE TOO (DS RESPONSE RATIO)
library(dplyr)
library("colorRamps")

#Arrange toplot based on performance (DSRR), then animal and and then firing rate to DS ENTRY (ZDS)
toplot2 <- arrange(toplot, as.numeric(DSRR), desc(as.numeric(rat)), as.numeric(NSRR), as.numeric(ZDS))
toplot2[,(ncol(toplot2)-3):ncol(toplot2)] <- round(toplot2[,(ncol(toplot2)-3):ncol(toplot2)], 1)

#I am going to use Firingmat in Zsc (ZDS) as my dataset. I want to rearrange the columns (units), based on performance, I can use the toplot2$Idx
ZscData <- ZDS
UnitOrder <- toplot2$Idx

toplot3 <- c()
for(i in 1:length(UnitOrder)){
        unit <- ZscData[,UnitOrder[i]]
        toplot3 <- cbind(toplot3, unit)
}
colnames(toplot3)<-UnitOrder
ZscData <- toplot3


# Set limits of FR
boxplot(ZscData) #Use this boxplot to set the max limit of FR so that outliers don't skew the distribution and therefore obscure the relationship
minFR <- -4
maxFR <- 4 #Anything above this will be categorized as this value

colValues <- seq(minFR, maxFR, by=0.1)
colValues <- round(colValues, 1)
mypalette <- colorRampPalette(colors=c("darkblue", "orangered"))(length(colValues))
NAcolor <- "white" #Color for NA cells

mypalette <- c(NAcolor, mypalette)
idx <- 1:ncol(ZscData)


############### PLOT ###################################
plot.new()
plot.window(xlim=c(1, 5.6), ylim=c(0, length(idx)))

xmin <- 1.20 #Choose the x position of the beginning of the heatmap (reference for everything else)
xmax <- 3
nbins <- nrow(ZscData)
xlocations <- seq(xmin, xmax, by=(xmax-xmin)/nbins)
#This loop makes the actual grid with colors
for(i in idx){
        for(j in 1:nrow(ZscData)){
                a <- xlocations[j]
                b <- xlocations[j+1]
                
                #Find Zsc values in the idx of all possible values. 
                #For some reason treating them as numbers gives me wacky results
                FRdatapoint <- round(ZscData[j,i], 1)
                if(is.na(FRdatapoint)) {k <- 1} else { #First color in mycolorpalette is white, so if value is NA, make it white by assigning k=1
                        if(FRdatapoint > maxFR){k <- length(mypalette)} #If FR value is bigger than maxFR (outlier), assign the last value in the palette
                        if(FRdatapoint < maxFR){k <- grep(paste("^", FRdatapoint, "$", sep=""), colValues)+1} #Plus one because the lowest value in "mypalette" is white (for NAs) 
                        if(FRdatapoint < minFR){k <- 2} #If FR value is smaller than minFR (outlier), assign the lowest value in the palette (that is not white, which is 1)
                }
                
                #if(!is.na(k)==TRUE){
                        colpick <- mypalette[k]
                        #} else {colpick <- "white"}
                
                rect(xleft=a, 
                     xright=b, 
                     ytop=i, 
                     ybottom=i-1, 
                     col=colpick, 
                     border=NA)    
        }
} 


#Neuron # axis
axis(side=2, pos=xmin-0.02, cex.axis=0.6, at=c(0, max(idx)-1), labels=c(min(idx), max(idx)), las=3, padj=1)
mtext(text="Neuron #", side=2, line=0.6, at=length(idx)/2, cex=0.8, font=2, padj=5.5)

#Legend
legendTop <- max(idx)
legendBottom <- min(idx)
legendStep <-(legendTop-legendBottom)/length(colValues)

for (i in 1:length(colValues)){
        if(i==min(1:length(colValues))){ybottom=legendBottom; ytop=legendBottom+legendStep} 
        else {ybottom=legendBottom+(legendStep*(i-1));ytop=legendBottom+(legendStep*(i-1))+legendStep}
        
        rect(xleft=xmin-0.3, xright=xmin-0.25,
             ybottom=ybottom, ytop=ytop,
             col=mypalette[i],
             border=NA)
}

#Legend axis
axis(side=2, pos=xmin-0.31, at=seq(legendBottom+legendStep, legendTop, by=5), 
     labels=round(seq(from=minFR, to=maxFR, length.out=length(seq(legendBottom+legendStep, legendTop, by=5))), 1), las=2, padj=-0.05)

mtext("Firing rate(Z sc)", side=3, line=-1, at=xmin-0.4, cex=0.8, font=2)
mtext(">=", side = 3, line=-1.9, at=xmin-0.7)
mtext("<=", side = 1, line=-2.7, at=xmin-0.7)

# Main title
mtext(text="FR around time of REWARDED ENTRY throughout learning (NO INFUSION)", side=3, line=2, cex=1.2, font=2, col="gray50")

#X axis
secs <- seq(-psthmin, psthmax, by=1) #psthmin and max are objects I created in the MV neural histograms TRAINING script
divs <- seq(xmin, xmax, by=(xmax-xmin)/length(secs))
axis(side=1, line=-1.2, at=divs, labels=c(secs, max(secs)+1), padj=-1, cex=0.8)
mtext("Time from ENTRY", side=1, line=0, at=xmin+(xmax-xmin)/2, font=2)

ZeroPoint <- divs[which(secs==0)]
abline(v=ZeroPoint, col="black")
mtext("NAcC FIRING RATE", side=3, line=-0.5, at=xmin+(xmax-xmin)/2, font=2, cex = 0.9)



### DIVISIONS BY RAT AND DS RESPONSE RATIO
ratsess <- c()
for(i in 1:nrow(toplot2)){
        ratsess <- c(ratsess, as.numeric(paste(toplot2$rat[i], toplot2$session[i], sep="")))
}

PerfValuesbyRAT <- rle(ratsess)
PerfLinesbyRAT <- cumsum(PerfValuesbyRAT$lengths)
for(i in 1:length(PerfLinesbyRAT)){
        segments(y0=PerfLinesbyRAT[i], y1=PerfLinesbyRAT[i],
                 x0=xmin, x1=xmin+4.3, col="gray")
}
segments(y0=0, y1=0, x0=xmin, x1=xmin+2.6, col="gray")

#Make a rectangle where the behavior is going to be plotted
rectMin <- xmin+1.9
rectMax <- xmin+3
rect(xleft=rectMin, xright=rectMax, ybottom=0, ytop=max(idx))

#Here I'm plotting changes in DSRR. Min possible is 0 and max is 1. Make 0.1 bins.
minBeh <- 0
maxBeh <- 1
binsBeh <- seq(minBeh, maxBeh, by=0.1)

sizeBin <- (rectMax-rectMin)/(length(binsBeh)-1) #Size of each bin on the plot (minus one because I don't want to count in bin 0)
coordBins <- c()
for(i in 1:length(binsBeh)){
        a <- rectMin + i*sizeBin
        coordBins <- c(coordBins, a)
}

#Make vertical lines with these bins
for(i in 0:(length(coordBins)-1)){
        if(i==0){x <- rectMin} else {x <- coordBins[i]}
        segments(x0=x, x1=x, y0=min(idx)-1.5, y1=max(idx), col="gray90")
        if(i==(length(coordBins)-1)/2){
                segments(x0=x, x1=x, y0=min(idx)-1.5, y1=max(idx), col="gray")
        }
}

rect(xleft=rectMin, xright=rectMax, ybottom=0, ytop=max(idx)) #redo rectangle bc it's obscured

#Tickmarks w/o taking into account different rats with same performance
#axis(side=4, pos=rectMax, at=c(0, cumsum(PerfValues$lengths)), labels=FALSE)

#Tickmarks taking into account different rats with same performance
#axis(side=4, pos=rectMax, at=c(0, cumsum(PerfValuesbyRAT$lengths)), labels=FALSE)

#Labels (so that they appear between tickmarks) NOT BY RAT, JUST PERFORMANCE
# midTicks <- c()
# for(i in 1:length(PerfLines)){
#   if(i == 1){a <- (PerfLines[i]-0)/2} 
#   else{a <- PerfLines[i-1]+(PerfLines[i]-PerfLines[i-1])/2}
#   midTicks <- c(midTicks, a)
# }
# axis(side=4, pos=rectMax, at=midTicks, tick=F, labels=length(PerfLines):1, cex.axis=0.6, padj=-3)
# mtext("Ranking of session based on performance (DS response ratio)", side=4, line=-11, cex.axis=0.7, font=2)

# #Labels to show in between tickmarks taking into account each rat
midTicks <- c()
for(i in 1:length(PerfLinesbyRAT)){
        if(i == 1){a <- (PerfLinesbyRAT[i]-0)/2}
        else{a <- PerfLinesbyRAT[i-1]+(PerfLinesbyRAT[i]-PerfLinesbyRAT[i-1])/2}
        midTicks <- c(midTicks, a)
}
# axis(side=4, pos=rectMax, at=midTicks, tick=F, labels=length(PerfLinesbyRAT):1, cex.axis=0.6, padj=-3)
# mtext("Ranking of session based on performance (DS response ratio)", side=4, line=-11, cex.axis=0.7, font=2)



# #Determine the length of the DSRR bars
# PerfPerSessPERCENT <- round(PerfValues$values, 2)*100
# lengthPerPoint <- (rectMax-rectMin)/100 #Determine length per DSRR point (0.01)

#Determine the length of the DSRR bars BY RAT
PerfPerSessPERCENT <- round(toplot2$DSRR[PerfLinesbyRAT], 2)*100
NSPerSessPERCENT <- round(toplot2$NSRR[PerfLinesbyRAT], 2)*100
lengthPerPoint <- (rectMax-rectMin)/100 #Determine length per DSRR point (0.01)


#Draw the bars of DSRR performance
offset<- 0.4 #y distance from midtickmark of the bars
b <- c()
for(i in 1:length(PerfPerSessPERCENT)){
        b <- PerfPerSessPERCENT[i]*lengthPerPoint
        segments(x0=rectMin, x1=rectMin+b, 
                 y0=midTicks[i]+offset, y1=midTicks[i]+offset, col="darkgoldenrod4", lwd=2)
}

#Draw the bars of NSRR performance
b <- c()
for(i in 1:length(NSPerSessPERCENT)){
        b <- NSPerSessPERCENT[i]*lengthPerPoint
        segments(x0=rectMin, x1=rectMin+b, 
                 y0=midTicks[i]-offset, y1=midTicks[i]-offset, col="darkgoldenrod2", lwd=2)
}

#Legend
legend(x=rectMin+0.45, y=20.5, legend=c("S+", "S-"), fill=c("darkgoldenrod4", "darkgoldenrod2"))

# X Axis of beh plot
axis(side=1, line=-1.2, at=c(rectMin, rectMin+(rectMax-rectMin)/2, rectMax), labels=c(minBeh, (maxBeh-minBeh)/2, maxBeh), padj=-1, cex=0.8)
mtext("PERFORMANCE", side=3, line=-0.5, at=rectMin+(rectMax-rectMin)/2, font=2, cex=0.9)
#mtext("(Response ratio)", line=-1.3, at=rectMin+(rectMax-rectMin)/2, font=3, cex=0.8)
mtext("Response ratio", side=1, line=0, at=rectMin+(rectMax-rectMin)/2, font=2)

### LAST RECTANGLE WITH DISTRIBUTION OF SESSION BY ANIMAL
#Make a rectangle where the animals are going to be plotted
rectMin2 <- xmin+3.1
rectMax2 <- xmin+4.3
rect(xleft=rectMin2, xright=rectMax2, ybottom=0, ytop=max(idx))

#MAke columns for each animal
ridx <- sort(unique(toplot2$rat))
minrat <- min(ridx)
maxrat <- max(ridx)
#ratCols <- seq(minrat, maxrat, by=1)
#ratCols <- c(minrat, maxrat)
ratCols <- unique(ridx)

sizeCol <- (rectMax2-rectMin2)/(length(ratCols)) #Size of each rat column on the plot (minus one because I don't want to count in column 0)
coordCols <- c()
for(i in 1:length(ratCols)){
        a <- rectMin2 + i*sizeCol
        coordCols <- c(coordCols, a)
}

mtext("DATA SOURCE", cex=0.9, side=3, line=-0.5, at=rectMin2+(rectMax2-rectMin2)/2, font=2)
mtext("(Session # by animal)", line=-1.3, at=rectMin2+(rectMax2-rectMin2)/2, font=3, cex=0.8)

#Make vertical lines with these columns
for(i in 0:(length(coordCols)-1)){
        if(i==0){x <- rectMin2} else {x <- coordCols[i]}
        segments(x0=x, x1=x, y0=min(idx)-1.5, y1=max(idx))
        if(i==(length(coordCols)-1)/2){
                segments(x0=x, x1=x, y0=min(idx)-1.5, y1=max(idx))
        }
}

## POPULATE THESE COLUMNS
coordCols2 <- c(rectMin2, coordCols[-length(coordCols)])
for(i in unique(toplot2$rat)){
        byrat <- filter(toplot2, rat==i)
        for(j in byrat$session){
                ybottom=min(idx[(toplot2$rat==i)&(toplot2$session==j)])-1
                ytop=max(idx[(toplot2$rat==i)&(toplot2$session==j)])
                
                whichRats <- unique(toplot2$rat)
                
                rect(xleft=coordCols2[whichRats==i], xright=coordCols2[whichRats==i]+sizeCol, 
                     ybottom=ybottom,
                     ytop=ytop,
                     col="darkred")
                text(x=coordCols2[whichRats==i]+sizeCol/2, y=ybottom+(ytop-ybottom)/2, labels=j, col="pink", srt=90, cex=1, font=2) 
        }
}

axis(side=1, line=-2.2, at=coordCols2+sizeCol/2, labels=whichRats, padj=0.5, tick=F, las=2, font=3, srt=270, cex.axis=0.8)
mtext("Animal ID", side=1, line=0, at=rectMin2+(rectMax2-rectMin2)/2, font=2)




