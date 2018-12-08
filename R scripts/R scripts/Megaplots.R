### MAKING HEATMAPS ORGANIZED BY PERFORMANCE IN THE SESSION
setwd("C:/Users/mvega/Dropbox/DISSERTATION/")
setwd("C:/Users/Mercedes/Dropbox/DISSERTATION/")
load(file=paste(getwd(), "/Experiment 4/Data for R/allFRandBEHinZsc.rdat", sep=""))

toplot <- behANDneu #I wasn't able to load the object, so I'll just use toplot which was the original object

graphFolder <- paste(getwd(),"/Experiment 4/Graphs/", sep="")
neurographFolder <- paste(getwd(),"/Experiment 4/Graphs/Neuronal data/", sep="")
behgraphFolder <- paste(getwd(),"/Experiment 4/Graphs/Behavioral data/", sep="")
mixedgraphFolder <- paste(getwd(),"/Experiment 4/Graphs/Mixed graphs/", sep="")


#### MEGAPLOT WITH OVERALL EVOLUTION OF PERFORMANCE TOO (DS RESPONSE RATIO)
library(dplyr)
library("colorRamps")

#Arrange toplot based on performance (DSRR), then animal and and then firing rate to DS (ZDS)
toplot2 <- arrange(toplot, as.numeric(DSRR), desc(as.numeric(rat)), as.numeric(NSRR), as.numeric(ZDSresponded))
toplot2[,(ncol(toplot2)-3):ncol(toplot2)] <- round(toplot2[,(ncol(toplot2)-3):ncol(toplot2)], 1)
ZscData <- toplot2[,(ncol(toplot2)-3):ncol(toplot2)] #subset firing data

boxplot(ZscData) #Use this boxplot to set the max limit of FR so that outliers don't skew the distribution and therefore obscure the relationship
minFR <- -1.3
maxFR <- 5.5 #Anything above this will be categorized as this value

colValues <- seq(minFR, maxFR, by=0.1)
colValues <- round(colValues, 1)
mypalette <- colorRampPalette(colors=c("darkblue", "orangered"))(length(colValues))
NAcolor <- "white" #Color for NA cells

mypalette <- c(NAcolor, mypalette)
idx <- 1:nrow(toplot2)

plot.new()
plot.window(xlim=c(1, 4), ylim=c(1, length(idx)))

xmin <- 1.20 #Choose the x position of the beginning of the heatmap (reference for everything else)

#This loop makes the actual grid with colors
for(i in idx){
  for(j in 1:ncol(ZscData)){
    if(j==1){a <- xmin}
    if(j==2){a <- xmin+0.25}
    if(j==3){a <- xmin+0.51}
    if(j==4){a <- xmin+0.76}
    
    #Find Zsc values in the idx of all possible values. 
    #For some reason treating them as numbers gives me wacky results
    FRdatapoint <- ZscData[i,j] 
    if(is.na(FRdatapoint)) {k <- 1} else { #First color in mycolorpalette is white, so if value is NA, make it white by assigning k=1
      if(FRdatapoint > maxFR){k <- length(colValues)} #If FR value is bigger than maxFR (outlier), assign the last value in the palette
      if(FRdatapoint < maxFR){k <- grep(paste("^", FRdatapoint, "$", sep=""), colValues)+1} #Plus one because the lowest value in "mypalette" is white (for NAs)    
    }
    
    if(!is.na(k)==TRUE){colpick <- mypalette[k]} else {colpick <- "white"}
    
    rect(xleft=a, 
         xright=a+0.25, 
         ytop=i, 
         ybottom=i-1, 
         col=colpick, 
         border=NA)    
  }
} 

#Neuron # axis
axis(side=2, pos=xmin-0.02, cex.axis=0.6, at=c(0, max(idx)-1), labels=c(min(idx), max(idx)), las=3, padj=1)
mtext(text="Neuron #", side=2, at=length(idx)/2, cex=0.8, font=2, padj=5.5)

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
axis(side=2, pos=xmin-0.31, at=seq(legendBottom+legendStep, legendTop, by=10), 
     labels=round(seq(from=minFR, to=maxFR, length.out=length(seq(legendBottom+legendStep, legendTop, by=10))), 1), las=2, padj=-0.05)

mtext("Firing rate(Z sc)", side=3, line=-1, at=xmin-0.4, cex=0.8, font=2)
mtext(">=", side = 3, line=-2.2, at=xmin-0.58)

# Main title
mtext(text="Post-cue (400ms) firing rate throughout learning", side=3, line=2, cex=1.2, font=2, col="gray50")

#X axis
mtext("NAcC FIRING RATE PER EVENT", side=3, line=-0.5, at=xmin+0.5, font=2, cex = 0.9)
#mtext("(sorted by DS response ratio and DS-evoked firing)", line=-1.3, at=xmin+0.5, font=3, cex=0.8)
mtext("Responded", side=1, line=0, at=xmin+0.25, font=2)
mtext("Missed", side=1, line=0, at=xmin+0.76, font=2)
mtext("DS", side=1, line=-1, at=xmin+0.12, font=3)
mtext("NS", side=1, line=-1, at=xmin+0.38, font=3)
mtext("DS", side=1, line=-1, at=xmin+0.63, font=3)
mtext("NS", side=1, line=-1, at=xmin+0.89, font=3)

# ### Divisions in DS Response ratio w/o taking into account data from diff rats with same performance
# PerfValues <- rle(toplot2$DSRR)
# 
# PerfLines <- cumsum(PerfValues$lengths)
# for(i in 1:length(PerfLines)){
#   segments(y0=PerfLines[i], y1=PerfLines[i],
#            x0=xmin, x1=xmin+1.85, col="gray")
# }
# segments(y0=0, y1=0, x0=xmin, x1=xmin+1.85, col="gray")

### DIVISIONS BY RAT AND DS RESPONSE RATIO
ratsess <- c()
for(i in 1:nrow(toplot2)){
  ratsess <- c(ratsess, as.numeric(paste(toplot2$rat[i], toplot2$session[i], sep="")))
}

PerfValuesbyRAT <- rle(ratsess)
PerfLinesbyRAT <- cumsum(PerfValuesbyRAT$lengths)
for(i in 1:length(PerfLinesbyRAT)){
        segments(y0=PerfLinesbyRAT[i], y1=PerfLinesbyRAT[i],
                 x0=xmin, x1=xmin+2.6, col="gray")
}
segments(y0=0, y1=0, x0=xmin, x1=xmin+2.6, col="gray")

#Make a rectangle where the behavior is going to be plotted
rectMin <- xmin+1.1
rectMax <- xmin+1.8
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
           y0=midTicks[i]+offset, y1=midTicks[i]+offset, col="darkgoldenrod4", lwd=3)
}

#Draw the bars of NSRR performance
b <- c()
for(i in 1:length(NSPerSessPERCENT)){
        b <- NSPerSessPERCENT[i]*lengthPerPoint
        segments(x0=rectMin, x1=rectMin+b, 
                 y0=midTicks[i]-offset, y1=midTicks[i]-offset, col="darkgoldenrod2", lwd=3)
}

#Legend
legend(x=rectMin+0.3, y=20, legend=c("DS", "NS"), fill=c("darkgoldenrod4", "darkgoldenrod2"))

# X Axis of beh plot
axis(side=1, line=-1.2, at=c(rectMin, rectMin+(rectMax-rectMin)/2, rectMax), labels=c(minBeh, (maxBeh-minBeh)/2, maxBeh), padj=-1, cex=0.8)
mtext("PERFORMANCE", side=3, line=-0.5, at=rectMin+(rectMax-rectMin)/2, font=2, cex=0.9)
#mtext("(Response ratio)", line=-1.3, at=rectMin+(rectMax-rectMin)/2, font=3, cex=0.8)
mtext("Response ratio", side=1, line=0, at=rectMin+(rectMax-rectMin)/2, font=2)

### LAST RECTANGLE WITH DISTRIBUTION OF SESSION BY ANIMAL
#Make a rectangle where the animals are going to be plotted
rectMin2 <- xmin+2
rectMax2 <- xmin+2.6
rect(xleft=rectMin2, xright=rectMax2, ybottom=0, ytop=max(idx))

#MAke columns for each animal
ridx <- sort(unique(toplot2$rat))
minrat <- min(ridx)
maxrat <- max(ridx)
ratCols <- seq(minrat, maxrat, by=1)

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
coordCols2 <- c(rectMin2, coordCols)
for(i in 1:length(unique(toplot2$rat))){
  byrat <- filter(toplot2, rat==i)
  for(j in 1:length(byrat$session)){
    ybottom=min(idx[(toplot2$rat==i)&(toplot2$session==j)])-1
    ytop=max(idx[(toplot2$rat==i)&(toplot2$session==j)])
    rect(xleft=coordCols2[i], xright=coordCols2[i]+sizeCol, 
         ybottom=ybottom,
         ytop=ytop,
         col="darkred")
    text(x=coordCols2[i]+sizeCol/2, y=ybottom+(ytop-ybottom)/2, labels=j, col="pink", srt=90, cex=0.7, font=2) 
  }
}

axis(side=1, line=-1.8, at=(coordCols2[-length(coordCols2)])+sizeCol/2, labels=LETTERS[1:length(ridx)], padj=0.5, tick=F, las=2, font=1)
mtext("Animal ID", side=1, line=0, at=rectMin2+(rectMax2-rectMin2)/2, font=2)








#### MEGAPLOT WITH OVERALL EVOLUTION OF PERFORMANCE (TASK ACCURACY) TOO
library(dplyr)
library("colorRamps")

toplot2 <- arrange(toplot, as.numeric(Accuracy), as.numeric(ZDSresponded))
toplot2[,9:12] <- round(toplot2[,9:12], 1)
ZscData <- toplot2[,9:12]

boxplot(ZscData) #Use this boxplot to set the max limit of FR so that outliers don't skew the distribution and therefore obscure the relationship
minFR <- -1.3
maxFR <- 5.5 #Anything above this will be chategorized as this value

colValues <- seq(minFR, maxFR, by=0.1)
colValues <- round(colValues, 1)
mypalette <- colorRampPalette(colors=c("darkblue", "orangered"))(length(colValues))
NAcolor <- "white" #Color for NA cells

mypalette <- c(NAcolor, mypalette)
idx <- 1:nrow(toplot2)

plot.new()
plot.window(xlim=c(1, 4), ylim=c(1, length(idx)))

xmin <- 1.20 #Choose the x position of the beginning of the heatmap (reference for everything else)

#This loop makes the actual grid with colors
for(i in idx){
  for(j in 1:ncol(ZscData)){
    if(j==1){a <- xmin}
    if(j==2){a <- xmin+0.25}
    if(j==3){a <- xmin+0.51}
    if(j==4){a <- xmin+0.76}
    
    #Find Zsc values in the idx of all possible values. 
    #For some reason treating them as numbers gives me wacky results
    FRdatapoint <- ZscData[i,j] 
    if(is.na(FRdatapoint)) {k <- 1} else { #First color in mycolorpalette is white, so if value is NA, make it white by assigning k=1
      if(FRdatapoint > maxFR){k <- length(colValues)} #If FR value is bigger than maxFR (outlier), assign the last value in the palette
      if(FRdatapoint < maxFR){k <- grep(paste("^", FRdatapoint, "$", sep=""), colValues)+1} #Plus one because the lowest value in "mypalette" is white (for NAs)    
    }
    
    if(!is.na(k)==TRUE){colpick <- mypalette[k]} else {colpick <- "white"}
    
    rect(xleft=a, 
         xright=a+0.25, 
         ytop=i, 
         ybottom=i-1, 
         col=colpick, 
         border=NA)    
  }
}

#Neuron # axis
axis(side=2, pos=xmin-0.02, cex.axis=0.6, at=c(0, max(idx)-1), labels=c(min(idx), max(idx)), las=3, padj=1)
mtext(text="Neuron #", side=2, at=length(idx)/2, cex=0.8, font=2, padj=5.5)

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
axis(side=2, pos=xmin-0.31, at=seq(legendBottom+legendStep, legendTop, by=10), 
     labels=round(seq(from=minFR, to=maxFR, length.out=length(seq(legendBottom+legendStep, legendTop, by=10))), 1), las=2, padj=-0.05)

mtext("Firing rate(Z sc)", side=3, line=-1, at=xmin-0.4, cex=0.8, font=2)
mtext(">=", side = 3, line=-2.2, at=xmin-0.58)

# Main title
mtext(text="Post-cue (400ms) firing rate throughout learning", side=3, line=2, cex=1.2, font=2, col="gray50", adj=0.1)

#X axis
mtext("Firing rate per event", side=3, line=-0.5, at=xmin+0.5, font=2)
mtext("DS", side=1, line=0, at=xmin+0.25, font=2)
mtext("NS", side=1, line=0, at=xmin+0.76, font=2)
mtext("Resp.", side=1, line=-1, at=xmin+0.12, font=3)
mtext("Miss.", side=1, line=-1, at=xmin+0.38, font=3)
mtext("Resp.", side=1, line=-1, at=xmin+0.63, font=3)
mtext("Miss.", side=1, line=-1, at=xmin+0.89, font=3)

### Divisions in TASK ACCURACY
PerfValues <- rle(toplot2$Accuracy)

PerfLines <- cumsum(PerfValues$lengths)
for(i in 1:length(PerfLines)){
  segments(y0=PerfLines[i], y1=PerfLines[i],
           x0=xmin, x1=xmin+1.85, col="gray")
}
segments(y0=0, y1=0, x0=xmin, x1=xmin+1.85, col="gray") #Make a line for the 0 y coord. too

#Make a rectangle where the behavior is going to be plotted
rectMin <- xmin+1.1
rectMax <- xmin+1.8
rect(xleft=rectMin, xright=rectMax, ybottom=0, ytop=max(idx))

#Here I'm plotting changes in DSRR. Min possible is 0 and max is 1. Make 0.1 bins.
minBeh <- -0.2
maxBeh <- 0.7
binsBeh <- seq(minBeh, maxBeh, by=0.1)

sizeBin <- (rectMax-rectMin)/(length(binsBeh)-1) #Size of each bin on the plot (minus one because I don't want to count in bin 0)
coordBins <- c()
for(i in 1:length(binsBeh)){
  if(i==1){a <- rectMin} else {a <- rectMin + (i-1)*sizeBin}
  coordBins <- c(coordBins, a)
}

#Make vertical lines with these bins
for(i in 0:(length(coordBins)-1)){
  if(i==0){x <- rectMin} else {x <- coordBins[i]}
  segments(x0=x, x1=x, y0=min(idx)-1.5, y1=max(idx), col="gray90")
}

zeroPoint <- coordBins[binsBeh==0]
segments(x0=zeroPoint, x1=zeroPoint, y0=min(idx)-1.5, y1=max(idx), col="gray")


rect(xleft=rectMin, xright=rectMax, ybottom=0, ytop=max(idx)) #redo rectangle bc it's obscured

#Tickmarks
axis(side=4, pos=rectMax, at=c(0, cumsum(PerfValues$lengths)), labels=FALSE)

#Labels (so that they appear between tickmarks)
midTicks <- c()
for(i in 1:length(PerfLines)){
  if(i == 1){a <- (PerfLines[i]-0)/2} 
  else{a <- PerfLines[i-1]+(PerfLines[i]-PerfLines[i-1])/2}
  midTicks <- c(midTicks, a)
}
axis(side=4, pos=rectMax, at=midTicks, tick=F, labels=length(PerfLines):1, cex.axis=0.6, padj=-3)
mtext("Ranking of session based on performance (Accuracy (DS-NS))", side=4, line=-11, cex.axis=0.7, font=2)

#Determine the length of the ACCURACY bars
PerfPerSessPERCENT <- PerfValues$values*100
lengthPerPoint <- (rectMax-rectMin)/100 #Determine length per ACCURACY point (0.01)

#Draw the bars of performance
b <- c()
for(i in 1:length(PerfPerSessPERCENT)){
  b <- PerfPerSessPERCENT[i]*lengthPerPoint
  segments(x0=zeroPoint, x1=zeroPoint+b, 
           y0=midTicks[i], y1=midTicks[i], col="darkred", lwd=3)
}

# X Axis of beh plot
axis(side=1, line=-1, at=c(rectMin, zeroPoint, rectMax), labels=c(minBeh, 0, maxBeh), padj=-1)
mtext("Performance", side=3, line=-0.5, at=rectMin+(rectMax-rectMin)/2, font=2)
mtext("(Accuracy (DS-NS))", line=-1.3, at=rectMin+(rectMax-rectMin)/2, font=3)






#### MEGAPLOT WITH OVERALL EVOLUTION OF PERFORMANCE TOO (TOTAL RESPONDING)
library(dplyr)
library("colorRamps")

toplot2 <- arrange(toplot, as.numeric(TOTALresp), as.numeric(ZDSresponded))
toplot2[,9:12] <- round(toplot2[,9:12], 1)
ZscData <- toplot2[,9:12]

boxplot(ZscData) #Use this boxplot to set the max limit of FR so that outliers don't skew the distribution and therefore obscure the relationship
minFR <- -1.3
maxFR <- 5.5 #Anything above this will be chategorized as this value

colValues <- seq(minFR, maxFR, by=0.1)
colValues <- round(colValues, 1)
mypalette <- colorRampPalette(colors=c("darkblue", "orangered"))(length(colValues))
NAcolor <- "white" #Color for NA cells

mypalette <- c(NAcolor, mypalette)
idx <- 1:nrow(toplot2)

plot.new()
plot.window(xlim=c(1, 4), ylim=c(1, length(idx)))

xmin <- 1.20 #Choose the x position of the beginning of the heatmap (reference for everything else)

#This loop makes the actual grid with colors
for(i in idx){
  for(j in 1:ncol(ZscData)){
    if(j==1){a <- xmin}
    if(j==2){a <- xmin+0.25}
    if(j==3){a <- xmin+0.51}
    if(j==4){a <- xmin+0.76}
    
    #Find Zsc values in the idx of all possible values. 
    #For some reason treating them as numbers gives me wacky results
    FRdatapoint <- ZscData[i,j] 
    if(is.na(FRdatapoint)) {k <- 1} else { #First color in mycolorpalette is white, so if value is NA, make it white by assigning k=1
      if(FRdatapoint > maxFR){k <- length(colValues)} #If FR value is bigger than maxFR (outlier), assign the last value in the palette
      if(FRdatapoint < maxFR){k <- grep(paste("^", FRdatapoint, "$", sep=""), colValues)+1} #Plus one because the lowest value in "mypalette" is white (for NAs)    
    }
    
    if(!is.na(k)==TRUE){colpick <- mypalette[k]} else {colpick <- "white"}
    
    rect(xleft=a, 
         xright=a+0.25, 
         ytop=i, 
         ybottom=i-1, 
         col=colpick, 
         border=NA)    
  }
}

#Neuron # axis
axis(side=2, pos=xmin-0.02, cex.axis=0.6, at=c(0, max(idx)-1), labels=c(min(idx), max(idx)), las=3, padj=1)
mtext(text="Neuron #", side=2, at=length(idx)/2, cex=0.8, font=2, padj=5.5)

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
axis(side=2, pos=xmin-0.31, at=seq(legendBottom+legendStep, legendTop, by=10), 
     labels=round(seq(from=minFR, to=maxFR, length.out=length(seq(legendBottom+legendStep, legendTop, by=10))), 1), las=2, padj=-0.05)

mtext("Firing rate(Z sc)", side=3, line=-1, at=xmin-0.4, cex=0.8, font=2)
mtext(">=", side = 3, line=-2.2, at=xmin-0.58)

# Main title
mtext(text="Post-cue (400ms) firing rate throughout learning", side=3, line=2, cex=1.2, font=2, col="gray50", adj=0.1)

#X axis
mtext("Firing rate per event", side=3, line=-0.5, at=xmin+0.5, font=2)
mtext("DS", side=1, line=0, at=xmin+0.25, font=2)
mtext("NS", side=1, line=0, at=xmin+0.76, font=2)
mtext("Resp.", side=1, line=-1, at=xmin+0.12, font=3)
mtext("Miss.", side=1, line=-1, at=xmin+0.38, font=3)
mtext("Resp.", side=1, line=-1, at=xmin+0.63, font=3)
mtext("Miss.", side=1, line=-1, at=xmin+0.89, font=3)

### Divisions in DS Response ratio
PerfValues <- rle(toplot2$DSRR)

PerfLines <- cumsum(PerfValues$lengths)
for(i in 1:length(PerfLines)){
  segments(y0=PerfLines[i], y1=PerfLines[i],
           x0=xmin, x1=xmin+1.85, col="gray")
}
segments(y0=0, y1=0, x0=xmin, x1=xmin+1.85, col="gray")

#Make a rectangle where the behavior is going to be plotted
rectMin <- xmin+1.1
rectMax <- xmin+1.8
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

#Tickmarks
axis(side=4, pos=rectMax, at=c(0, cumsum(PerfValues$lengths)), labels=FALSE)

#Labels (so that they appear between tickmarks)
midTicks <- c()
for(i in 1:length(PerfLines)){
  if(i == 1){a <- (PerfLines[i]-0)/2} 
  else{a <- PerfLines[i-1]+(PerfLines[i]-PerfLines[i-1])/2}
  midTicks <- c(midTicks, a)
}
axis(side=4, pos=rectMax, at=midTicks, tick=F, labels=length(PerfLines):1, cex.axis=0.6, padj=-3)
mtext("Ranking of session based on performance (Total response ratio)", side=4, line=-11, cex.axis=0.7, font=2)

#Determine the length of the DSRR bars
PerfPerSessPERCENT <- round(PerfValues$values, 2)*100
lengthPerPoint <- (rectMax-rectMin)/100 #Determine length per DSRR point (0.01)

#Draw the bars of performance
b <- c()
for(i in 1:length(PerfPerSessPERCENT)){
  b <- c(b, PerfPerSessPERCENT[i]*lengthPerPoint)
  segments(x0=rectMin, x1=rectMin + b, 
           y0=midTicks, y1=midTicks, col="darkred", lwd=3)
}

# X Axis of beh plot
axis(side=1, line=-1, at=c(rectMin, rectMin+(rectMax-rectMin)/2, rectMax), labels=c(minBeh, (maxBeh-minBeh)/2, maxBeh), padj=-1)
mtext("Performance", side=3, line=-0.5, at=rectMin+(rectMax-rectMin)/2, font=2)
mtext("(DS response ratio)", line=-1.3, at=rectMin+(rectMax-rectMin)/2, font=3)
