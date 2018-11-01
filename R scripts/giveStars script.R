giveStars <- function(p.vals){
        sapply(seq(1, length(p.vals)), function(x){
                if(p.vals[x]>0.05){a <- " "}
                if(p.vals[x]<=0.05){a <- "*"}
                if(p.vals[x]<=0.01){a <- "**"}
                if(p.vals[x]<=0.001){a <- "***"}
                return(a)
        })
}

save(giveStars, file="E:/Dropbox/NMDA/R functions/giveStars.R")
save(giveStars, file=paste("E:/Dropbox/NMDA/EXP3_NAc FR acquisition/R functions/giveStars.R"))
save(giveStars, file=paste("E:/Dropbox/NMDA/EXP4_Unilateral AP5/R functions/giveStars.R"))
