is.even <- function(x){
        sapply(seq(1, length(x)), function(i){
                if(round(x[i]/2, 0)==x[i]/2)
                {a <- TRUE} else {a <- FALSE}
                return(a)
        })
}

save(is.even, file="E:/Dropbox/NMDA/R Functions/is.even.R")

