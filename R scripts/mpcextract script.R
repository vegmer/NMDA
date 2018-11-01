### Function to extract data from MedPC files

mpcextract <- function(filename){
        rawin = scan(filename, what = "raw", sep = "\n")
        singlevars = grep("^[A-Z]{1}:[[:space:]]{0,}[0-9]", rawin)
        arrayvars = grep("^[A-Z]{1}:$", rawin)
        listout = list(info = as.vector(rawin[c(1:10)]))
        if( length(singlevars) > 0) {
                svarlist = as.numeric(substr(rawin[singlevars], start= 3, stop = 100))
                svarlist = as.list(svarlist)
                names(svarlist) = tolower(substr(rawin[singlevars], start= 1, stop = 1))
        }
        
        if( length(singlevars) == 0){svarlist = list(singlevars = "No single variables found")}
        if(length(arrayvars) > 0){
                for(j in 1:length(arrayvars)){
                        if(j < length(arrayvars)) {iAB = (arrayvars[j]+1):(arrayvars[j+1] - 1)}
                        if(j == length(arrayvars)) {iAB = (arrayvars[j]+1):(length(rawin))}
                        arrayblock = rawin[iAB]
                        avars = substr(arrayblock, start = 8, stop = 100)
                        avars = unlist(strsplit(avars, split = "[[:space:]]{1,}"))
                        avars = as.numeric(avars)
                        avars = avars[!is.na(avars)]
                        if(j == 1 ){avarlist = list(avars)}
                        else {avarlist[[j]] = avars}
                } #End for(j . . .
                names(avarlist) = tolower(substr(rawin[arrayvars], start= 1, stop = 1))
        }# End if(length(arrayvars . . .
        if(length(arrayvars) == 0){avarlist = list(arrayvars = "No array variables found")}
        listout = c(listout, svarlist, avarlist)
        return(listout)
}


save(mpcextract, file="E:/Dropbox/NMDA/R Functions/mpcextract.R")
