#This was adapted from https://hpc-forge.cineca.it/files/CoursesDev/public/2015/Data_Analysis_School/script%20R/modelBased_genData.R
#library(e1071)
#library(parallel)
library(foreach)
library(doMC)
library(mclust)

parallelMcl <- function (
                         data,
                         ncores=4,
                         G=2:6,
                         modelNames=c("EII","EEI","EEE","VVV")
                         ) {


    param=expand.grid(G=G,modelNames=modelNames) 

    parallel.function <- function(i) { Mclust( data=data, G=as.character(param$G[i]), modelNames=as.character(param$modelNames[i]),prior=priorControl()) }

    registerDoMC(ncores)

    b=Sys.time()
    all <- foreach(i = 1:nrow(param)) %dopar% parallel.function(i)
    e=Sys.time()
    best = all[[which.max(sapply(all,function(x){x$bic}))]]
    time=e-b
    rm(e,b)
    print(time)

    print(ls(best))

    entropy=sum((table(best$cl)/length(best$cl))*log(table(best$cl)/length(best$cl)))
    print(entropy)
    return(list(all=all, best=best, entropy=entropy, time=time))
}

