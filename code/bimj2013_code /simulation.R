######################################################################
# This script:
# - generates R datasets in each of 16 simulation scenarios
# - save average performance metrics of PSM and PSW over R
#   for each propensity score (ps) model



######### load functions
source("functions.R")

######### SIMULATION STUDY  ############

######### set simulation  parameters
# set scenarios:
scenarios <- c("Aa", "Ad", "Ae", "Ag", "Ca", "Cd", "Ce", "Cg", "Fa", "Fd", "Fe", "Fg", "Ga", "Gd", "Ge", "Gg")


# set target parameter
par <- "ATE"

# set performance metrics
# perfmetrics<-names(funsim( funcov(500,"A","a"),"logit",par))
perfmetrics <- names(funsim_opt(funcov(1000, "A", "a"), "logit", par, NULL))

# set number of replications (set to 100 to speed up)
R <- 500


####### Simulation cycle
for (size in c(500, 1000, 2000)) {
    # to speed up run the simulation only for a size (e.g. 500)


    ## ps model results

    #  logit
    lgresults <- matrix(0, length(perfmetrics), 1)
    timestart <- Sys.time()
    for (i in 1:length(scenarios)) {
        pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
        partialresults <- matrix(0, length(perfmetrics), 1)
        for (j in 1:R)
        {
            partialresults <- partialresults + try(as.matrix(simplify2array(funsim(
                funcov(size, substr(scenarios[i], 1, 1), substr(scenarios[i], 2, 2)), "logit", par
            ))))
            setTxtProgressBar(pb, i)
        }
        lgresults <- cbind(lgresults, partialresults / R)
    }
    lgresults <- lgresults[, -1]
    colnames(lgresults) <- paste(scenarios)


    timeend <- Sys.time()
    exectime <- timeend - timestart
    exectime

    # write.table(lgresults, file = file.path("results",paste("logit","R",R,"size",size,".txt",sep="")))


    # bag
    bagresults <- matrix(0, length(perfmetrics), 1)
    for (i in 1:length(scenarios)) {
        pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
        partialresults <- matrix(0, length(perfmetrics), 1)
        for (j in 1:R)
        {
            partialresults <- partialresults + try(as.matrix(simplify2array(funsim(
                funcov(size, substr(scenarios[i], 1, 1), substr(scenarios[i], 2, 2)), "bag", par
            ))))
            setTxtProgressBar(pb, i)
        }
        bagresults <- cbind(bagresults, partialresults / R)
    }
    bagresults <- bagresults[, -1]
    colnames(bagresults) <- paste(scenarios)

    # write.table(bagresults, file = file.path("results",paste("bag","R",R,"size",size,".txt",sep="")))

    # nb
    nbresults <- matrix(0, length(perfmetrics), 1)
    for (i in 1:length(scenarios)) {
        pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
        partialresults <- matrix(0, length(perfmetrics), 1)
        for (j in 1:R)
        {
            partialresults <- partialresults + try(as.matrix(simplify2array(funsim(
                funcov(size, substr(scenarios[i], 1, 1), substr(scenarios[i], 2, 2)), "nb", par
            ))))
            setTxtProgressBar(pb, i)
        }
        nbresults <- cbind(nbresults, partialresults / R)
    }
    nbresults <- nbresults[, -1]
    colnames(nbresults) <- paste(scenarios)


    # write.table(nbresults, file = file.path("results",paste("nb","R",R,"size",size,".txt",sep="")))


    ##############################################################
    ################  cycle using optimal tuning for tree, rf and nn
    # note:  use funsim_opt instead of funsim

    # scenarios<-c("Aa" ,"Ad", "Ae",
    # "Ag", "Ca", "Cd", "Ce", "Cg","Fa", "Fd", "Fe", "Fg","Ga", "Gd", "Ge", "Gg")

    # set optimal tuning parameters
    # c(mtry, cp, size, decay, size, decay)
    tunopt <- list(
        c(2, 0.002, 5, 2, 1, 0.007), c(5, 0.02, 1, 0.0001, 1, 0.0001), c(1, 0.0012, 5, 0.38, 5, 0.38),
        c(5, 0.0015, 11, 3.17, 1, 0.001),
        c(2, 0.037, 3, 0.52, 3, 0.52), c(2, 0.001, 7, 2.77, 13, 0.43), c(2, 0.034, 1, 3, 15, 0.83),
        c(2, 0.034, 3, 0.52, 3, 0.52), c(3, 0.002, 1, 0.02, 9, 0.35), c(3, 0.002, 1, 0.02, 9, 0.35), c(3, 0.0023, 11, 1.23, 9, 0.35),
        c(2, 0.015, 3, 0.05, 3, 0.05), c(2, 0.03, 3, 0.26, 3, 0.26), c(2, 0.05, 1, 0.3, 1, 0.3), c(2, 0.05, 1, 0.3, 1, 0.3),
        c(7, 0.002, 1, 0.005, 1, 0.005)
    )


    # tree
    treeresults <- matrix(0, length(perfmetrics), 1)
    timestart <- Sys.time()
    for (i in 1:length(scenarios)) {
        pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
        # partialresults<-matrix(0,length(perfmetrics),1)
        partialresults <- list()
        for (j in 1:R) {
            set.seed(j)
            prova <- funcov(size, substr(scenarios[i], 1, 1), substr(scenarios[i], 2, 2))

            partialresults[[j]] <- tryCatch(funsim_opt(prova, "tree", par, tunopt[[i]]),
                error = function(e) paste("tree error at ", j, sep = "")
            )
            setTxtProgressBar(pb, i)
        }

        partialresults <- partialresults[lapply(partialresults, function(x) class(x)) == "list"]
        partialresults <- sapply(partialresults, simplify2array)
        # partialresults <- sapply(partialresults,simplify2array)
        treeresults <- cbind(treeresults, apply(as.matrix(partialresults), 1, mean))
    }
    treeresults <- treeresults[, -1]
    colnames(treeresults) <- paste(scenarios)

    # rf

    rfresults <- matrix(0, length(perfmetrics), 1)
    timestart <- Sys.time()
    for (i in 1:length(scenarios)) {
        pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
        # partialresults<-matrix(0,length(perfmetrics),1)
        partialresults <- list()
        for (j in 1:R) {
            set.seed(j)
            prova <- funcov(size, substr(scenarios[i], 1, 1), substr(scenarios[i], 2, 2))

            partialresults[[j]] <- tryCatch(funsim_opt(prova, "randomforest", par, tunopt[[i]]),
                error = function(e) paste("rf error at ", j, sep = "")
            )
            setTxtProgressBar(pb, i)
        }

        partialresults <- partialresults[lapply(partialresults, function(x) class(x)) == "list"]
        partialresults <- sapply(partialresults, simplify2array)
        # partialresults <- sapply(partialresults,simplify2array)
        rfresults <- cbind(rfresults, apply(as.matrix(partialresults), 1, mean))
    }
    rfresults <- rfresults[, -1]
    colnames(rfresults) <- paste(scenarios)
    # rfresults
    # timeend<-Sys.time();exectime<i-timeend-timestart
    # exectime


    nnresults <- matrix(0, length(perfmetrics), 1)
    for (i in 1:length(scenarios)) {
        pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
        partialresults <- list()
        for (j in 1:R) {
            set.seed(j)
            prova <- funcov(size, substr(scenarios[i], 1, 1), substr(scenarios[i], 2, 2))

            partialresults[[j]] <- tryCatch(funsim_opt(prova, "nn", par, tunopt[[i]]),
                error = function(e) paste("net error at ", j, sep = "")
            )
            setTxtProgressBar(pb, i)
        }
        partialresults <- partialresults[lapply(partialresults, function(x) class(x)) == "list"]
        partialresults <- sapply(partialresults, simplify2array)
        # partialresults <- sapply(partialresults,simplify2array)
        nnresults <- cbind(nnresults, apply(as.matrix(partialresults), 1, mean))
    }
    nnresults <- nnresults[, -1]
    colnames(nnresults) <- paste(scenarios)


    #### save workspace before twang  (twang is about 10 times slower than the others)

    save.image(file = file.path("results", paste("R", R, " size", size, " par", par, ".Rdata", sep = "")))

    ### gbm twang

    twresults <- matrix(0, length(perfmetrics), 1)
    for (i in 1:length(scenarios)) {
        pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
        partialresults <- matrix(0, length(perfmetrics), 1)
        for (j in 1:R)
        {
            partialresults <- partialresults + try(as.matrix(simplify2array(
                funsim(
                    funcov(size, substr(scenarios[i], 1, 1), substr(scenarios[i], 2, 2)), "gbmtwang", par
                )
            )))
            setTxtProgressBar(pb, i)
        }
        twresults <- cbind(twresults, partialresults / R)
    }
    twresults <- twresults[, -1]
    colnames(twresults) <- paste(scenarios)


    # write.table(twresults, file = file.path("results",paste("twang","R",R,"size",size,".txt",sep="")))

    #### save workspace after gbm twang

    save.image(file = file.path("results", paste("R", R, " size", size, " par", par, ".Rdata", sep = "")))
} # end loop on size


######### The end ####################