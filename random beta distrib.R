
## ss: site-species matrix
## envdata: environmental data
random.beta <- function(ss, envdata, sampid, nit, runmod = T, partout = NULL) {

    require(future.apply)
    require(mgcv)
    require(ranger)
    require(diptest)
    require(pdp)
    require(fields)

    varall <- names(envdata)[-1] # assume all columns of envdata will be used
                                # for modeling except first one
    varall_sd <- varall # define predictor variables that for SD

    ## drop xslope and pctshrb2019ws from SD computation
    ## because mean and SD for these variables are highly correlated
    dlist <- c("xslope")
    for (i in dlist)
        varall_sd <- varall_sd[-which(i == varall_sd)]

    tnames <- names(ss)[-1] # assume all columns of ss are taxon names
                                        #except for first one

    ## convert ss to presence absence
    mat0 <- as.matrix(ss[, tnames])
    mat1 <- mat0 >0
    dim(mat1) <- dim(mat0)
    ss2 <- data.frame(ss[,1], mat1)
    names(ss2) <- c(sampid, tnames)

    ## merge ss with envdata
    ss <- merge(ss2, envdata, by = sampid)
    cat("Number of samples:", nrow(ss), "\n")

    ## output number of ASVs
    cat("Total number of asv:", sum(apply(ss[, tnames], 2, any)), "\n")

    ## get separation distance between samples
    dist0 <- rdist.earth(ss[, c("lon.dd83", "lat.dd83")], miles = F)
    diag(dist0) <- NA
    ## these are commands to save dist0 to disk to reduce
    ## computation when running the code repeatedly
#   save(dist0, file = "dist0.rda")
#    load("dist0.rda")

    ## functions to compute summaries of environment variables in regions
    getenv <- function(isamp, ss) apply(ss[isamp,], 2, mean)
    getenv_sd <- function(isamp, ss) apply(ss[isamp,], 2, sd)

    ## set random seed here
    set.seed(21)
    ## select subsample (1/10 of available data)
    nsub <- round(nrow(ss)/12)
    cat("N sub:", nsub, "\n")

    ## identify nsite nearest neighbors  for each site
    nsite <- 11 # number of nearest neighbors

    isamp <- as.list(rep(NA, times = nit))
    maxdist <- as.list(rep(NA, times = nit))

    ## save nearest neighbor selections and spatial grain size
    for (k in 1:nit){
        ipick <- sample(nrow(ss), nsub, replace = F)
        isamp[[k]] <- matrix(NA, nrow = nsite+1, ncol = nsub)
        maxdist[[k]] <- rep(NA, times = nsub)
        for (i in 1:nsub) {
            x <- order(dist0[ipick[i],])
            isamp[[k]][,i] <- c(x[1:nsite], ipick[i])
            maxdist[[k]][i] <- dist0[ipick[i],x[nsite]]
        }
    }

    ## uncomment to print example values of alpha, beta, gamma, and
    ## spatial grain size
#    getbeta2 <- function(isamp, ss) {
#        1 - mean(apply(ss[isamp,], 1, sum))/
#            sum(apply(ss[isamp, ], 2, any))}
#    getgamma2 <- function(isamp, ss) {
#        sum(apply(ss[isamp, ], 2, any))}
#    getalpha2 <- function(isamp, ss) {
#        mean(apply(ss[isamp,], 1, sum))}
#    print(summary(maxdist[[1]]))
#    print(summary(apply(isamp[[1]], 2, getalpha2, ss = ss[, tnames])))
#    print(summary(apply(isamp[[1]], 2, getgamma2, ss = ss[, tnames])))
#    print(summary(apply(isamp[[1]], 2, getbeta2, ss = ss[, tnames])))
#    stop()

    ## run one iteration of sampling to set up prediction matrix
    np <- 40  # number of points for computing partial dependence

    ## compute mean environmental conditions
    ip <- 1
    env.mn <- t(apply(isamp[[ip]], 2, getenv, ss=ss[, varall]))
    env.mn <- cbind(env.mn, maxdist[[1]])
    dimnames(env.mn)[[2]][length(dimnames(env.mn)[[2]])] <- "maxdist"

    ## compute SD environmental conditions
    env.mn_sd <- t(apply(isamp[[ip]], 2,getenv_sd,ss= ss[, varall_sd]))
    names0 <- dimnames(env.mn_sd)[[2]]
    dimnames(env.mn_sd)[[2]] <- paste(names0, "sd", sep = "_")

    ## combine mean and SD environmental conditions in one matrix
    alldat <- cbind(env.mn, env.mn_sd)
    print(ncol(alldat))

    ## calculate median and range of each predictor variable
    mn0 <- apply(alldat, 2, median)
    r0 <- apply(alldat, 2, range)

    ## set up matrices of predictors in which only one variable varies
    ## across its range and the rest are set to median values
    ## these matrices are used to compute partial dependence relationships
    mat0 <- rep(mn0, times = rep(np, times = length(mn0)))
    dim(mat0) <- c(np, length(mn0))
    dimnames(mat0)[[2]] <- names(mn0)
    predmat <- as.list(rep(NA, times = length(mn0)))
    names(predmat) <- names(mn0)
    for (i in 1:length(mn0)) {
        mat1 <- mat0
        mat1[, names(mn0)[i]] <- seq(r0[1, names(mn0)[i]],
                                     r0[2, names(mn0)[i]], length = np)
        predmat[[i]] <- mat1
    }

    ## subroutine to calculate partial dependence relationships
    getpartials <- function(isamp0, maxdist0, ss, predmat, varall,varall_sd,
                            tnames) {
        getgamma <- function(isamp, ss) sum(apply(ss[isamp, ], 2, any))
        getalpha <- function(isamp, ss)  mean(apply(ss[isamp,], 1, sum))
        getbeta <- function(isamp, ss) {
            (1 - mean(apply(ss[isamp,], 1, sum))/
                sum(apply(ss[isamp, ], 2, any)))}
        getenv <- function(isamp, ss) apply(ss[isamp,], 2, mean)
        getenv_sd <- function(isamp, ss) apply(ss[isamp,], 2, sd)

        ## set up envdata for subsampled data
        env.mn <- t(apply(isamp0, 2, getenv, ss=ss[, varall]))
        env.mn <- cbind(env.mn, maxdist0)
        dimnames(env.mn)[[2]][length(dimnames(env.mn)[[2]])] <- "maxdist"
        env.mn_sd <- t(apply(isamp0, 2, getenv_sd, ss= ss[, varall_sd]))
        names0 <- dimnames(env.mn_sd)[[2]]
        dimnames(env.mn_sd)[[2]] <- paste(names0, "sd", sep = "_")

        ## uncomment below to get partial dependence relationships
        ## for beta, alpha, or gamma
        beta <- apply(isamp0, 2, getbeta,ss= ss[, tnames])
        ## subbing in alpha diversity
#        beta <-  apply(isamp0, 2, getalpha,ss= ss[, tnames])
        ## subbing in gamma
#        beta <-  apply(isamp0, 2, getgamma,ss= ss[, tnames])

        alldat <- cbind(beta, env.mn, env.mn_sd)
        dimnames(alldat)[[2]][1] <- "beta"

        ## fit random forest model
        mod <- ranger(data = alldat, dependent.variable.name = "beta",
                      num.trees = 2000) #mtry = round((ncol(alldat)-1)/2),
                    #  importance = "permutation")

        ## calculate partial dependence relationships
        partout <- matrix(NA, ncol = length(predmat), nrow = nrow(predmat[[1]]))
        for (i in 1:length(predmat)) {
            partout[,i] <- predict(mod, data = predmat[[i]])$predictions
        }
        attr(partout, "r2") <- mod$r.squared
        return(partout)
    }

    ## single function call for testing. uncomment to
    ## test code before running in parallel over multiple processors
#    getpartials(isamp[[1]], maxdist[[1]], ss, predmat, varall, tnames)

    if (runmod) {
        ## set up to run on multiple cores. Here 15 are used
        plan(multisession, workers = 15)
        partout <- future_mapply(getpartials,isamp, maxdist,
                                 MoreArgs = list(ss = ss, predmat = predmat,
                                     varall = varall,
                                     varall_sd = varall_sd,
                                     tnames = tnames),
                                 future.seed = TRUE)
        dim(partout) <- c(np, length(varall_sd) + length(varall) + 1, nit)
        plan(sequential)
        return(partout)
    }
    else {
        ## do post processing...make plots and calculate summary statistics

        ## example of land use relationships with mean and SD of environmental
        luenv <- F
        if (luenv) {
#            png(width = 6, height = 5, pointsize = 9, units = "in", res = 600,
#                file = "luvar.png")
#            dev.new()
            par(mar = c(4,4,1,1), mfrow = c(2,2), bty = "l", mgp = c(2.3,1,0))
            plot(alldat[, "pctcrop2019ws"], alldat[, "anc.result"], axes = F,
                 xlab = "% row crop", ylab = expression(Mean~ANC~(mu*eq/cm)),
                 pch = 21, col = "grey39", bg = "white")
            logtick.exp(0.001, 10, 2, c(F,F))
            axis(1)
            plot(alldat[, "pctcrop2019ws"], alldat[, "ntl.diss"], axes = F,
                 xlab = "% row crop", ylab = expression(Mean~N[diss]~(mu*g/L)),
                 pch = 21, col = "grey39", bg = "white")
            logtick.exp(0.001, 10, 2, c(F,F))
            axis(1)
            plot(alldat[, "pctcrop2019ws"], alldat[, "anc.result_sd"], axes = T,
                 xlab = "% row crop", ylab = expression(SD~ln(ANC)),
                 pch = 21, col = "grey39", bg = "white")
            plot(alldat[, "pctcrop2019ws"], log(alldat[, "ppt_sd"]), axes = T,
                 xlab = "% row crop", ylab = expression(SD~ln(N[diss])),
                 pch = 21, col = "grey39", bg = "white")
            dev.off()
            stop()
        }

        ## show example of subsampled sites in a map
        mapexamp <- F
        if (mapexamp) {
            ## project locations
            require(maps)
            require(mapproj)
            require(plotrix)
#            png(width = 5, height = 3, file = "mapex.png",
#                units = "in", res = 600, pointsize = 7)
            tiff(width = 5, height = 3, file = "mapex.tif",
                 units = "in", res = 600, pointsize = 7,
                 compression = "none")
#            dev.new()
            par(mar = c(0,0,0,0))
            map("state", proj = "albers", par = c(30,40))
            pout0 <- mapproject(ss$lon.dd83, ss$lat.dd83, proj = "")
            points(pout0$x, pout0$y, col = "grey", pch = 21, bg = "white")
            points(pout0$x[isamp[[1]][12,]], pout0$y[isamp[[1]][12,]], col = "grey", pch = 21, bg = "black")

            ip <- 9
            points(pout0$x[isamp[[1]][1:11,ip]], pout0$y[isamp[[1]][1:11,ip]], col = "grey", pch = 21, bg = "blue")
            points(pout0$x[isamp[[1]][12,ip]], pout0$y[isamp[[1]][12,ip]], col = "grey", pch = 21, bg = "red")

            dist1 <- as.matrix(dist(cbind(pout0$x[isamp[[1]][,ip]], pout0$y[isamp[[1]][,ip]])))
            max1 <- max(dist1[12,])
            draw.circle(pout0$x[isamp[[1]][12,ip]], pout0$y[isamp[[1]][12,ip]],
                        radius = max1)

            dev.off()
            stop()
        }

        ## calculate significance and importance of predictors

        ## set up storage locations for predictor stats
        varp <- names(mn0)
        mag <- rep(NA, times = length(varp))
        names(mag) <- varp
        sig <- rep(NA, times = length(varp))
        names(sig) <- varp
        dir <- rep(NA, times = length(varp))
        names(dir) <- varp

        for (j in 1:length(varp)) {
            xmat <- partout[, which(names(mn0) == varp[j]),]
            ## center each realtionship by mean value
            mnpred <- apply(xmat, 2, mean)
            for (i in 1:ncol(xmat)) xmat[,i] <- xmat[,i] - mnpred[i]

            ## compute upper and lower percentiles and median
            bnds <- matrix(NA, ncol = 3, nrow = nrow(xmat))
            for (i in 1:nrow(xmat))
                bnds[i,] <- quantile(xmat[i,], probs = c(0.005, 0.5, 0.995))
            ## test for significant trend that differs from zero slope
            sig[j] <- any(min(bnds[,3]) < bnds[,1]) |
                any(max(bnds[,1]) > bnds[,3])
            mag[j] <- diff(range(bnds[,2]))
            dir[j] <- as.numeric((bnds[nrow(bnds),2] - bnds[1,2])>0)
        }
        ## summary of results in a data frame
        dfout <- data.frame(var = names(mag), mag = as.vector(mag),
                            sig = as.vector(sig), dir = as.vector(dir))

#        imp0 <- apply(impout, 1, mean)
#        dfout <- merge(dfout, data.frame(var = names(imp0),
#                                         imp = as.vector(imp0)),
#                       by = "var")

        ## select significant predictors
        dfout2 <- dfout[dfout$sig,]
        ## sort by importance
        dfout2 <- dfout2[rev(order(dfout2$mag)),]
        print(dfout2)

        ## uncomment here to output selected variables
        ## to varsel
        ##return(dfout2)

        ## number of env heterogeneity variables that show
        ## positive and negative relationships
        incvec <- regexpr("_sd", dfout2$var) != -1
        print(table(dfout2$dir[incvec]))
        print(table(dfout2$dir[!incvec]))

        print("EH variables that decrease:")
        print(dfout2$var[incvec & dfout2$dir == 0])
        print("Mean variables that decrease")
        print(dfout2$var[!incvec & dfout2$dir == 0])
        print("Mean variables that increase")
        print(dfout2$var[!incvec & dfout2$dir == 1])

        w <- regexpr("_sd", dfout2$var)
        selvec <- w != -1
        root1 <- dfout2$var
        root1[selvec] <- substring(root1[selvec], 1, w[selvec]-1)
        root0 <- sort(unique(root1))
        lab0 <- rep(NA, times = length(root0))

        ## importance plot
        doplot1 <- F
        if (doplot1) {
            lab0 <- c("Air temp",  "ANC",
                      "Cond-adjusted", "DOC",
                       "Elevation",
                      "Wet N deposition",  "Manure",
                      "Water temperature", "Nitrate-nitrite",
                       "% row crop",
                       "Precipitation",
                      "Bedrock depth",
                      "Runoff", "SiO2",
                      "Wet days",
                      "Catchment area", "Slope")
            names(lab0) <- root0
            print(lab0)
            ## grp labels are the type of predictor
            grp <- c("C", "Q", "Q",
                     "Q", "T","A",
                     "A", "C", "Q",
                     "A",
                     "C", "G", "C",
                     "G", "C",
                      "T", "T")
            names(grp) <- root0
            print(grp)

            print(table(grp))

            labsav <- lab0[root1]
            labsav[selvec] <- paste(labsav[selvec], "(SD)")

            grpsav <- grp[root1]
            print(grpsav)

            library(RColorBrewer)
            palette <- brewer.pal(5, "Accent")
            print(palette)
            names(palette) <- sort(unique(grpsav))
            print(palette)
            colsav <- palette[grpsav]
            print(colsav)

#            png(width = 3, height = 4, pointsize = 6, units = "in",
#                "impplot.png", res = 600)
            tiff(width = 3, height = 4, pointsize = 6, units = "in",
                "impplot.tif", res = 600)

            par(mar = c(4,11,1,1), las = 1, mgp = c(2.3,1,0))
            barplot(dfout2$mag, horiz = T, names.arg = labsav,
                    col = colsav, xlab = "Importance")

            dev.off()

        }

        varnames <- dfout2$var

        cat("Number of sig:", nrow(dfout2), "\n")
        cat("Number of sd var:", sum(selvec), "\n")
        dfout3 <- dfout2[!selvec,]
        print(dfout3)

        ## partial dependence plots with uncertainty interval
        grey.t <- adjustcolor("grey", alpha.f = 0.5)

        meanplot <- T
        chemplot <- F
        if (meanplot) {
#            png(width = 6.5, height = 2, pointsize = 10, units = "in",
#                res = 600, file = "srplot.png")
            tiff(width = 6.5, height = 2, pointsize = 10, units = "in",
                res = 600, file = "srplotalpha.tif")
            par(mar = c(4,4,1,1), mfrow = c(1,3), mgp = c(2.4,1,0), bty = "l")

            varplot <- c("pctcrop2019ws",
                         "inorgnwetdep.2008ws",
                         "manurews")
            lab0 <- c(expression(symbol('%')~row~crop),
                      expression(Wet~N~deposition~(kg~N/ha/yr)),
                      expression(Manure~(kg~N/ha/yr)))
            logt <- c(F, F, F)
        }
        else {
            if (chemplot) {
#                png(width = 6.5, height = 2, pointsize = 10, units = "in",
#                    res = 600, file = "srplotchem.png")
                tiff(width = 6.5, height = 2, pointsize = 10, units = "in",
                    res = 600, file = "srplotchem.tif")
                par(mar = c(4,4,1,1), mfrow = c(1,3), mgp = c(2.4,1,0),
                    bty = "l")
                varplot <- c("anc.result", "nox", "cond2")
                lab0 <- c(expression(ANC~(mu*eq/L)),
                          expression(NO[x]~(mu*g/L)),
                          expression(Cond-adjusted~(mu*S/cm)))
                logt <- c(T,T, T)
            }
            else {
#                png(width = 5, height = 2, pointsize = 9, units = "in",
#                    res = 600, file = "srplotsd.png")
                tiff(width = 5, height = 2, pointsize = 9, units = "in",
                    res = 600, file = "srplotsd.tif")
                par(mar = c(4,4,1.5,1), mfrow = c(1,2), mgp = c(2.4,1,0), bty = "l")
                varplot <- c("anc.result_sd",
                             "ppt_sd")
                lab0 <- c(expression(SD~ln(ANC)),
                          expression(SD~Precipitation))
                logt <- c(F,F)
#                varplot <- dfout2$var[19:26]
#                lab0 <- varplot
#                logt <- rep(F, times = length(varplot))
            }
        }

        ## get plot limits
        ylim <- numeric(0)
        for (j in varplot) {
            isel <- which(j == varplot)
            xmat <- partout[, which(names(mn0) == j),]
            mnpred <- apply(xmat, 2, mean)
            for (i in 1:ncol(xmat)) xmat[,i] <- xmat[,i] - mnpred[i]
            bnds <- matrix(NA, ncol = 3, nrow = nrow(xmat))
            for (i in 1:nrow(xmat))
                bnds[i,] <- quantile(xmat[i,], probs = c(0.005, 0.5, 0.995))
            ylim <- range(c(ylim, bnds))
        }

        for (j in varplot) {
            isel <- which(j == varplot)
            xmat <- partout[, which(names(mn0) == j),]
            mnpred <- apply(xmat, 2, mean)
            for (i in 1:ncol(xmat)) xmat[,i] <- xmat[,i] - mnpred[i]
            bnds <- matrix(NA, ncol = 3, nrow = nrow(xmat))
            for (i in 1:nrow(xmat))
                bnds[i,] <- quantile(xmat[i,], probs = c(0.005, 0.5, 0.995))
            plot(predmat[[j]][, j], bnds[,2],
                 type = "n",
                 ylim = ylim, xlab = lab0[isel],
                 ylab = expression(beta~diversity),
                 axes = F)
            if (logt[isel]) logtick.exp(0.001, 10, c(1), c(F,F))
            else   axis(1)
            axis(2)
            box(bty ="l")
            polygon(c(predmat[[j]][, j],
                      rev(predmat[[j]][, j])),
                    c(bnds[,1], rev(bnds[,3])), col = grey.t, border = NA)
#                lines(predmat[[dfout3$var[j]]][, dfout3$var[j]],bnds[,1], lty = "dashed")
                                        #                lines(predmat[[dfout3$var[j]]][, dfout3$var[j]],bnds[,3], lty = "dashed")
            lines(predmat[[j]][, j],bnds[,2])
        }
        dev.off()
    }
    return()
}


## this function sets up ss and envdata data frames from
## raw data frames
selvars <- function(ss, envdata, runmod) {

    varscat <- c(names(envdata)[53:136], "lat.dd83", "lon.dd83",
                 "elevation", "airtemp")# "nox", "amm")

    ## drop some predictors that are not relevant
    dlist <- c("kffact.pt", "superfunddensws", "superfunddenswsrp100",
               "tridenswsrp100", "pctconif2019ws",
               "caows", "fstfrz.pt", "pcthay2019ws", "pctfire2010ws",
               "elev.delt", "rhmean.pt", "pctshrb2019ws")
    for (i in dlist) {
        ip <- which(varscat == i)
        if (length(ip) > 0) varscat <- varscat[-ip]
        else print(i)
    }

    ## select taxon names from ss
    tnames <- names(ss)[regexpr("ASV", names(ss)) != -1]
    ss <- ss[, c("uid", tnames)]

    ## select mutable variables
    chemvar <- c("doc", "tss.result", "pct.safn", "cond2", "ptl.diss",
                 "nox","don", "anc.result")
    varall <- c(varscat, chemvar)

    ## set zero TSS to half min detection level and log transform
    incvec <- envdata$tss.result == 0
    envdata$tss.result[incvec] <- 0.05
    envdata$tss.result <- log(envdata$tss.result)

    ## set missing amm to zero. Not included as a separate variable
    ## so not worried about log transforming zero
    incvec <- is.na(envdata$amm)
    envdata$amm[incvec] <- 0
    print(summary(envdata$amm))

    ## calculating DON as the difference between DIN, NOx, and AMM
    envdata$don <- envdata$ntl.diss - envdata$nox - envdata$amm
    ## set don <= 0 to minval*0.5
    incvec <- envdata$don <= 0
    ## 33 estimates of don are less than or equal to zero
    ## set those to half the minimum positive value
    incvec[is.na(incvec)] <- F
    minval <- min(envdata$don[!incvec], na.rm = T)
    print(minval)
    envdata$don[incvec] <- 0.5*minval
    envdata$don <- log(envdata$don)

    incvec <- envdata$nox == 0
    envdata$nox[incvec] <- 0.1
    envdata$nox <- log(envdata$nox)
    #envdata$amm <- log(envdata$amm)
    envdata$doc <- log(envdata$doc)
    envdata$cond <- log(envdata$cond)
    envdata$ptl.diss <- log(envdata$ptl.diss)
    envdata$ntl.diss <- log(envdata$ntl.diss)
    print(summary(envdata$ntl.diss))
    envdata$wsareasqkm <- log(envdata$wsareasqkm)

    incvec <- envdata$anc.result <= 0
    incvec[is.na(incvec)] <- F
    print(sum(incvec))
    minval <- min(envdata$anc.result[!incvec], na.rm = T)
    envdata$anc.result[incvec] <- 0.5*minval
    envdata$anc.result <- log(envdata$anc.result)

    ## calculate cond adjusted for ANC
    envdata$cond2 <- exp(envdata$cond) - exp(-2.4)*exp(envdata$anc.result)

    ## select wadeable streams
    envdata <- subset(envdata, protocol =="WADEABLE")

    ## drop missing biological data
    incvec <- !is.na(ss[, tnames[1]])
    ss <- ss[incvec,]
    print(nrow(ss))

    ## drop missing envir data
    incvec <- rep(F, times = nrow(envdata))
    for (i in varall) {
        incvec <- incvec | is.na(envdata[, i])
    }
    envdata <- envdata[!incvec,]

    ## set negative cond2 to minimum positive value
    print(summary(envdata$cond2))
    incvec <- envdata$cond2 < 0
    minval <- 0.5*min(envdata$cond2[!incvec])
    envdata$cond2[incvec] <- minval
    envdata$cond2 <- log(envdata$cond2)

    ## drop repeat visits
    envdata <- envdata[!duplicated(envdata$site.id),]

#    write.table(envdata[, c("uid", varall)], sep = ",", row.names = F,
#                file = "envdata.csv")
#    write.table(ss[, c("uid", tnames)], sep = ",", row.names = F,
#                file = "diatoms.asv.csv")
#    stop()

    print(length(varall))
    print(sort(varall))

    if (runmod) {
        partout <- random.beta(ss[, c("uid", tnames)], envdata[, c("uid", varall)], sampid = "uid", nit = 500, runmod = runmod)
        return(partout)
    }
    else {
        varsel <- random.beta(ss[, c("uid", tnames)], envdata[, c("uid", varall)], sampid = "uid", nit = 2, runmod = runmod, partout = partout.beta)
        return(varsel)
    }

}

#partout.alpha <- selvars(ss.diatom, envdata1819, runmod = T)
#save(partout.alpha,file = "partout.alpha.rda")
varsel <- selvars(ss.diatom, envdata1819, runmod = F)
