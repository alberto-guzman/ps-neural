############## funzione: funsim (x, psmethod)

############## funsim uses psmethod to estimate ps; then use 1) weighting and 2) matching to estimate treatment effect. Finally, it calculates performance metrics.

# inputs: x = a simulated dataset created by funcov, a psmethod
# output: a list containing performance metrics (bias and var of estimate, covariate balance diagnostics)



funsim <- function(x, psmethod, par = "ATE") {

  # ~~ estimate ps


  if (psmethod == "truelogit") {
    mod <- glm(T ~ V1 + V2 + V3 + V4, data = x, family = binomial)
    ps <- mod$fitted
  } else if (psmethod == "logit") {
    mod <- glm(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = x, family = binomial)
    ps <- mod$fitted
  } else if (psmethod == "tree") {
    mod <- rpart(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, method = "class", data = x)
    ps <- predict(mod)[, 2]
  } else if (psmethod == "randomforest") {
    mod <- randomForest(factor(T) ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = x, ntree = 500)
    ps <- predict(mod, type = "prob")[, 2]
  } else if (psmethod == "gbm") {
    mod <- gbm(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = x, distribution = "bernoulli", interaction.depth = 1, n.trees = 100)
    ps <- predict.gbm(mod, data = x, n.trees = 100, type = "response")
  } else if (psmethod == "gbmtwang") {
    mod <- ps(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = x, n.trees = 10000, interaction.depth = 3, verbose = FALSE, shrinkage = 0.0005)
    ps <- as.vector(mod$ps[, 1])
  } else if (psmethod == "bag") {
    mod <- bagging(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = x)
    ps <- predict(mod, newdata = x, type = "prob")
  } else if (psmethod == "nnet") {
    mod <- nnet(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10,
      data = x, entropy = T, size = 10, decay = 0, maxit = 2000, trace = F
    )
    ps <- as.numeric(predict(mod, type = "raw")) # nb: anche predizioni fuori da [0,1]
    #   ps=exp(ps)/(1+exp(ps))
  } else if (psmethod == "nb") {
    mod <- naiveBayes(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = x)
    ps <- predict(mod, newdata = x, type = "raw")[, 2]
  }

  #### fit measure for ps

  auc <- auc(sensitivity(ps, factor(x$T), perc.rank = TRUE))

  #### true (sample) att:

  g <- mean(x$indeff[x$T == 1])

  ### true (sample) ate

  if (par == "ATE") {
    g <- mean(x$indeff)
  }

  #### estimating ATT via propensity score weighting

  weights <- ifelse(x$T == 1, 1, ps / (1 - ps))

  #### estimating ATE via propensity score weighting

  if (par == "ATE") {
    weights <- ifelse(x$T == 1, 1 / ps, 1 / (1 - ps))
  }

  hatg <- wtd.mean(x$Y[x$T == 1], weights = weights[x$T == 1]) - wtd.mean(x$Y[x$T == 0], weights = weights[x$T == 0])

  absrbias <- abs((hatg - g) / g) * 100
  varhatg <- (hatg - g)^2

  modw <- lm(Y ~ T, data = x, weights = weights)
  hatgsew <- summary(modw)$coefficients[c("T"), c("Std. Error")]
  covw <- ifelse(g > hatg - 2 * hatgsew & g < hatg + 2 * hatgsew, 1, 0)



  # estimating ATT via matching with caliper

  rr <- Match(Y = x$Y, Tr = x$T, X = ps, caliper = 0.25, M = 1, replace = TRUE, ties = FALSE)

  if (par == "ATE") {
    rr <- Match(Y = x$Y, Tr = x$T, X = ps, caliper = 0.25, M = 1, estimand = par, replace = TRUE, ties = FALSE)
  }

  # estimating ATE via matching with caliper

  allmdata <- rbind(x[rr$index.treated, ], x[rr$index.control, ])

  hatgm <- rr$est
  absrbiasm <- abs((hatgm - g) / g) * 100
  varhatgm <- (hatgm - g)^2

  hatgsem <- rr$se.standard
  covm <- ifelse(g > hatgm - 2 * hatgsem & g < hatgm + 2 * hatgsem, 1, 0)
  # = mean(mdata$Y[mdata$Tr==1])-mean(mdata$Y[mdata$Tr==0])

  # ~  size of matched dataset

  orig.nobs <- rr$orig.nobs
  orig.tnobs <- rr$orig.treated.nobs
  # match.nobs   <-rr$nobs # this is wrong! it is the nobs inthe unmatched dataset
  match.nobs <- length(rr$mdata$Tr)
  match.tnobs <- length(rr$mdata$Tr[rr$mdata$Tr == 1])




  # Balance BEFORE for V1 ... V10

  bb <- MatchBalance(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10
    + V1 * V2 + V1 * V3 + V1 * V4 + V1 * V5 + V1 * V6 + V1 * V7 + V1 * V8 + V1 * V9 + V1 * V10
    + V2 * V3 + V2 * V4 + V2 * V5 + V2 * V6 + V2 * V7 + V2 * V8 + V2 * V9 + V2 * V10
    + V3 * V4 + V3 * V5 + V3 * V6 + V3 * V7 + V3 * V8 + V3 * V9 + V3 * V10
    + V4 * V5 + V4 * V6 + V4 * V7 + V4 * V8 + V4 * V9 + V4 * V10
    + V5 * V6 + V5 * V7 + V5 * V8 + V5 * V9 + V5 * V10
    + V6 * V7 + V6 * V8 + V6 * V9 + V6 * V10
    + V7 * V8 + V7 * V9 + V7 * V10
    + V8 * V9 + V8 * V10
    + V9 * V10, data = x, ks = FALSE, nboots = 0, print.level = 0)
  # e.g.


  # ASAM
  asb_b <- vector()
  for (i in 1:10) {
    asb_b[[i]] <- bb$BeforeMatching[[i]]$sdiff
  }

  # ASAM with interactions
  interasbb <- vector()
  for (i in 1:(length(bb$BeforeMatching))) {
    interasbb[[i]] <- bb$BeforeMatching[[i]]$sdiff
  }

  # qq mean difference raw
  vecqqmeanrawb <- vector()
  for (i in 1:10) {
    vecqqmeanrawb[[i]] <- bb$BeforeMatching[[i]]$qqsummary.raw$meandiff
  }

  # qq mean diff
  vecqqmeanb <- vector()
  for (i in 1:10) {
    vecqqmeanb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff
  }
  interqqb <- vector()
  for (i in 1:length(bb$BeforeMatching)) {
    interqqb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$meandiff
  }

  # qq max difference
  vecqqmaxb <- vector()
  for (i in 1:10) {
    vecqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff
  }

  # qq max with difference
  interqqmaxb <- vector()
  for (i in 1:(length(bb$BeforeMatching))) {
    interqqmaxb[[i]] <- bb$BeforeMatching[[i]]$qqsummary$maxdiff
  }

  # variance ratio
  # i.e. before log (var ratio) for ps; reference: Imbens-Rubin Ch 10, eq. 14.5 and tobacco litigation
  varratio_b <- abs(log(var(ps[x$T == 1]) / var(ps[x$T == 0])))

  # ~~ performance metrics: balance AFTER weighting for w1 ...w10

  ba <- MatchBalance(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10
    + V1 * V2 + V1 * V3 + V1 * V4 + V1 * V5 + V1 * V6 + V1 * V7 + V1 * V8 + V1 * V9 + V1 * V10
    + V2 * V3 + V2 * V4 + V2 * V5 + V2 * V6 + V2 * V7 + V2 * V8 + V2 * V9 + V2 * V10
    + V3 * V4 + V3 * V5 + V3 * V6 + V3 * V7 + V3 * V8 + V3 * V9 + V3 * V10
    + V4 * V5 + V4 * V6 + V4 * V7 + V4 * V8 + V4 * V9 + V4 * V10
    + V5 * V6 + V5 * V7 + V5 * V8 + V5 * V9 + V5 * V10
    + V6 * V7 + V6 * V8 + V6 * V9 + V6 * V10
    + V7 * V8 + V7 * V9 + V7 * V10
    + V8 * V9 + V8 * V10
    + V9 * V10, data = x, weights = weights, ks = FALSE, nboots = 0, print.level = 0)

  # ASAM
  asb_a <- vector()
  for (i in 1:10) {
    asb_a[[i]] <- ba$BeforeMatching[[i]]$sdiff
  }

  # ASAM with interactions
  interasba <- vector()
  for (i in 1:(length(ba$BeforeMatching))) {
    interasba[[i]] <- ba$BeforeMatching[[i]]$sdiff
  }


  # QQ mean raw
  # # mean difference of empirical quantiles for w1 with weights applied
  # # (cannot use weighths with matchBalance forr qq measure)
  covnames <- colnames(x[, -which(names(x) %in% c("Y", "T"))])
  vecqqmeanrawaw <- vector()
  for (i in 1:length(covnames)) {
    vecqqmeanrawaw[i] <- abs(
      mean(wtd.quantile(x[x$T == 1, ] [, covnames[i]], probs = seq(0.01, 1 - 0.01, 0.01), weights = weights[x$T == 1]) - mean(wtd.quantile(x[x$T == 0, ] [, covnames[i]], probs = seq(0.01, 1 - 0.01, 0.01), weights = weights[x$T == 0])))
    )
  }

  # # equivalent to:
  # W1qqrawa=abs(mean(wtd.quantile(x[x$T==1,]$w1,probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==1])-mean(wtd.quantile(x[x$T==0,]$w1,probs=seq(0.01,1-0.01,0.01),weights=weights[x$T==0]))))

  # mean and max differences of ecdf after weighting (nb: non si possono usare pesi con matchbalance per qq)
  covnames <- colnames(x[, -which(names(x) %in% c("Y", "T"))])
  vecqqmeanaw <- vecqqmaxaw <- vector()
  for (i in 1:length(covnames)) {
    vals <- sort(unique(x[, covnames[i]]))
    wt <- stepfun(wtd.Ecdf(x[x$T == 1, ][, covnames[i]], weights = weights[x$T == 1])$x[-1], wtd.Ecdf(x[x$T == 1, ][, covnames[i]], weights = weights[x$T == 1])$ecdf)
    swt <- wt(vals)
    wc <- stepfun(wtd.Ecdf(x[x$T == 0, ][, covnames[i]], weights = weights[x$T == 0])$x[-1], wtd.Ecdf(x[x$T == 0, ][, covnames[i]], weights = weights[x$T == 0])$ecdf)
    swc <- wc(vals)
    vecqqmeanaw[i] <- mean(abs(swt - swc))
    vecqqmaxaw[i] <- max(abs(swt - swc))
  }


  # variance ratio
  varratio_a <- abs(log(wtd.var(ps[x$T == 1], weights = weights[x$T == 1]) / wtd.var(ps[x$T == 0], weights = weights[x$T == 0]))) # log (var ratio) after weighting; reference in Inbens-Rubin Ch 10, eq. 14.5 and tobacco litigation




  # ~~ performance metrics: balance AFTER matching for w1 ...w10

  bam <- MatchBalance(T ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10
    + V1 * V2 + V1 * V3 + V1 * V4 + V1 * V5 + V1 * V6 + V1 * V7 + V1 * V8 + V1 * V9 + V1 * V10
    + V2 * V3 + V2 * V4 + V2 * V5 + V2 * V6 + V2 * V7 + V2 * V8 + V2 * V9 + V2 * V10
    + V3 * V4 + V3 * V5 + V3 * V6 + V3 * V7 + V3 * V8 + V3 * V9 + V3 * V10
    + V4 * V5 + V4 * V6 + V4 * V7 + V4 * V8 + V4 * V9 + V4 * V10
    + V5 * V6 + V5 * V7 + V5 * V8 + V5 * V9 + V5 * V10
    + V6 * V7 + V6 * V8 + V6 * V9 + V6 * V10
    + V7 * V8 + V7 * V9 + V7 * V10
    + V8 * V9 + V8 * V10
    + V9 * V10, data = allmdata, ks = FALSE, nboots = 0, print.level = 0)

  # ASAM after matching
  asb_am <- vector()
  for (i in 1:10) {
    asb_am[[i]] <- bam$BeforeMatching[[i]]$sdiff
  }

  # ASAM with interactions after matching
  interasbam <- vector()
  for (i in 1:(length(bam$BeforeMatching))) {
    interasbam[[i]] <- bam$BeforeMatching[[i]]$sdiff
  }

  vecqqmeanrawam <- vector()
  for (i in 1:10) {
    vecqqmeanrawam[[i]] <- bam$BeforeMatching[[i]]$qqsummary.raw$meandiff
  }

  # mean diff in ecdf after matching
  vecqqmeanam <- vector()
  for (i in 1:10) {
    vecqqmeanam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff
  }

  interqqam <- vector()
  for (i in 1:(length(bam$BeforeMatching))) {
    interqqam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$meandiff
  }

  # standardized max difference of ecdfs after matching
  vecqqmaxam <- vector()
  for (i in 1:10) {
    vecqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff
  }

  interqqmaxam <- vector()
  for (i in 1:(length(bam$BeforeMatching))) {
    interqqmaxam[[i]] <- bam$BeforeMatching[[i]]$qqsummary$maxdiff
  }

  # variance ratio
  varratio_am <- abs(log(var(ps[rr$index.treated]) / var(ps[rr$index.control]))) # abs log (var ratio) after matching

  output <- list(
    "auc" = auc,
    "hatgw" = hatg, "absrbiasw" = absrbias, "varhatgw" = varhatg,
    "hatgm" = hatgm, "absrbiasm" = absrbiasm, "varhatgm" = varhatgm,
    "hatgsem" = hatgsem, "hatgsew" = hatgsew,
    "covw" = covw, "covm" = covm,
    #  mean (over covariates) of balance summaries :
    # asam
    "masb_b" = mean(abs(asb_b)),
    "masb_aw" = mean(abs(asb_a)),
    "masb_am" = mean(abs(asb_am)),
    # mean asam>20
    "over20masb_b" = sum(abs(asb_b[which(abs(asb_b) > 20)]) - 20) / length(asb_b),
    "over20masb_aw" = sum(abs(asb_a[which(abs(asb_a) > 20)]) - 20) / length(asb_a),
    "over20masb_am" = sum(abs(asb_am[which(abs(asb_am) > 20)]) - 20) / length(asb_am),
    # mean asam>20
    "over10masb_b" = sum(abs(asb_b[which(abs(asb_b) > 10)]) - 10) / length(asb_b),
    "over10masb_aw" = sum(abs(asb_a[which(abs(asb_a) > 10)]) - 10) / length(asb_a),
    "over10masb_am" = sum(abs(asb_am[which(abs(asb_am) > 10)]) - 10) / length(asb_am),
    # asam with interactions
    "masbinter_b" = mean(abs(interasbb)),
    "masbinter_aw" = mean(abs(interasba)),
    "masbinter_am" = mean(abs(interasbam)),
    #
    "qqmeanraw_b" = mean(vecqqmeanrawb),
    "qqmeanraw_aw" = mean(vecqqmeanrawaw),
    "qqmeanraw_am" = mean(vecqqmeanam),
    #
    "qqmean_b" = mean(vecqqmeanb),
    "qqmean_aw" = mean(vecqqmeanaw),
    "qqmean_am" = mean(vecqqmeanam),
    # "qqmeaninter_b"=mean(abs(interqqb)),"qqmeaninter_am"=mean(abs(interqqam)),
    "qqmax_b" = mean(vecqqmaxb),
    "qqmax_aw" = mean(vecqqmaxaw),
    "qqmax_am" = mean(vecqqmaxam),
    # variance ratio
    "varratio_b" = varratio_b,
    "varratio_aw" = varratio_a,
    "varratio_am" = varratio_am,
    #
    "match.tnobs" = match.tnobs, "match.nobs" = match.nobs
  )
  return(output)
}



funsim(data, "truelogit", par = "ATT")