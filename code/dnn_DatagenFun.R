datagen <- function(n, scenarioT, scenarioY) {
  #####################
  # DATA GENERATION CODE
  #####################

  # creating matrix
  set.seed(5)
  mat <- rcorrmatrix(10)

  # test to check it is positive definite
  is.positive.definite(mat)

  # generating mixed data
  # n <- 50
  num_pois <- 2
  num_bin <- 1
  num_ord <- 2
  num_norm <- 5
  lamvec <- sample(10, 2)
  pbin <- runif(1)
  pord <- list(c(0.3, 0.7), c(0.2, 0.3, 0.5))
  nor.mean <- c(3.1, 4, 12, 3, 8)
  nor.var <- c(0.85, .9, .2, .6, .8)


  intmat <- intermat(num_pois, num_bin, num_ord, num_norm, corr_mat = mat, pbin, pord, lamvec, nor.mean, nor.var)

  # gen the data
  dat <- genPBONdata(n, num_pois, num_bin, num_ord, num_norm, intmat, lamvec, pbin, pord, nor.mean, nor.var)

  # gen dataframe
  dat <- as.data.frame(dat$data)
  summary(dat)
  cor(dat, use = "complete.obs", method = "pearson")

  # ~~ scenarios for data generation models
  # A: model with additivity and linearity
  # B: mild non-linearity
  # C: moderate non-linearity
  # D: mild non-additivity
  # E: mild non-additivity and non-linearity
  # F: moderate non-additivity
  # G: moderate non-additivity and non-linearity
  # binary exposure modeling

  # ~~~~~~~~~~~~~~Global Variables~~~~~~~~~~~~~~~~~~~~#

  # ~~ coefficients for the treatment and potential outcomes equations

  # ~~ coefficients for treatment equation (from Setoguchi (reference 1))
  b0 <- 0
  b1 <- 0.8
  b2 <- -0.25
  b3 <- 0.6
  b4 <- -0.4
  b5 <- -0.8
  b6 <- -0.5
  b7 <- 0.7



  # ~~ coefficients for potential outcomes equation Y(0) and Y(1)
  # ~~ coefficients for Y(0): from Setoguchi (reference 1)
  # ~~ coefficients for Y(1): same as above plus delta. In particular
  # ~~ the intercept in Y(1) has been fixed to obtain Y(1)-Y(0)=-0.4 in each scenario
  a00 <- -3.85
  a01 <- 0.3
  a02 <- -0.36
  a03 <- -0.73
  a04 <- -0.2
  a05 <- 0.71
  a06 <- -0.19
  a07 <- 0.26


  # ~~ scenarios for treatment assignment
  attach(dat)

  if (scenarioT == "A") {
    trueps <- (1 + exp(-(0 + b1 * V1 + b2 * V2 + b3 * V3 + b4 * V4 + b5 * V5 + b6 * V6 + b7 * V7)))^-1
  } else
  if (scenarioT == "B") {
    trueps <- (1 + exp(-(0 + b1 * V1 + b2 * V2 + b3 * V3 + b4 * V4 + b5 * V5 + b6 * V6 + b7 * V7
      + b2 * V2 * V2)))^-1
  } else
  if (scenarioT == "C") {
    trueps <- (1 + exp(-(0 + b1 * V1 + b2 * V2 + b3 * V3 + b4 * V4 + b5 * V5 + b6 * V6 + b7 * V7
      + b2 * V2 * V2 + b4 * V4 * V4 + b7 * V7 * V7)))^-1
  } else
  if (scenarioT == "D") {
    trueps <- (1 + exp(-(0 + b1 * V1 + b2 * V2 + b3 * V3 + b4 * V4 + b5 * V5 + b6 * V6 + b7 * V7
      + b1 * 0.5 * V1 * V3 + b2 * 0.7 * V2 * V4 + b4 * 0.5 * V4 * V5 + b5 * 0.5 * V5 * V6)))^-1
  } else
  if (scenarioT == "E") {
    trueps <- (1 + exp(-(0 + b1 * V1 + b2 * V2 + b3 * V3 + b4 * V4 + b5 * V5 + b6 * V6 + b7 * V7
      + b2 * V2 * V2 + b1 * 0.5 * V1 * V3 + b2 * 0.7 * V2 * V4 + b4 * 0.5 * V4 * V5 + b5 * 0.5 * V5 * V6)))^-1
  } else
  if (scenarioT == "F") {
    trueps <- (1 + exp(-(0 + b1 * V1 + b2 * V2 + b3 * V3 + b4 * V4 + b5 * V5 + b6 * V6 + b7 * V7
      + b1 * 0.5 * V1 * V3 + b2 * 0.7 * V2 * V4 + b3 * 0.5 * V3 * V5 + b4 * 0.7 * V4 * V6 + b5 * 0.5 * V5 * V7
      + b1 * 0.5 * V1 * V6 + b2 * 0.7 * V2 * V3 + b3 * 0.5 * V3 * V4 + b4 * 0.5 * V4 * V5 + b5 * 0.5 * V5 * V6)))^-1
  } else {
    # scenario G
    trueps <- (1 + exp(-(0 + b1 * V1 + b2 * V2 + b3 * V3 + b4 * V4 + b5 * V5 + b6 * V6 + b7 * V7
      + b2 * V2 * V2 + b4 * V4 * V4 + b7 * V7 * V7 + b1 * 0.5 * V1 * V3 + b2 * 0.7 * V2 * V4 + b3 * 0.5 * V3 * V5
      + b4 * 0.7 * V4 * V6 + b5 * 0.5 * V5 * V7 + b1 * 0.5 * V1 * V6 + b2 * 0.7 * V2 * V3 + b3 * 0.5 * V3 * V4
      + b4 * 0.5 * V4 * V5 + b5 * 0.5 * V5 * V6)))^-1
  }

  # ~~ binary treatment T
  unif1 <- runif(n, 0, 1)
  T <- ifelse(trueps > unif1, 1, 0) # there is a probability of unif1 that T=1

  # ~~ scenarios for outcome

  if (scenarioY == "a") {
    Y0 <- a00 + a01 * V1 + a02 * V2 + a03 * V3 + a04 * V4 +
      a05 * V8 + a06 * V9 + a07 * V10
    Y1 <- (a00 - 0.4) + (a01 + 0) * V1 + (a02 + 0) * V2 + (a03 + 0) * V3 + (a04 + 0) * V4 + (a05 + 0) * V8 + (a06 + 0) * V9 + (a07 + 0) * V10
  } else
  if (scenarioY == "d") { # 1 inter
    Y0 <- a00 + a01 * V1 + a02 * V2 + a03 * V3 + a04 * V4 +
      a05 * V8 + a06 * V9 + a07 * V10
    Y1 <- (a00 - 0.4) + (a01 + 0) * V1 + (a02 + 0.25 * a02) * V2 + (a03 + 0) * V3 + (a04 + 0) * V4 + (a05 + 0) * V8 + (a06 + 0) * V9 + (a07 + 0) * V10
  } else
  if (scenarioY == "f") { # 2 inter
    Y0 <- a00 + a01 * V1 + a02 * V2 + a03 * V3 + a04 * V4 +
      a05 * V8 + a06 * V9 + a07 * V10
    Y1 <- (a00 - 0.4) + (a01 + 0) * V1 + (a02 + 0.25 * a02) * V2 + (a03 + 0) * V3 + (a04 + 0.25 * a04) * V4 + (a05 + 0) * V8 + (a06 + 0) * V9 + (a07 + 0) * V10
  } else
  if (scenarioY == "g") { # (2 quad, 2 inter)
    Y0 <- a00 + a01 * V1 + a02 * V2 + a03 * V3 + a04 * V4 +
      a05 * V8 + a06 * V9 + a07 * V10
    +0.1 * a02 * V2^2 + 0.1 * a04 * V4^2
    Y1 <- (a00 - 0.4) + (a01 + 0) * V1 + (a03 + 0) * V3 + (a05 + 0) * V8 + (a06 + 0) * V9 + (a07 + 0) * V10
    +0.25 * 0.1 * a02 * V2^2 + 0.25 * 0.1 * a04 * V4^2
    +(a02 + 0.25 * a02) * V2 + (a04 + 0.25 * a04) * V4
  } else
  if (scenarioY == "b") { # (1 quad terms)
    Y0 <- a00 + a01 * V1 + a02 * V2 + a03 * V3 + a04 * V4 +
      a05 * V8 + a06 * V9 + a07 * V10
    +0.5 * a02 * V2^2
    Y1 <- (a00 - 0.2) + (a01 + 0) * V1 + (a02 + 0) * V2 + (a03 + 0) * V3 + (a04 + 0) * V4 + (a05 + 0) * V8 + (a06 + 0) * V9 + (a07 + 0) * V10
    +a02 * V2^2
  } else
  if (scenarioY == "c") { # (2 quad terms)
    Y0 <- a00 + a01 * V1 + a02 * V2 + a03 * V3 + a04 * V4 +
      a05 * V8 + a06 * V9 + a07 * V10
    +0.5 * a02 * V2^2 + 0.5 * a04 * V4^2
    Y1 <- (a00 - 0.1) + (a01 + 0) * V1 + (a02 + 0) * V2 + (a03 + 0) * V3 + (a04 + 0) * V4 + (a05 + 0) * V8 + (a06 + 0) * V9 + (a07 + 0) * V10
    +a02 * V2^2 + a04 * V4^2
  } else
  if (scenarioY == "e") {
    # scenario e (1 quad term + 1 interaction)
    Y0 <- a00 + a01 * V1 + a02 * V2 + a03 * V3 + a04 * V4 +
      a05 * V8 + a06 * V9 + a07 * V10
    +0.5 * a02 * V2^2
    Y1 <- (a00 - 0.2) + (a01 + a01) * V1 + (a02 + a02) * V2 + (a03 + a03) * V3 + (a04 + a04) * V4 + (a05 + 0) * V8 + (a06 + 0) * V9 + (a07 + 0) * V10
    +a02 * V2^2 + (a04 + 0.25 * a04) * V4
  }

  # continuous outcome Y

  Y <- T * Y1 + (1 - T) * Y0

  # true individual effect of T on Y

  indeff <- Y1 - Y0


  # create simulation dataset

  sim <- as.data.frame(cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, T, Y, indeff))
  return(sim)
}
# will save out a dataframe data so I can impliment the other thing in keras
data <- datagen(500, "G", "e")
#write.csv(data, "data_for_keras.csv", row.names = FALSE)




#group agregation 
aggregate(x = x$indeff,                # Specify data column
          by = list(x$T),              # Specify group indicator
          FUN = mean)                           # Specify function (i.e. mean)