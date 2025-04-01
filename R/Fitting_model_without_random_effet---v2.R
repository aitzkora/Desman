########################
### CLEANING SESSION ###
########################

rm(list = ls())

#################
### LIBRARIES ###
#################

library(dplyr) 
library(gtools)

###################
### DATA IMPORT ###
###################

Df1 <- read.csv2(file = "./data/survcalibfeb2024_1.csv")
Df2 <- read.csv2(file = "./data/survcalibaug2024_1.csv")
Df <- bind_rows(Df1, Df2)

Mat.Cov <- read.csv(file = "./data/Latrine_cov.csv", header = TRUE, sep = ";", dec = ",", row.names = 1)
nb.cov <- ncol(x = Mat.Cov) -1

Df <- full_join(Df, Mat.Cov, by = join_by(siteID))
id.col.cov <- (ncol(Df)-nb.cov+1):ncol(Df)

idx <- which(is.na(Df$durationg))
Df$durationg[idx] <- Inf

##########################
### FITTING ALL MODELS ###
### NO RANDOM EFFECT   ###
##########################

LogVraisemblance.NoCov <- function(par) {
  lalpha <- par[1]
  lbeta <- par[2]
  res <- vector(mode = "numeric", length = nrow(Df))
  for (i in 1:nrow(Df)) {
    aux1 <- pweibull(q = Df$durationd[i], scale = exp(lalpha), shape = exp(lbeta), lower.tail = FALSE)
    aux2 <- pweibull(q = Df$durationg[i], scale = exp(lalpha), shape = exp(lbeta), lower.tail = FALSE)
    res[i] <- aux1 - aux2
  }
  return(sum(log(res)))
}

LogVraisemblance.SingleCov <- function(par) {
  lalpha <- par[1]
  lbeta <- par[2]
  gamma <- par[3]
  res <- vector(mode = "numeric", length = nrow(Df))
  for (i in 1:nrow(Df)) {
    aux1 <- pweibull(q = Df$durationd[i], lower.tail = FALSE, shape = exp(lbeta), 
                     scale = exp(lalpha + gamma[1]*Z[i])) 
    aux2 <- pweibull(q = Df$durationg[i], lower.tail = FALSE, shape = exp(lbeta), 
                     scale = exp(lalpha + gamma[1]*Z[i])) 
    res[i] <- aux1 - aux2
  }
  return(sum(log(res)))
}

LogVraisemblance.MultCov <- function(par) {
  nb.cov <- ncol(Z)
  lalpha <- par[1]
  lbeta <- par[2]
  gamma <- par[3:(3+ncol(Z)-1)]
  res <- vector(mode = "numeric", length = nrow(Df))
  for (j in 1:nrow(Df)) {
    aux1 <- pweibull(q = Df$durationd[j], lower.tail = FALSE, shape = exp(lbeta), 
                     scale = exp(lalpha + sum(Z[j, ]*gamma))) 
    aux2 <- pweibull(q = Df$durationg[j], lower.tail = FALSE, shape = exp(lbeta), 
                     scale = exp(lalpha + sum(Z[j, ]*gamma))) 
    res[j] <- aux1 - aux2
  }
  return(sum(log(res)))
}


##########################
### FITTING ALL MODELS ###
### NO RANDOM EFFECT   ###
##########################

# Initial point for the iterative algo to maximize the likelihood
par.init.base <- c(mean((Df$durationd + Df$durationg)[-idx]/2), 1)
lpar.init.base <- log(par.init.base)

# Matrix of all possible sub-models
Mat.Sel.Var <- permutations(n = 2, r = nb.cov, v = c(TRUE, FALSE), repeats = TRUE)
nb.sous.model <- nrow(Mat.Sel.Var)

# Preparing the table for the results
Res <- matrix(nrow = nb.sous.model, ncol = 3+nb.cov)
Res <- as.data.frame(Res)
colnames(Res) <- c("AIC", "alpha", "beta", colnames(Mat.Cov)[-1])
rownames(Res) <- 1:nb.sous.model

# Fitting all sub-models
for (i in 1:nb.sous.model) {
  Sel.Var <- Mat.Sel.Var[i, ]
  Z <- Df[, id.col.cov[Sel.Var]]
  cat("*** MODEL WITH COVARIATES:", paste(colnames(Mat.Cov)[Sel.Var], collapse = " - "), "***", "\n")
  if (sum(Sel.Var) == 0) {
    lpar.init <- lpar.init.base
    fit <- optim(par =  lpar.init, fn = LogVraisemblance.NoCov, control = list(fnscale = -1), method = "L-BFGS-B")
    LogVraisemblance.NoCov(lpar.init)
    par <- fit$par
    aic <- -2*fit$value + 2*length(par)
    Res[i, 1] <- aic
    Res[i, 2:3] <- par[1:2]
  }
  if (sum(Sel.Var) == 1) {
    lpar.init <- c(lpar.init.base, 0)
    fit <- optim(par =  lpar.init, fn = LogVraisemblance.SingleCov, control = list(fnscale = -1), method = "L-BFGS-B")
    par <- fit$par
    aic <- -2*fit$value + 2*length(par)
    Res[i, 1] <- aic
    Res[i, 2:3] <- par[1:2]
    Res[i, (4:(3+nb.cov))[Sel.Var]] <- par[-(1:2)]
  }
  if (sum(Sel.Var) > 1) {
    lpar.init <- c(lpar.init.base, rep(x = 0, times = sum(Sel.Var)))
    fit <- optim(par =  lpar.init, fn = LogVraisemblance.MultCov, control = list(fnscale = -1), method = "L-BFGS-B")
    par <- fit$par
    aic <- -2*fit$value + 2*length(par)
    Res[i, 1] <- aic
    Res[i, 2:3] <- par[1:2]
    Res[i, (4:(3+nb.cov))[Sel.Var]] <- par[-(1:2)]
  }
}

ordre <- order(Res$AIC)
Res <- Res[ordre, ]
Res$alpha <- exp(Res$alpha)
Res$beta <- exp(Res$beta)

library(WriteXLS)
WriteXLS(x = Res, ExcelFileName = "Res_model_without_random_effet---v2.xlsx")

