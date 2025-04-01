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
nb.cov <- ncol(x = Mat.Cov) -1 # Why ? 

Df <- full_join(Df, Mat.Cov, by = join_by(siteID))
id.col.cov <- (ncol(Df)-nb.cov+1):ncol(Df)

idx <- which(is.na(Df$durationg))
Df$durationg[idx] <- Inf

#
#
###########################
#### FITTING ALL MODELS ###
#### NO RANDOM EFFECT   ###
###########################
#
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
#
#
###########################
#### FITTING ALL MODELS ###
#### RANDOM EFFECT      ###
###########################
#
Shared.LogVraisemblance.NoCov <- function(par) {
  lalpha <- par[1]
  lbeta <- par[2]
  ltheta <- par[3]
  res <- vector(mode = "numeric", length = n)
  for (i in 1:n) {
    idx <- which(Df$siteID == Latrine[i])
    Df.tmp <- Df[idx, ]
    faux <- function(w) {
      tmp <- vector(mode = "numeric", length = m[i])
      for (j in 1:m[i]) {
        aux1 <- pweibull(q = Df.tmp$durationd[j], scale = w*exp(lalpha), shape = exp(lbeta), lower.tail = FALSE)
        aux2 <- pweibull(q = Df.tmp$durationg[j], scale = w*exp(lalpha), shape = exp(lbeta), lower.tail = FALSE)
        tmp[j] <- aux1 - aux2
      }
      aux <- prod(tmp)*dgamma(x = w, shape = exp(ltheta), scale = 1/exp(ltheta))
      return(aux)
    }
    res[i] <- integrate(f = function(w){sapply(w,faux)}, lower = 0, upper = Inf)$value
    #cat(i, " ", res[i], "\n")
  }
  return(sum(log(res)))
}
#
#Shared.LogVraisemblance.SingleCov <- function(par) {
#  lalpha <- par[1]
#  lbeta <- par[2]
#  ltheta <- par[3]
#  gamma <- par[4]
#  res <- vector(mode = "numeric", length = n)
#  for (i in 1:n) {
#    idx <- which(Df$siteID == Latrine[i])
#    Df.tmp <- Df[idx, ]
#    faux <- function(w) {
#      tmp <- vector(mode = "numeric", length = m[i])
#      for (j in 1:m[i]) {
#        aux1 <- pweibull(q = Df.tmp$durationd[j], lower.tail = FALSE, shape = exp(lbeta), 
#                         scale = w*exp(lalpha + gamma[1]*Z[j])) 
#        aux2 <- pweibull(q = Df.tmp$durationg[j], lower.tail = FALSE, shape = exp(lbeta), 
#                         scale = w*exp(lalpha + gamma[1]*Z[j])) 
#        tmp[j] <- aux1 - aux2
#      }
#      aux <- prod(tmp)*dgamma(x = w, shape = exp(ltheta), scale = 1/exp(ltheta))
#      return(aux)
#    }
#    res[i] <- integrate(f = function(w){sapply(w,faux)}, lower = 0, upper = Inf)$value
#  }
#  return(sum(log(res)))
#}
#
Shared.LogVraisemblance.MultCov <- function(par) {
  nb.cov <- ncol(Z)
  lalpha <- par[1]
  lbeta <- par[2]
  ltheta <- par[3]
  gamma <- par[4:(4+ncol(Z)-1)]
  res <- vector(mode = "numeric", length = n)
  for (i in 1:n) {
    idx <- which(Df$siteID == Latrine[i])
    Df.tmp <- Df[idx, ]
    faux <- function(w) {
      tmp <- vector(mode = "numeric", length = m[i])
      for (j in 1:m[i]) {
        aux1 <- pweibull(q = Df.tmp$durationd[j], lower.tail = FALSE, shape = exp(lbeta), scale = w*exp(lalpha + sum(Z[j, ]*gamma))) 
        aux2 <- pweibull(q = Df.tmp$durationg[j], lower.tail = FALSE, shape = exp(lbeta), scale = w*exp(lalpha + sum(Z[j, ]*gamma))) 
        tmp[j] <- aux1 - aux2
      }
      aux <- prod(tmp)*dgamma(x = w, shape = exp(ltheta), scale = 1/exp(ltheta))
      return(aux)
    }
    res[i] <- integrate(f = function(w){sapply(w,faux)}, lower = 0, upper = Inf)$value
  }
  return(sum(log(res)))
}

#
#
#
#
###########################
#### FITTING ALL MODELS ###
#### RANDOM EFFECT      ###
###########################
#
# Initial point for the iterative algo to maximize the likelihood
par.init.base <- c(mean((Df$durationd + Df$durationg)[-idx]/2), 1, 1)

# Occupency of the lattrines
Latrine <- sort(unique(Df$siteID))
n <- length(Latrine)
m <- vector(mode = "integer", length = n)
for (i in 1:n) {
  m[i] <- length(which(Df$siteID == Latrine[i]))
}

## Matrix of all possible sub-models
Mat.Sel.Var <- permutations(n = 2, r = nb.cov, v = c(TRUE, FALSE), repeats = TRUE)
nb.sous.model <- nrow(Mat.Sel.Var)

# Preparing the table for the results
Res <- matrix(nrow = nb.sous.model, ncol = 4+nb.cov)
Res <- as.data.frame(Res)
colnames(Res) <- c("AIC", "alpha", "beta", colnames(Mat.Cov)[-1], "theta")
rownames(Res) <- 1:nb.sous.model
Shared.LogVraisemblance.NoCov(par.init.base)
# Fitting all sub-models
for (i in 1:nb.sous.model) {
  Sel.Var <- Mat.Sel.Var[i, ]
  Z <- Df[, id.col.cov[Sel.Var]]
  cat(id.col.cov[Sel.Var], "\n")
  cat("*** MODEL WITH COVARIATES:", paste(colnames(Mat.Cov[-1])[Sel.Var], collapse = " - "), "***", "\n")
  #if (sum(Sel.Var) == 0) {
  #  par.init <- log(par.init.base)
  #  fit <- try(optim(par =  par.init, fn = Shared.LogVraisemblance.NoCov, control = list(fnscale = -1), method = "L-BFGS-B"))
  #  par <- exp(fit$par)
  #  aic <- -2*fit$value + 2*length(par)
  #  Res[i, 1] <- aic
  #  Res[i, 2:3] <- par[1:2]
  #  Res[i, 4+nb.cov] <- par[3]
  #}
  #if (sum(Sel.Var) == 1) {
  #  par.init <- c(par.init.base, 0)
  #  fit <- try(optim(par =  par.init, fn = Shared.LogVraisemblance.SingleCov, control = list(fnscale = -1)))
  #  par <- fit$par
  #  aic <- -2*fit$value + 2*length(par)
  #  Res[i, 1] <- aic
  #  Res[i, 2:3] <- par[1:2]
  #  Res[i, (4:(3+nb.cov))[Sel.Var]] <- par[3]
  #  Res[i, 4+nb.cov] <- par[4]
  #}
  if (sum(Sel.Var) > 1) {
    par.init <- c(par.init.base, rep(x = 0, times = sum(Sel.Var)))
    #fit <- try(optim(par = par.init, fn = Shared.LogVraisemblance.MultCov, control = list(fnscale = -1), method = "L-BFGS-B"))
    #par <- fit$par
    
    #aic <- -2*fit$value + 2*length(par)
    #Res[i, 1] <- aic
    #Res[i, 2:3] <- par[1:2]
    #Res[i, (4:(3+nb.cov))[Sel.Var]] <- par[-(1:2)]
  }
  ordre <- order(Res$AIC)
  Res <- Res[ordre, ]
  
  #library(WriteXLS)
  #WriteXLS(x = Res, ExcelFileName = "Res_model_with_random_effet.xlsx")
  
}

#
