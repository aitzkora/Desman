Cov.Mat <- read.csv(file = "./data/Latrine_cov.csv", header = TRUE, sep = ";", dec = ",")
nb.sites <- nrow(Cov.Mat)

Case.Selected <- 1
if (Case.Selected == 1) {
  parms <- c(3.77, 0.55, 0.47, 0.001, 18.36, 1.07, 1.14, 1.03)
}

set.seed(123)

Durées.Exactes <- vector(mode = "list", length = nrow(Cov.Mat))
for (i in 1:nb.sites) {
  Si <- rgamma(n = 1, shape = parms[3], scale = 1/parms[3])
  tps <- 0
  j <- 0
  Tps <- NULL
  tmax <- 25
  s <- Si*parms[3]*sum(Cov.Mat[i, -(1:2)]*parms[4:8])
  while (tps < tmax) {
    j <- j+1
    Vj <- rweibull(n = 1, scale = parms[1], shape = parms[2])
    tps <- tps+Vj
    Tps <- c(Tps, tps)
  }
  Durées.Exactes[[i]] <- Tps
}

# Nb.Occ <- unlist(lapply(Durées.Exactes, length))
# Level.Occ.tmp <- cut(Nb.Occ, breaks = c(0,5, 14))
# Level.Occ <- vector(mode = "character", length = nb.sites)
# Level.Occ[Level.Occ.tmp == "(0,5]"] <- "Low.Occ"
# Level.Occ[Level.Occ.tmp == "(5,14]"] <- "High.Occ"
# Level.Occ <- ordered(Level.Occ, levels = c("Low.Occ", "High.Occ"))
# 
# Cov.Mat.splitted <- split(x = Cov.Mat[, -(1:2)], f = Level.Occ)
# lapply(Cov.Mat.splitted, summary)
# for (i in 1:5) {
#   wtest <- wilcox.test(x = Cov.Mat.splitted[[1]][, i], y = Cov.Mat.splitted[[2]][, i])
#   print(wtest)
# }

delta <- 2
Df2 <- NULL
for (i in 1:nb.sites) {
  Tij <- Durées.Exactes[[i]]
  Xij <- delta*ceiling(Tij/delta)
  Xij <- unique(Xij)
  nb <- length(Xij)
  # if (nb == 2) {
  #   durationd <- c(0, Xij[1])
  #   durationg <- c(Xij[1], Inf)
  #   Df.tmp <- data.frame(Cov.Mat[i, ], durationg, durationd)
  #   Df2 <- rbind.data.frame(Df2, Df.tmp)
  # }
  if (nb >= 2) {
    durationd <- c(0, Xij[1:(nb-1)])
    durationg <- c(Xij[1:(nb-1)], Inf)
    Df.tmp <- data.frame(Cov.Mat[i, ], durationg, durationd)
    Df2 <- rbind.data.frame(Df2, Df.tmp)
  }
}

File.Name <- paste("Simulation-Case-",Case.Selected, ".csv", sep = "")
write.csv(x = Df2, file = File.Name)

