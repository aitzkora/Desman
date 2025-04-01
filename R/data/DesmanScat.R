rm(list = ls())

library(readxl)
Desman <- read_excel("Desman.xlsx", na = "NA")
Desman <- Desman[-c(16, 22, 37), ]
Dates <- as.Date(x = colnames(Desman)[-1])
Nb.Dates.Obs <- length(Dates)
Nb.Latrines <- nrow(Desman)
  
Df <- NULL

for (i in 1:Nb.Latrines) {
  idx <- which(Desman[i, ] == 1)
  r <- length(idx)
  date.0 <- Dates[1]
  left <- vector(mode = "numeric", length = r)
  right <- vector(mode = "numeric", length = r)
  for (k in 1:r) {
    date.left <- Dates[idx[k]-2]
    date.right <- Dates[idx[k]-1]
    left[k] <- as.numeric(difftime(time1 = date.left, time2 = date.0))
    right[k] <- as.numeric(difftime(time1 = date.right, time2 = date.0))
    # Desman[i, idx[k]] <- 0
    date.0 <- Dates[idx[k]-1]
  }
  # if (Desman[i, Nb.Dates.Obs+1] == 0) {
  #   left <- NA
  #   right <- as.numeric(difftime(time1 = Dates[Nb.Dates.Obs], time2 = date.0))
  #   r <- r+1
  # }
  latrine <- rep(x = Desman$Latrine_id[k], times = r)
  Df.new <- data.frame(latrine, left, right)
  Df <- rbind.data.frame(Df, Df.new)
}
write.csv(x = Df, file = "DesmanT.csv")

library(frailtypack)
# runShiny()
Desman.Surv <- SurvIC(t0 = 0, lower = Df$left, upper = Df$right, 
                      event = rep(1, length = nrow(Df)))

frailtyPenal(formula = Desman.Surv ~ cluster(latrine), data = Df, RandDist="Gamma",
             n.knots = 6, kappa = 5000)







