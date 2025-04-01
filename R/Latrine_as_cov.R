DesmanFeb24 <- read.csv2("./data/survcalibfeb2024_1.csv")

library(survival)

# idx <- which(DesmanFeb24$durationd == 0)
# DesmanFeb24$durationd[idx] <- NA
# 
# Desman.Surv <- with(data = DesmanFeb24, expr = 
#                       Surv(time = durationd, time2 = durationg,
#                            type = "interval2"))

# Avec les censures à droite

library(frailtypack)

idx <- which(is.na(DesmanFeb24$durationg))
event <- rep(1, length = nrow(DesmanFeb24))
event[idx] <- 0
DesmanFeb24$durationg[idx] <- DesmanFeb24$durationd[idx]
DesmanFeb24 <- cbind.data.frame(DesmanFeb24, event)

write.csv(x = DesmanFeb24, file = "data/survcalibfeb2024_2.csv")

Desman.Surv <- SurvIC(t0 = 0, lower = DesmanFeb24$durationd, upper = DesmanFeb24$durationg, 
                      event = event)

frailtyPenal(formula = Desman.Surv ~ cluster(siteID), data = DesmanFeb24, 
             RandDist="Gamma", n.knots = 6, kappa = 5000)
#Avec Df

idx <- which(is.na(Df$durationg))
event <- rep(1, length = nrow(Df))
event[idx] <- 0
Df$durationg[idx] <- Df$durationd[idx]
Df <- cbind.data.frame(Df, event)

write.csv(x = Df, file = "data/survcalibfeb2024_2.csv")

Desman.Surv <- SurvIC(t0 = 0, lower = Df$durationd, upper = Df$durationg, 
                      event = event)

frailtyPenal(formula = Desman.Surv ~ cluster(siteID), data = Df, 
             RandDist="Gamma", n.knots = 6, kappa = 5000)

# Sans les censures à droite

DesmanFeb24.v2 <- DesmanFeb24[-idx, ]            
event.v2 <- event[-idx]

Desman.Surv.v2 <- SurvIC(t0 = 0, lower = DesmanFeb24.v2$durationd, 
                         upper = DesmanFeb24.v2$durationg, 
                         event = event.v2)

frailtyPenal(formula = Desman.Surv.v2 ~ cluster(siteID), data = DesmanFeb24.v2, 
             RandDist="Gamma", n.knots = 6, kappa = 5000)

# Jeu de données dupliqué

DesmanFeb24.v3 <- rbind.data.frame(DesmanFeb24, DesmanFeb24)
event.v3 <- c(event, event)
Desman.Surv.v3 <- SurvIC(t0 = 0, lower = DesmanFeb24.v3$durationd, 
                         upper = DesmanFeb24.v3$durationg, 
                         event = event.v3)
frailtyPenal(formula = Desman.Surv.v3 ~ cluster(siteID), data = DesmanFeb24.v3, 
             RandDist="Gamma", n.knots = 6, kappa = 5000)



parfm(Desman.Surv ~ 1, cluster = "siteID", data = DesmanFeb24, 
      dist = "exponential", frailty = "gamma")


# Weibull sans cov

Desman.Fit1 <- survreg(formula = Desman.Surv ~ 1, dist = "weibull",
                        data = DesmanFeb24)
summary(Desman.Fit1)

# Weibull sans cov siteID

DesmanFeb24$siteID <- as.factor(DesmanFeb24$siteID)

Desman.Fit2 <- survreg(formula = Desman.Surv ~ siteID, dist = "weibull",
                       data = DesmanFeb24)
summary(Desman.Fit2)


