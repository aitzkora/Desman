# Paramters

par <- c(151.34,  0.55, 0.43)
alpha <- par[1] # baseline distribution = weibull : scale
beta <- par[2] # baseline distribution = weibull : shape
lambda <- par[3] # random effect = gamma : shape : lambda, scale : 1/lambda

# Cumulative distribution function of T
Frep <- function(t) {
  g <- function(s) {
    pweibull(q = t/s, scale = alpha, shape = beta)*dgamma(x = s, shape = lambda, scale = 1/lambda)
  }
  aux <- integrate(f = g, lower = 0, upper = Inf)
  res <- aux$value
  return(res)
}

# Simulation - case 1: single time
set.seed(123)
s <- rgamma(n = 1, shape = lambda, scale = 1/lambda)
tps <- rweibull(n = 1, shape = beta, scale = s*alpha)

# Simulation - case 2: multiple times
set.seed(123)
nb.sim <- 100
tps <- vector(mode = "numeric", length = nb.sim)
for (i in 1:nb.sim) {
  s <- rgamma(n = 1, shape = lambda, scale = 1/lambda)
  tps[i] <- rweibull(n = 1, shape = beta, scale = s*alpha)
}

plot(ecdf(tps), do.points = FALSE)
tmax <- max(tps)
t.set <- seq(from = 0, by = 0.1, to = tmax)
Ft <- NULL
for (t in t.set) {
  aux <- Frep(t)
  Ft <- c(Ft, aux)
}
lines(x = t.set, y = Ft, col = "red")
boxplot(tps)

# Simulation - case 3: up to time tmax
set.seed(12)
tmax <- 200

tc <- 0
tps <- NULL
s <- rgamma(n = 1, shape = lambda, scale = 1/lambda)
cst <- pweibull(q = 10, shape = beta, scale = s*alpha, lower.tail = FALSE) 
while (tc < tmax) {
  u <- runif(n = 1)
  aux.tps <- qweibull(p = u*cst, shape = beta, scale = s*alpha)
  tc <- tc + aux.tps
  tps <- c(tps, tc)
} 




