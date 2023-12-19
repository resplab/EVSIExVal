z <- 0.02

phi <- 0.086

sn <- 0.9534

sp <- 0.678337

sNB <- (phi*sn-(1-phi)*(1-sp)*z/(1-z))/phi

w <- (1-phi)/phi*z/(1-z)

SE_NB <- 0.051

N <- 1/SE_NB^2*(sn*(1-sn)/phi + w^2*sp*(1-sp)/(1-phi) + w^2*(1-sp)^2/phi/(1-phi))

SE_ln_OE <- 0.245

(1-phi)/phi/SE_ln_OE^2