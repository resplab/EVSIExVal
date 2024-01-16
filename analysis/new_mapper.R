a <- 1

delta <- 1

m <- 0.5

x <- (1:1000)/1001

plot(pbeta(x, a, a*(1-m)/m+1), type='l')
lines(pbeta(x, (a+delta), (a+delta)*(1-m)/m+1), type='l', col='blue')

lines(dbeta(x, (a+delta)+1, (a+delta)*(1-m)/m), type='l', col='red')



#################

m <- 0.6

a1 <- 2
b1 <- a1*(1-m)/m

a2 <- 3
b2 <- a2*(1-m)/m

x <- (1:1000)/1001

Y <- 2*x*dbeta(x,a1,b1)*pbeta(x,a1,b1)
Z <- 2*x*dbeta(x,a2,b2)*pbeta(x,a2,b2)

mean(pbeta(x,a1,b1))
mean(pbeta(x,a2,b2))

plot(pbeta(x,a1,b1)^1,type='l')
lines(pbeta(x,a2,b2)^1,col='red')

sum(Y)/length(x)
sum(Z)/length(x)

plot(Y, type='l')
lines(Z, col='blue')


###############2024.01.14

logit <- function(x) {log(x/(1-x))}
expit <- function(x) {1/(1+exp(-x))}


pdf <- function(x, mu, sigma) {1/sigma/sqrt(2*pi)/x/(1-x)*exp(-(logit(x)-mu)^2/(2*sigma^2))}
cdf <- function(x, mu, sigma) {pnorm(logit(x),mu,sigma)}

mu <- -0.4515322
sigma <- 0.7134626

x <- (1:1000)/1001
plot(x, pdf(x,mu,sigma), type='l')
plot(x, cdf(x,mu,sigma), type='l')

n <- 10^7
x<-rnorm(n,mu,sigma)
mean(1/(1+exp(-(x))))
rm(x)

K <- 10000
i <- 1:(K-1)
sum(expit(qnorm(i/K,mu,sigma)))/(K-1)

1-integrate(cdf,0,1,mu,sigma)$value

integrate(function(x,mu,sigma){x*pdf(x,mu,sigma)},0,1,mu,sigma)$value

#########################################################

m <- 0.4
c <- 0.7

F1 <- 1-m

D <- sqrt((1-2*c)*m*m - (1-2*c)*m + 0.25)

F2 <- 1-c(m-0.5+D, m-0.5-D)[1]

pdf <- function(x, mu, sigma) {1/sigma/sqrt(2*pi)/x/(1-x)*exp(-(logit(x)-mu)^2/(2*sigma^2))}
cdf <- function(x, mu, sigma) {pnorm(logit(x),mu,sigma)}
cdf2 <- function(x, mu, sigma) {pnorm(logit(x),mu,sigma)^2}


f <- function(x)
{
  f1 <- integrate(cdf,0,1,mu=x[1],sigma=x[2])$value
  f2 <- integrate(cdf2,0,1,mu=x[1],sigma=x[2])$value

  (f1-F1)^2+(f2-F2)^2
}


optim(c(0,1), f)
