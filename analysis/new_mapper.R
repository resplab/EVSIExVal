
#########################################################

logit <- function(x) {log(x/(1-x))}
expit <- function(x) {1/(1+exp(-x))}


dlogitnormal <- function(x, mu, sigma) {1/sigma/sqrt(2*pi)/x/(1-x)*exp(-(logit(x)-mu)^2/(2*sigma^2))}
plogitnormal <- function(x, mu, sigma) {pnorm(logit(x),mu,sigma)}

dprobitnormal <- function(x, mu, sigma) {dnorm(qnorm(x),mu,sigma)/dnorm(qnorm(x))}
pprobitnormal <- function(x, mu, sigma) {pnorm(qnorm(x),mu,sigma)}




solve_logitnormal <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- solve_logitnormal(c(1-m,c), init)
    return(c(-tmp[1],tmp[2]))
  }

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0#expit(x[1]-6*x[2])
    U <- 1#expit(x[1]+6*x[2])
    f1 <- integrate(plogitnormal, L, U , mu=x[1],sigma=(x[2]))$value
    f2 <- integrate(function(x, mu, sigma) {plogitnormal(x,mu,sigma)^2} , L, U ,mu=x[1], sigma=(x[2]))$value

    (f1-F1)^2+(f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(logit(m), (0.2))
  }
  res <- optim(init, f, method="L-BFGS-B", lower=c(-10,0), upper=c(10,10)) #, control=list(trace=100))

  if(res$convergence==0)
    c(mu=res$par[1], sigma=(res$par[2]))
  else
    NULL
}



solve_probitnormal <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]
  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0#expit(x[1]-6*x[2])
    U <- 1#expit(x[1]+6*x[2])
    f1 <- integrate(pprobitnormal, L, U , mu=x[1],sigma=(x[2]))$value
    f2 <- integrate(function(x, mu, sigma) {pprobitnormal(x,mu,sigma)^2} , L, U ,mu=x[1], sigma=(x[2]))$value

    (f1-F1)^2+(f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(logit(m), (0.2))
  }
  res <- optim(init, f, method="L-BFGS-B", lower=c(-10,0), upper=c(10,10)) #, control=list(trace=100))

  if(res$convergence==0)
    c(mu=res$par[1], sigma=(res$par[2]))
  else
    NULL
}







solve_probitnormal2 <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0
    U <- 1
    f2 <- integrate(function(x, mu) {pprobitnormal(x, mu, sqrt((mu/qnorm(m))^2-1))^2}, L, U, mu=x)$value

    (f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(sign(m-0.5))
  }
  if(m<0.5)
  {
    L <- -10
    U <- qnorm(m)
  }
  else
  {
    L <- qnorm(m)
    U <- 10
  }

  res <- optim((L+U)/2, f, method="Brent", lower=L, upper=U) #, control=list(trace=100))

  if(res$convergence==0)
    c(mu=res$par[1], sigma=sqrt((res$par[1]/qnorm(m))^2-1))
  else
    NULL
}




solve_probitnormal3 <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0
    U <- 1
    f2 <- integrate(function(x, sigma) {pprobitnormal(x,qnorm(m)*sqrt(1+sigma^2),sigma)^2}, L, U, sigma=x)$value

    (f2-F2)^2
  }

  L=0
  U=10

  if(is.null(init))
  {
    init <- 1
  }

  res <- optim(init, f, method="Brent", lower=L, upper=U) #, control=list(trace=100))

  if(res$convergence==0)
    c(mu=qnorm(m)*sqrt(1+res$par[1]^2), sigma=res$par[1])
  else
    NULL
}








solve_beta <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]
  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0#expit(x[1]-6*x[2])
    U <- 1#expit(x[1]+6*x[2])
    f1 <- integrate(pbeta, L, U , x[1], x[2])$value
    f2 <- integrate(function(x, alpha, beta) {pbeta(x, alpha, beta)^2} , L, U , x[1], x[2])$value

    (f1-F1)^2+(f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(logit(m), (0.2))
  }
  res <- optim(init, f, method="L-BFGS-B", lower=c(0.0001,0.0001), upper=c(10000,10000)) #, control=list(trace=100))

  if(res$convergence==0)
    c(alpha=res$par[1], beta=res$par[2])
  else
    NULL
}

solve_beta2 <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0
    U <- 1
    f2 <- integrate(function(x, alpha) {pbeta(x, alpha, alpha*(1-m)/m)^2}, L, U, alpha=x)$value

    (f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(1)
  }
  res <- optim(init, f, method="Brent", lower=c(0.0001), upper=c(10000)) #, control=list(trace=100))

  if(res$convergence==0)
    c(alpha=res$par[1], beta=res$par[1]*(1-m)/m)
  else
    NULL
}



m <- 0.086
c <- 0.81

res_logitnormal <- solve_logitnormal(c(m,c))
res_probitnormal <- solve_probitnormal3(c(m,c))
res_beta <- solve_beta2(c(m,c))

x <- (0:(10000-1))/10000
plot(x, dlogitnormal(x,res_logitnormal[1],res_logitnormal[2]), type='l')
lines(x, dprobitnormal(x,res_probitnormal[1],res_probitnormal[2]), type='l', col='blue')
lines(x, dbeta(x,res_beta[1],res_beta[2]), type='l', col='red')

plot(x, plogitnormal(x,res_logitnormal[1],res_logitnormal[2]), type='l')
lines(x, pprobitnormal(x,res_probitnormal[1],res_probitnormal[2]), type='l', col='blue')
lines(x, pbeta(x,res_beta[1],res_beta[2]), type='l', col='red')


###################################################################
simulate <- function(type=c('beta','logitnormal','probitnormal'), args=c(1,1), n_sim=10^6)
{
  if(type=="logitnormal")
  {
    pi <- 1/(1+exp(-rnorm(n_sim,args[1],args[2])))
    Y <- rbinom(n_sim,1,pi)
  }

  if(type=="probitnormal")
  {
    pi <- pnorm(rnorm(n_sim,args[1],args[2]))
    Y <- rbinom(n_sim,1,pi)
  }

  if(type=="beta")
  {
    pi <- rbeta(n_sim,args[1],args[2])
    Y <- rbinom(n_sim,1,pi)
  }

  require(pROC)

  c(m=mean(pi), c=pROC::roc(Y~pi,)$auc)
}

simulate("logitnormal", res_logitnormal)
simulate("probitnormal", res_probitnormal)
simulate("beta", res_beta)


####################

out <- data.frame(m=double(), c=double(), type=character(), parm1=double(), parm2=double())

ms <- (1:99)/100
cs <- (55:95)/100

N <- length(ms)*length(cs)*3

out[N,1] <- 0

index <- 1
for(m in ms)
  for(c in cs)
  {
    cat(c(m,c),"|");

    res <- solve_logitnormal(c(m,c)); if(is.null(res)) {res<-c(NA,NA); message("Bad")}
    out[index,] <- list(m,c,"logitnormal",res[1],res[2])
    index <- index+1

    res <- solve_probitnormal3(c(m,c)); if(is.null(res)) {res<-c(NA,NA); message("Bad")}
    out[index,] <- list(m,c,"probitnormal",res[1],res[2])
    index <- index+1

    res <- solve_beta2(c(m,c)); if(is.null(res)) {res<-c(NA,NA); message("Bad")}
    out[index,] <- list(m,c,"beta",res[1],res[2])
    index <- index+1
  }
