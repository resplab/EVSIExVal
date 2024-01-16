
#########################################################

logit <- function(x) {log(x/(1-x))}
expit <- function(x) {1/(1+exp(-x))}



pdf <- function(x, mu, sigma) {1/sigma/sqrt(2*pi)/x/(1-x)*exp(-(logit(x)-mu)^2/(2*sigma^2))}
cdf <- function(x, mu, sigma) {pnorm(logit(x),mu,sigma)}


solve <- function(target=c(m=0.5,c=0.75), type="logit-normal", init=NULL)
{
  m <- target[1]
  c <- target[2]
  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0#expit(x[1]-6*x[2])
    U <- 1#expit(x[1]+6*x[2])
    f1 <- integrate(cdf, L, U , mu=x[1],sigma=(x[2]))$value
    f2 <- integrate(function(x, mu, sigma) {cdf(x,mu,sigma)^2} , L, U ,mu=x[1], sigma=(x[2]))$value

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

m <- 0.99
c <- 0.90

res <- solve(c(m,c), init=c(5.18,1.10))


x <- (1:1000)/1001
plot(x, pdf(x,res[1],res[2]), type='l')
plot(x, cdf(x,res[1],res[2]), type='l')


###################################################################
simulate <- function(type=c('beta','logit-normal','probit-normal'), args=c(1,1), n_sim=10^6)
{
  if(type=="logit-normal")
  {
    pi <- 1/(1+exp(-rnorm(n_sim,args[1],args[2])))
    Y <- rbinom(n_sim,1,pi)
  }

  require(pROC)

  c(m=mean(pi), c=pROC::roc(Y~pi,)$auc)
}

simulate("logit-normal", res)
