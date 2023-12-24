dca <- function(val_data, zs=(0:99)/100, weights=NULL)
{
  n <- dim(val_data)[1]
  
  NB_model <- NB_all <- rep(0, length(zs))
  
  if(is.null(weights)) weights <- rep(1,n)
    
  for(j in 1:length(zs))
  {
    NB_model[j] <- sum(weights*(val_data$pi>zs[j])*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
    NB_all[j] <- sum(weights*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
  }
  
  return(data.frame(z=zs, NB_model=NB_model, NB_all=NB_all))
}



voi_ex_glm <- function(model, val_data, method=c("bootstrap","model_based_ll","model_based_bs","asymptotic"), Bayesian_bootstrap=F, n_sim=1000, zs=(0:99)/100, weights=NULL)
{
  n <- dim(val_data)[1]
  
  if(method=="asymptotic")
  {
    if(is.null(weights)) weights <- rep(1,n)
        
    ENB_perfect <- ENB_current <- rep(0, length(zs))
    for(j in 1:length(zs))
    {
      NB_model <- sum(weights*(val_data$pi>zs[j])*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      NB_all <- sum(weights*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      parms <- NB_BVN(val_data$Y,val_data$pi,zs[j], weights)
      # print(paste0(n,"_",zs[j],"_",parms[5]))
      if(is.na(parms[5])){
        ENB_perfect[j] <- ENB_current[j] <- max(0,NB_model,NB_all)
        } else{
      if(parms[5]>0.999999) parms[5]<- 0.999999
      if(parms[5]< -0.999999) parms[5]<- -0.999999
  
    
      tryCatch(
        {ENB_perfect[j] <- do.call(mu_max_truncated_bvn,as.list(parms))}
      , error=function(cond) {
        return(NULL)
      })
      ENB_current[j] <- max(0,NB_model,NB_all)
        }
    }
    return(data.frame(z=zs, ENB_model=NB_model, ENB_all=NB_all , ENB_perfect=ENB_perfect, ENB_current=ENB_current, EVPIv=ENB_perfect-ENB_current, p_useful=NA))
  }
  
  NB_model <- NB_all <- matrix(0, n_sim, ncol=length(zs))
  
  if(method=="bootstrap")
  {
    for(i in 1:n_sim)
    {
      w_x <- bootstrap(n, Bayesian_bootstrap, weights = weights)
      for(j in 1:length(zs))
      {
        NB_model[i,j] <- sum(w_x*(val_data$pi>zs[j])*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
        NB_all[i,j] <- sum(w_x*(val_data$Y-(1-val_data$Y)*zs[j]/(1-zs[j])))/n
      }
      #if(i%%100 ==0) {plot(zs,NB_model[i,]); lines(zs,NB_all[i,])}
    }
  }
  else
  {
    stop("Method ",method," is not recognized.")
  }
  
  ENB_model <- ENB_all <- ENB_perfect <- ENB_current <- EVPIv <- p_useful <- rep(NA,length(zs))
  for(i in 1:length(zs))
  {
    ENB_model[i] <- mean(NB_model[,i])
    ENB_all[i] <- mean(NB_all[,i])
                                              
    ENB_perfect[i] <- mean(pmax(NB_model[,i],NB_all[,i],0))
    ENB_current[i] <- max(ENB_model[i],ENB_all[i],0)
    EVPIv[i] <- ENB_perfect[i] - ENB_current[i]
    p_useful[i] <- mean((pmax(NB_model[,i],NB_all[,i],0)-NB_model[,i])==0)
  }
  
  data.frame(z=zs, ENB_model=ENB_model, ENB_all=ENB_all, ENB_current=ENB_current, ENB_perfect=ENB_perfect, EVPIv=EVPIv, p_useful=p_useful)
}







NB_BVN <- function(y,pi,z,weights=NULL){
  # set up
  n <- length(y)
  a <- as.numeric(pi>z)
  if(is.null(weights))
  {
    rho <- mean(y)
    TPR <- mean(y*a)/rho
    FPR <- mean(a*(1-y))/(1-rho)
  }
  else
  {
    rho <- weighted.mean(y,w = weights)
    TPR <- weighted.mean(y*a, w = weights)/rho
    FPR <- weighted.mean(a*(1-y), w=weights)/(1-rho) 
  }
    tz <- z/(1-z)
  
  # mean
  mu <- c(rho*TPR-tz*(1-rho)*FPR,rho-tz*(1-rho))
  
  # var
  sig_Y <- rho*(1-rho)/n
  sig_TPR <- rho*TPR*(1-rho*TPR)/n
  sig_FPR <- (1-rho)*FPR*(1-(1-rho)*FPR)/n
  
  sig11 <-  sig_TPR + tz^2 * sig_FPR + 2 * tz * rho * (1-rho) * TPR * FPR / n
  sig22 <- (1+tz)^2 * sig_Y
  sig12 <- sig_Y * (1+tz)* (TPR + tz * FPR)
  cor_coefficient <- sig12/(sqrt(sig11)*sqrt(sig22))
  
  # mean of NB_model, mean of NB_all,
  # var of NB_model, var of NB_all,
  # correlation of NB_model and NB_all
  return(c(mu,sig11,sig22,cor_coefficient))
}


mu_max_truncated_bvn <-
  function(mu1, mu2, sig1sq, sig2sq, rho) {
    sig1 <- sqrt(sig1sq)
    sig2 <- sqrt(sig2sq)
    f1 <-  function(mu1, mu2, sig1, sig2, rho) {
      tmp1 <- sig1 - rho * sig2
      tmp2 <-
        (-sig1 * mu2 + rho * sig2 * mu1) / (sig1 * sig2 * sqrt(1 - rho ^
                                                                 2))
      mu1 * (as.numeric(tmp1 > 0) +  0 * as.numeric(tmp1 == 0) * pnorm(tmp2)) -
        pnorm(tmp2) * (-sig1 * dnorm(-mu1 / sig1) + mu1 * pnorm(-mu1 / sig1))
    }
    
    f2 <- function(mu1, mu2, sig1, sig2, rho) {
      tmp1 <- (sig1 - rho * sig2)
      alpha_num <- (sig1 * mu2 - rho * sig2 * mu1)
      beta_num <- sig1 * sig2 * sqrt(1 - rho ^ 2)
      alpha <- alpha_num / tmp1
      beta <- beta_num / tmp1
      
      if (tmp1 > 0) {
        a <- (tmp1 * mu1 - alpha_num) / beta_num
        b <- tmp1 * sig1 / beta_num
        T1_rho <- -1 / sqrt(1 + b ^ 2)
        
        if (tmp1 < 1e-10) {
          T1_1 <- pnorm(alpha_num / beta_num)
        } else{
          T1_1 <-  pnorm((-a / b) / sqrt(1 + (1 / b) ^ 2))
        }
        
        T1 <-
          mu1 * (T1_1  - as.numeric(pmvnorm(
            lower = c(-Inf,-Inf),
            upper = c(-a / sqrt(1 + b ^ 2),-alpha_num / beta_num),
            mean =
              c(0, 0),
            sigma = matrix(c(1, T1_rho, T1_rho, 1), 2)
          )[[1]]))
        a <- (alpha - mu1) / sig1
        b <- beta / sig1
        T2_t <- sqrt(1 + b ^ 2)
        
        if (tmp1 < 1e-10) {
          T2 <-
            -sig1 / b * dnorm(alpha_num / beta_num) * (1 - pnorm(-mu1 / sig1))
        } else{
          T2 <-
            -sig1 / T2_t * dnorm(a / T2_t) *
            (1 - pnorm(-T2_t * alpha_num / beta_num + a * b / T2_t))
        }
        
        return(T1 + T2)
      }
      else{
        beta <- abs(beta)
        
        a <- (abs(tmp1) * mu1 + alpha_num) / beta_num
        b <- abs(tmp1) * sig1 / beta_num
        T1_rho <- -1 / sqrt(1 + b ^ 2)
        
        if (abs(tmp1) < 1e-10) {
          T1_1 <- pnorm(-alpha_num / beta_num)
        } else{
          T1_1 <-  pnorm((-a / b) / sqrt(1 + (1 / b) ^ 2))
        }
        
        T1 <-
          -mu1 * (T1_1  - as.numeric(pmvnorm(
            lower = c(-Inf,-Inf),
            upper = c(-a / sqrt(1 + b ^ 2), alpha_num / beta_num),
            mean =
              c(0, 0),
            sigma = matrix(c(1, T1_rho, T1_rho, 1), 2)
          )[[1]]))
        
        a <- (alpha - mu1) / sig1
        b <- beta / sig1
        T2_t <- sqrt(1 + b ^ 2)
        
        if (abs(tmp1) < 1e-10) {
          T2 <-
            sig1 / b * dnorm(-alpha_num / beta_num) * (1 - pnorm(-mu1 / sig1))
        } else{
          T2 <-
            sig1 / T2_t * dnorm(a / T2_t) * (1 - pnorm(T2_t * alpha_num / beta_num +
                                                         a * b / T2_t))
        }
        
        return(T1 + T2)
        
      }
    }
    
    f1(mu1, mu2, sig1, sig2, rho) + f1(mu2, mu1, sig2, sig1, rho) -
      f2(mu1, mu2, sig1, sig2, rho) - f2(mu2, mu1, sig2, sig1, rho)
  }


bootstrap <- function (n, Bayesian=F, weights=NULL)
{
  if(Bayesian)
  {
    u <- c(0,sort(runif(n-1)),1)
    u <- (u[-1] - u[-length(u)])*n
    if(!is.null(weights)) u <- u*weights*n/sum(u*weights)
  }
  else
  {
    if(is.null(weights)) weights <-rep(1/n,n)
    u <- rmultinom(1,n,weights)
  }
  as.vector(u)
}



EVSI <- function(model, val_data, z, future_sample_sizes, n_sim=10^6, prior=list(prev=c(1,1), sn=c(1,1), sp=c(1,1)))
{
  #browser()
  n <- dim(val_data)[1]
  data <- list(prev=c(sum(val_data$Y), n - sum(val_data$Y)),
               sn=c(sum(val_data$Y*(val_data$pi>z)), sum(val_data$Y*(1-(val_data$pi>z)))),
               sp=c(sum((1-val_data$Y)*(val_data$pi<=z)), sum((1-val_data$Y)*(1-(val_data$pi<=z))))
  )
  
  posterior <- list(prev=prior$prev+data$prev, sn=prior$sn+data$sn, sp=prior$sp+data$sp)
  
  prev <- posterior$prev[1]/sum(posterior$prev)
  sn <- posterior$sn[1]/sum(posterior$sn)
  sp <- posterior$sp[1]/sum(posterior$sp)
  
  cur_NBs <- c(0, 
               prev*sn-(1-prev)*(1-sp)*z/(1-z),
               prev-(1-prev)*z/(1-z)
  )
  
  pop_data <- data.frame(prev=rbeta(n_sim, posterior$prev[1], posterior$prev[2]),
                         sn=rbeta(n_sim, posterior$sn[1], posterior$sn[2]),
                         sp=rbeta(n_sim, posterior$sp[1], posterior$sp[2]))
  true_NBs <- cbind(0, 
                    pop_data$prev*pop_data$sn-(1-pop_data$prev)*(1-pop_data$sp)*z/(1-z),
                    pop_data$prev-(1-pop_data$prev)*z/(1-z)
  )
  
  EVPI <- mean(apply(true_NBs,1,max)) - max(cur_NBs)
  
  if(length(future_sample_sizes)>0)
  {
    EVSI <- SE <- future_sample_sizes*NA
    
    require(progress)
    pb <- progress_bar$new(total=length(future_sample_sizes))
    
    for(i in 1:length(future_sample_sizes))
    {
      future_sample_size <- future_sample_sizes[i]
      future_D <- rbinom(n_sim, size=future_sample_size, prob=pop_data$prev)
      future_TP <- rbinom(n_sim, size=future_D, prob=pop_data$sn)
      future_TN <- rbinom(n_sim, size=future_sample_size-future_D, prob=pop_data$sp)
      
      pooled_data <- list(
        prev=cbind(future_D+data$prev[1],future_sample_size-future_D+data$prev[2]),
        sn=cbind(future_TP+data$sn[1],future_D-future_TP+data$sn[2]),
        sp=cbind(future_TN+data$sp[1],future_sample_size-future_D-future_TN+data$sp[2])
      )
      
      posterior2 <- list(prev=prior$prev+pooled_data$prev, sn=prior$sn+pooled_data$sn, sp=prior$sp+pooled_data$sp)
      
      prev2 <- posterior2$prev[,1]/rowSums(posterior2$prev)
      sn2 <- posterior2$sn[,1]/rowSums(posterior2$sn)
      sp2 <- posterior2$sp[,1]/rowSums(posterior2$sp)
      
      future_NBs <- cbind(0, 
                          prev2*sn2-(1-prev2)*(1-sp2)*z/(1-z),
                          prev2-(1-prev2)*z/(1-z)
      )
      
      winners <- apply(future_NBs, 1, which.max)
      
      #winning_NBs <- future_NBs[cbind(1:n_sim,winners)] /
      winning_NBs <- true_NBs[cbind(1:n_sim,winners)] 
      
      SE[i] <- sqrt(var(winning_NBs)/length(winning_NBs))
      
      EVSI[i] <-  mean(winning_NBs) - max(cur_NBs)
      
      pb$tick()
    }
  }
  else EVSI=NULL
  
  list(EVPI=EVPI, EVSI=EVSI)
}



EVSI_ag <- function(evidence=list(prev=c(1,1), sn=c(1,1), sp=c(1,1)), z, future_sample_sizes, n_sim=10^6, prior=list(prev=c(1,1), sn=c(1,1), sp=c(1,1)))
{
  n <-0
  
  data <- evidence
  
  posterior <- list(prev=prior$prev+data$prev, sn=prior$sn+data$sn, sp=prior$sp+data$sp)
  
  prev <- posterior$prev[1]/sum(posterior$prev)
  sn <- posterior$sn[1]/sum(posterior$sn)
  sp <- posterior$sp[1]/sum(posterior$sp)
  
  cur_NBs <- c(0,
               prev*sn-(1-prev)*(1-sp)*z/(1-z),
               prev-(1-prev)*z/(1-z)
  )
  
  pop_data <- data.frame(prev=rbeta(n_sim, posterior$prev[1], posterior$prev[2]),
                         sn=rbeta(n_sim, posterior$sn[1], posterior$sn[2]),
                         sp=rbeta(n_sim, posterior$sp[1], posterior$sp[2]))
  true_NBs <- cbind(0,
                    pop_data$prev*pop_data$sn-(1-pop_data$prev)*(1-pop_data$sp)*z/(1-z),
                    pop_data$prev-(1-pop_data$prev)*z/(1-z)
  )
  p_best <- as.data.frame(table(apply(true_NBs,1,which.max))/nrow(true_NBs))
  colnames(p_best)  <- c("Decision", "p_best")
  levels(p_best$Decision) <- c("1"="Treat none", "2"="Use Model", "3"="treat All") 
  
  EVPI <- mean(apply(true_NBs,1,max)) - max(cur_NBs)
  
  if(length(future_sample_sizes)>0)
  {
    EVSI <- SE <- future_sample_sizes*NA
    
    require(progress)
    pb <- progress_bar$new(total=length(future_sample_sizes))
    
    for(i in 1:length(future_sample_sizes))
    {
      future_sample_size <- future_sample_sizes[i]
      future_D <- rbinom(n_sim, size=future_sample_size, prob=pop_data$prev)
      future_TP <- rbinom(n_sim, size=future_D, prob=pop_data$sn)
      future_TN <- rbinom(n_sim, size=future_sample_size-future_D, prob=pop_data$sp)
      
      pooled_data <- list(
        prev=cbind(future_D+data$prev[1],future_sample_size-future_D+data$prev[2]),
        sn=cbind(future_TP+data$sn[1],future_D-future_TP+data$sn[2]),
        sp=cbind(future_TN+data$sp[1],future_sample_size-future_D-future_TN+data$sp[2])
      )
      
      posterior2 <- list(prev=prior$prev+pooled_data$prev, sn=prior$sn+pooled_data$sn, sp=prior$sp+pooled_data$sp)
      
      prev2 <- posterior2$prev[,1]/rowSums(posterior2$prev)
      sn2 <- posterior2$sn[,1]/rowSums(posterior2$sn)
      sp2 <- posterior2$sp[,1]/rowSums(posterior2$sp)
      
      future_NBs <- cbind(0,
                          prev2*sn2-(1-prev2)*(1-sp2)*z/(1-z),
                          prev2-(1-prev2)*z/(1-z)
      )
      
      winners <- apply(future_NBs, 1, which.max)
      
      winning_NBs <- true_NBs[cbind(1:n_sim,winners)]
      
      SE[i] <- sqrt(var(winning_NBs)/length(winning_NBs))
      
      EVSI[i] <-  mean(winning_NBs) - max(cur_NBs)
      
      pb$tick()
    }
  }
  else EVSI=NULL
  
  list(EVPI=EVPI, EVSI=EVSI, p_best=p_best)
}



EVSI_generic <- function(samples=data.frame(prev=rbeta(1000,1,1), sn=rbeta(1000,1,1), sp=rbeta(1000,1,1)), z, future_sample_sizes, n_sim=10^3)
{
  M <- nrow(samples)
  
  prev <- mean(samples$prev)
  sn <- mean(samples$sn)
  sp <- mean(samples$sp)
  
  cur_NBs <- c(0,
               prev*sn-(1-prev)*(1-sp)*z/(1-z),
               prev-(1-prev)*z/(1-z)
  )
  
  true_NBs <- cbind(0,
                    samples$prev*samples$sn-(1-samples$prev)*(1-samples$sp)*z/(1-z),
                    samples$prev-(1-samples$prev)*z/(1-z)
  )
  
  bests <- apply(true_NBs,1,which.max)
  
  p_best <- as.data.frame(table(bests)/nrow(true_NBs))
  colnames(p_best)  <- c("Decision", "p_best")
  levels(p_best$Decision) <- c("1"="Treat none", "2"="Use Model", "3"="treat All") 
  
  EVPI <- mean(apply(true_NBs,1,max)) - max(cur_NBs)
  
  if(length(future_sample_sizes)>0)
  {
    EVSI <- SE <- future_sample_sizes*NA
    
    require(progress)
    pb <- progress_bar$new(total=length(future_sample_sizes))
    
    S <- do.call(rbind, replicate(n_sim, samples, simplify=FALSE))
    TNB <- do.call(rbind, replicate(n_sim, true_NBs, simplify=FALSE))
    
    for(i in 1:length(future_sample_sizes))
    {
      future_sample_size <- future_sample_sizes[i]
      
      future_D <- rbinom(nrow(S), size=future_sample_size, prob=S$prev)
      future_TP <- rbinom(nrow(S), size=future_D, prob=S$sn)
      future_TN <- rbinom(nrow(S), size=future_sample_size-future_D, prob=S$sp)
      
      ws <<- rep(rbinom(1,1,0.5), nrow(samples))
      updated_NBs <<- matrix(double(1), nrow=nrow(samples),ncol=3)
      
      find_winner <- function(D_TP_TN) #Takes one realization of future study and updates the evidence and eclares the winner
      {
        ws <<- dbinom(D_TP_TN[1], future_sample_size, samples[,1])*
          dbinom(D_TP_TN[2], D_TP_TN[1], samples[,2])*
          dbinom(D_TP_TN[3], future_sample_size-D_TP_TN[1], samples[,3])
        
        updated_NBs <<- cbind(0, 
                              ws*(samples$prev*samples$sn-(1-samples$prev)*(1-samples$sp)*z/(1-z))/sum(ws),
                              ws*(samples$prev-(1-samples$prev)*z/(1-z))/sum(ws)
        )
        which.max(colMeans(updated_NBs))
      }
      
      winners <- apply(cbind(future_D, future_TP, future_TN), 1, find_winner)
      
      winning_NBs <- TNB[cbind(1:nrow(TNB),winners)]
      
      SE[i] <- sqrt(var(winning_NBs)/length(winning_NBs))
      
      EVSI[i] <-  mean(winning_NBs) - max(cur_NBs)
      
      pb$tick()
    }
  }
  else EVSI=NULL
  
  list(EVPI=EVPI, EVSI=EVSI, p_best=p_best)
}




