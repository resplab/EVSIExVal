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


##' @export
EVSI <- function(model, val_data, z, future_sample_sizes, n_sim=10^6, prior=list(prev=c(1,1), se=c(1,1), sp=c(1,1)))
{
  n <- dim(val_data)[1]
  data <- list(prev=c(sum(val_data$Y), n-sum(val_data$Y)),
               se=c(sum(val_data$Y*(val_data$pi>z)), sum(val_data$Y*(1-(val_data$pi>z)))),
               sp=c(sum((1-val_data$Y)*(val_data$pi<=z)), sum((1-val_data$Y)*(1-(val_data$pi<=z))))
  )

  posterior <- list(prev=prior$prev+data$prev, se=prior$se+data$se, sp=prior$sp+data$sp)

  prev <- posterior$prev[1]/sum(posterior$prev)
  se <- posterior$se[1]/sum(posterior$se)
  sp <- posterior$sp[1]/sum(posterior$sp)

  cur_NBs <- c(0,
               prev*se-(1-prev)*(1-sp)*z/(1-z),
               prev-(1-prev)*z/(1-z)
  )

  pop_data <- data.frame(prev=rbeta(n_sim, posterior$prev[1], posterior$prev[2]),
                         se=rbeta(n_sim, posterior$se[1], posterior$se[2]),
                         sp=rbeta(n_sim, posterior$sp[1], posterior$sp[2]))
  true_NBs <- cbind(0,
                    pop_data$prev*pop_data$se-(1-pop_data$prev)*(1-pop_data$sp)*z/(1-z),
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
      future_TP <- rbinom(n_sim, size=future_D, prob=pop_data$se)
      future_TN <- rbinom(n_sim, size=future_sample_size-future_D, prob=pop_data$sp)

      pooled_data <- list(
        prev=cbind(future_D+data$prev[1],future_sample_size-future_D+data$prev[2]),
        se=cbind(future_TP+data$se[1],future_D-future_TP+data$se[2]),
        sp=cbind(future_TN+data$sp[1],future_sample_size-future_D-future_TN+data$sp[2])
      )

      posterior2 <- list(prev=prior$prev+pooled_data$prev, se=prior$se+pooled_data$se, sp=prior$sp+pooled_data$sp)

      prev2 <- posterior2$prev[,1]/rowSums(posterior2$prev)
      se2 <- posterior2$se[,1]/rowSums(posterior2$se)
      sp2 <- posterior2$sp[,1]/rowSums(posterior2$sp)

      future_NBs <- cbind(0,
                          prev2*se2-(1-prev2)*(1-sp2)*z/(1-z),
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



##' @export
EVSI_ag <- function(evidence=list(prev=c(1,1), se=c(1,1), sp=c(1,1)), z, future_sample_sizes, n_sim=10^6, prior=list(prev=c(1,1), se=c(1,1), sp=c(1,1)), do_EVSIp=TRUE, do_summary=TRUE, progress=TRUE)
{
  out <- list()

  n <-0

  data <- evidence

  posterior <- list(prev=prior$prev+data$prev, se=prior$se+data$se, sp=prior$sp+data$sp)

  prev <- posterior$prev[1]/sum(posterior$prev)
  se <- posterior$se[1]/sum(posterior$se)
  sp <- posterior$sp[1]/sum(posterior$sp)

  pop_data <- data.frame(prev=rbeta(n_sim, posterior$prev[1], posterior$prev[2]),
                         se=rbeta(n_sim, posterior$se[1], posterior$se[2]),
                         sp=rbeta(n_sim, posterior$sp[1], posterior$sp[2]))
  true_NBs <- cbind(0,
                    pop_data$prev*pop_data$se-(1-pop_data$prev)*(1-pop_data$sp)*z/(1-z),
                    pop_data$prev-(1-pop_data$prev)*z/(1-z)
  )

  cur_NBs <- cbind("Treat none"=0,
                   "Use model"=mean(true_NBs[,2]),
                   "Treat all"=mean(true_NBs[,3])
  )

  #For output
  if(do_summary)
  {
    summary <- cur_NBs
    summary <- rbind(summary, apply(true_NBs, 2, quantile, c(0.025,0.975)))
    summary <- rbind(summary, c(0,0,0))
    p_best <- table(apply(true_NBs,1,which.max))/nrow(true_NBs)
    summary[4, as.numeric(colnames(t(as.matrix(p_best))))] <- p_best
    rownames(summary) <- c("Expected NB", "95% CI L", "95% CI H", "P_best")

    out$summary <- summary
  }

  EVPI <- mean(apply(true_NBs,1,max)) - max(cur_NBs)

  if(length(future_sample_sizes)>0)
  {
    EVSI <- SE <- EVSIp <- future_sample_sizes*NA

    if(progress)
    {
      require(progress)
      pb <- progress_bar$new(total=length(future_sample_sizes))
    }

    for(i in 1:length(future_sample_sizes))
    {
      future_sample_size <- future_sample_sizes[i]
      future_D <- rbinom(n_sim, size=future_sample_size, prob=pop_data$prev)
      future_TP <- rbinom(n_sim, size=future_D, prob=pop_data$se)
      future_TN <- rbinom(n_sim, size=future_sample_size-future_D, prob=pop_data$sp)

      pooled_data <- list(
        prev=cbind(future_D+data$prev[1],future_sample_size-future_D+data$prev[2]),
        se=cbind(future_TP+data$se[1],future_D-future_TP+data$se[2]),
        sp=cbind(future_TN+data$sp[1],future_sample_size-future_D-future_TN+data$sp[2])
      )

      posterior2 <- list(prev=prior$prev+pooled_data$prev, se=prior$se+pooled_data$se, sp=prior$sp+pooled_data$sp)

      prev2 <- posterior2$prev[,1]/rowSums(posterior2$prev)
      se2 <- posterior2$se[,1]/rowSums(posterior2$se)
      sp2 <- posterior2$sp[,1]/rowSums(posterior2$sp)

      future_NBs <- cbind(0,
                          prev2*se2-(1-prev2)*(1-sp2)*z/(1-z),
                          prev2-(1-prev2)*z/(1-z)
      )

      winners <- apply(future_NBs, 1, which.max)

      winning_NBs <- true_NBs[cbind(1:n_sim,winners)]

      SE[i] <- sqrt(var(winning_NBs)/length(winning_NBs))

      EVSI[i] <-  mean(winning_NBs) - max(cur_NBs)

      if(do_EVSIp)
      {
        true_winners <- apply(true_NBs, 1, which.max)
        EVSIp[i] <- sum(winners==true_winners)/n_sim
      }

      pb$tick()
    }
  }
  else EVSI=NULL

  out$EVPI <- EVPI
  out$EVSI <- EVSI

  if(do_EVSIp)
  {
    out$EVSIp <- EVSIp
  }

  out
}



##' @export
EVSI_g <- function(samples=data.frame(prev=rbeta(100,1,1), se=rbeta(100,1,1), sp=rbeta(100,1,1)), z, future_sample_sizes, n_sim=10^3)
{
  M <- nrow(samples)

  prev <- mean(samples$prev)
  se <- mean(samples$se)
  sp <- mean(samples$sp)

  #Turn samples to matrix for faster computerions, make sure order is (prev, se, sp)
  samples <- as.matrix(samples[,c('prev','se','sp')])


  true_NBs <- cbind(0,
                    samples[,1]*samples[,2]-(1-samples[,1])*(1-samples[,3])*z/(1-z),
                    samples[,1]-(1-samples[,1])*z/(1-z)
  )

  cur_NBs <- cbind("Treat none"=0,
               "Use model"=mean(true_NBs[,2]),
               "Treat all"=mean(true_NBs[,3])
  )

  #For output
  summary <- cur_NBs
  summary <- rbind(summary, apply(true_NBs, 2, quantile, c(0.025,0.975)))
  summary <- rbind(summary, c(0,0,0))
  p_best <- table(apply(true_NBs,1,which.max))/nrow(true_NBs)
  summary[4, as.numeric(colnames(t(as.matrix(p_best))))] <- p_best
  rownames(summary) <- c("Expected NB", "95% CI L", "95% CI H", "P_best")

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

      future_D <- rbinom(nrow(S), size=future_sample_size, prob=S[,1])
      future_TP <- rbinom(nrow(S), size=future_D, prob=S[,2])
      future_TN <- rbinom(nrow(S), size=future_sample_size-future_D, prob=S[,3])

      ws <<- rep(rbinom(1,1,0.5), nrow(samples))
      updated_NBs <<- matrix(double(1), nrow=nrow(samples),ncol=3)

      find_winner <- function(D_TP_TN) #Takes one realization of future study and updates the evidence and eclares the winner
      {
        ws <<- dbinom(D_TP_TN[1], future_sample_size, samples[,1])*
          dbinom(D_TP_TN[2], D_TP_TN[1], samples[,2])*
          dbinom(D_TP_TN[3], future_sample_size-D_TP_TN[1], samples[,3])

        updated_NBs <<- cbind(0,
                              ws*(samples[,1]*samples[,2]-(1-samples[,1])*(1-samples[,3])*z/(1-z))/sum(ws),
                              ws*(samples[,1]-(1-samples[,1])*z/(1-z))/sum(ws)
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

  list(summary=summary, EVPI=EVPI, EVSI=EVSI)
}



##' @export
EVSI_gf <- function(samples=data.frame(prev=rbeta(100,1,1), se=rbeta(100,1,1), sp=rbeta(100,1,1)), z, future_sample_sizes, n_sim=10^3)
{
  CEVSI(as.matrix(samples[,1:3]), z, future_sample_sizes, n_sim , F)
}


