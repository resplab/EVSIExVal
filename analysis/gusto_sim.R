setwd("C:/Users/msafavi/static/GitRepos/EVSIExVal")

source("include.R")

settings <- list(
  name="case_study",
  output_dir="M:/Projects/2023/Project.EVSIexval/Output/Results/WIP/",
  future_sample_sizes=c(500, 1000, 2000, 4000, 8000, 16000),
  val_sample_sizes=c(500, 1000, 2000, 4000, 8000, Inf),
  n_sim=10, #This one is the outer sim. Inner sim numbers are EVSI() default (10^6)
  zs=c(0.01,0.02)
)

out <- list() 
#out <- readRDS(paste0(settings$output_dir,"gusto_sim.RDS"))

library(doParallel)
library(parallel)

total_cores <- detectCores(logical = TRUE)  # returns the number of available hardware threads, and if it is FALSE, returns the number of physical cores
n_cores<- total_cores-1

cl <- makeCluster(n_cores)  

cl_settings <- settings
cl_settings$n_sim <- ceiling(settings$n_sim/n_cores)

registerDoParallel(cl)  


sim_EVSI<- function(model, big_val_data, settings) 
{
  EVSIs <- NULL
  pi <- predict(model, type="response", newdata=big_val_data)
  big_val_data$pi <- pi
  for(i in 1:settings$n_sim)
  {
    for(sample_size in settings$val_sample_sizes)
    {
      if(is.infinite(sample_size))
      {
        sample_size <- nrow(big_val_data)
      }
      tmp_data <- big_val_data[sample(1:nrow(big_val_data), sample_size, replace=T),]
      
      for(z in settings$zs)
      {
        cat(i, sample_size, z)
        tmp <- EVSI(model, tmp_data, z, settings$future_sample_sizes)
        EVSIs <- rbind(EVSIs, c(val_size=sample_size, z=z, unlist(tmp)))
      }
    }
  }
  EVSIs <- as.data.frame(EVSIs)
  EVSIs
}



clusterExport(cl,list('EVSI'))

res_cl <- clusterCall(cl, sim_EVSI, model=model, big_val_data=data_us, settings=cl_settings)

stopCluster(cl)


EVSIs <- res_cl[[1]]
for(i in 2:length(res_cl))
{
  EVSIs <- rbind(EVSIs, res_cl[[i]])
}

out$sim_results <- EVSIs

x <- sqldf("SELECT COUNT(*) AS N, val_size, z, AVG(EVPI) AS evpi, AVG(EVSI1) AS val1, AVG(EVSI2) AS val2, AVG(EVSI3) AS val3, AVG(EVSI4) AS val4, AVG(EVSI5) AS val5, AVG(EVSI6) AS val6 FROM EVSIs GROUP BY val_size, z")

z <- 0.01
y <- x[which(x$z==z),]
pdf(paste0(settings$output_dir,"EVSI_sim_",z,".pdf"), width=7, height=5)     
par(las=2)
plot(c(0,settings$future_sample_sizes), c(0,y[1,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='black', xlab="Future sample size", ylab="EVSI", lwd=2, xaxt = "n")
axis(1, at=settings$future_sample_sizes, labels=settings$future_sample_sizes)
lines(c(0,settings$future_sample_sizes), c(0,y[2,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='blue', lwd=2, lty=5)
lines(c(0,settings$future_sample_sizes), c(0,y[3,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='darkgreen', lwd=2, lty=2)
lines(c(0,settings$future_sample_sizes), c(0,y[4,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='orange', lwd=2, lty=4)
lines(c(0,settings$future_sample_sizes), c(0,y[5,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='darkred', lwd=2, lty=6)
lines(c(0,settings$future_sample_sizes), c(0,y[6,5:10]), type='l', ylim=c(0,max(y[,6:10])), col='red', lwd=2, lty=3)
#legend(-700, 1.03*max(y[,5:10]), legend=c("Development sample size\n", "500","1000","2000","4000","8000", paste0("n=",nrow(data_us))), lty=c(0, 1,5,2,4,6,3), col=c('white','black','blue','darkgreen','orange','darkred','red'), lwd=1, cex=0.6, bty="n")
legend(-700, 1.03*max(y[,5:10]), legend=c("Development sample sizes\n", "500","1000","2000","4000","8000"), lty=c(0,1,5,2,4,6), col=c('white','black','blue','darkgreen','orange','darkred'), lwd=1, cex=0.6, bty="n")
dev.off()


z <- 0.02
y <- x[which(x$z==z),]
pdf(paste0(settings$output_dir,"EVSI_sim_",z,".pdf"), width=7, height=5)     
par(las=2)
plot(c(0,settings$future_sample_sizes), c(0,y[1,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='black', xlab="Future sample size", ylab="EVSI", lwd=2, xaxt = "n")
axis(1, at=settings$future_sample_sizes, labels=settings$future_sample_sizes)
lines(c(0,settings$future_sample_sizes), c(0,y[2,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='blue', lwd=2, lty=5)
lines(c(0,settings$future_sample_sizes), c(0,y[3,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='darkgreen', lwd=2, lty=2)
lines(c(0,settings$future_sample_sizes), c(0,y[4,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='orange', lwd=2, lty=4)
lines(c(0,settings$future_sample_sizes), c(0,y[5,5:10]), type='l', ylim=c(0,max(y[,5:10])), col='darkred', lwd=2, lty=6)
lines(c(0,settings$future_sample_sizes), c(0,y[6,5:10]), type='l', ylim=c(0,max(y[,6:10])), col='red', lwd=2, lty=3)
#legend(-700, 1.03*max(y[,5:10]), legend=c("Development sample size\n", "500","1000","2000","4000","8000", paste0("n=",nrow(data_us))), lty=c(0, 1,5,2,4,6,3), col=c('white','black','blue','darkgreen','orange','darkred','red'), lwd=1, cex=0.6, bty="n")
legend(-700, 1.03*max(y[,5:10]), legend=c("Development sample size\n", "500","1000","2000","4000","8000"), lty=c(0,1,5,2,4,6), col=c('white','black','blue','darkgreen','orange','darkred'), lwd=1, cex=0.6, bty="n")
dev.off()
















#Not in the paper
z <- 0.05
y <- x[which(x$z==z),]
#pdf(paste0(settings$output_dir,"EVSI_sim_",z,".pdf"))     
plot(settings$future_sample_sizes, c(0,y[1,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='black', xlab="Future sample size", ylab="EVSI", lwd=2)
lines(settings$future_sample_sizes, c(0,y[2,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='blue', lwd=2, lty=5)
lines(settings$future_sample_sizes, c(0,y[3,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='darkgreen', lwd=2, lty=2)
lines(settings$future_sample_sizes, c(0,y[4,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='orange', lwd=2, lty=4)
lines(settings$future_sample_sizes, c(0,y[5,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='darkred', lwd=2, lty=3)
legend(-700, 1.03*max(y[,6:11]), legend=c("n=500","n=1000","n=2000","n=4000","n=8000"), lty=c(1,5,2,4,3), col=c('black','blue','darkgreen','orange','darkred'), lwd=1, cex=0.6)
#dev.off()

#Not in the paper
z <- 0.1
y <- x[which(x$z==z),]
#pdf(paste0(settings$output_dir,"EVSI_sim_",z,".pdf"))     
plot(settings$future_sample_sizes, c(0,y[1,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='black', xlab="Future sample size", ylab="EVSI", lwd=2)
lines(settings$future_sample_sizes, c(0,y[2,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='blue', lwd=2, lty=5)
lines(settings$future_sample_sizes, c(0,y[3,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='darkgreen', lwd=2, lty=2)
lines(settings$future_sample_sizes, c(0,y[4,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='orange', lwd=2, lty=4)
lines(settings$future_sample_sizes, c(0,y[5,6:11]), type='l', ylim=c(0,max(y[,6:11])), col='darkred', lwd=2, lty=3)
legend(-700, 1.03*max(y[,6:11]), legend=c("n=500","n=1000","n=2000","n=4000","n=8000"), lty=c(1,5,2,4,3), col=c('black','blue','darkgreen','orange','darkred'), lwd=1, cex=0.6)
#dev.off()


str_file <- paste0(settings$output_dir, "gusto_sim.RDS")
saveRDS(out,str_file)
print(paste("Results saved in ", str_file))

