setwd("C:/Users/msafavi/static/GitRepos/EVSIExVal/analysis")

source("include.R")

settings <- list(
  name="case_study",
  output_dir="M:/Projects/2023/Project.EVSIexval/Output/Results/WIP/",
  future_sample_sizes=c(500, 1000, 2000, 4000),
  val_sample_sizes=c(500, 1000, 2000, 4000, 8000),
  n_sim=100, #This one is the outer sim. Inner sim numbers are EVSI() default (10^6)
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
        evidence <- list(prev=c(sum(tmp_data$Y), nrow(tmp_data)-sum(tmp_data$Y)),
                         se=c(sum(tmp_data$Y*(tmp_data$pi>=z)),sum(tmp_data$Y)-sum(tmp_data$Y*(tmp_data$pi>=z))),
                         sp=c(sum((1-tmp_data$Y)*(tmp_data$pi<z)),sum(1-tmp_data$Y)-sum((1-tmp_data$Y)*(tmp_data$pi<z)))
        )
        tmp <- evsiexval::EVSI_ag(evidence, z, settings$future_sample_sizes)
        EVSIs <- rbind(EVSIs, c(val_size=sample_size, z=z, unlist(tmp)))
      }
    }
  }
  EVSIs <- as.data.frame(EVSIs)
  EVSIs
}


#clusterExport(cl,list('EVSI'))

res_cl <- clusterCall(cl, sim_EVSI, model=model, big_val_data=data_us, settings=cl_settings)

stopCluster(cl)


EVSIs <- res_cl[[1]]
for(i in 2:length(res_cl))
{
  EVSIs <- rbind(EVSIs, res_cl[[i]])
}

out$sim_results <- EVSIs


y <- sqldf("SELECT * FROM EVSIs ORDER BY z, val_size")
y$i <-0
cur_z <- 0; cur_val_size <-0
for(i in 1:nrow(y))
{
  if(y$z[i]==cur_z & y$val_size[i]==cur_val_size)
  {
    index <- index+1
  }
  else
  {
    index <- 1
    cur_z <- y$z[i]
    cur_val_size <- y$val_size[i]
  }
  y$i[i] <- index
}
if(max(y$i)>settings$n_sim) y <- y[-which(y$i>settings$n_sim),]


EVSIs <- y


###New graph (2023.12.29)
fss <- settings$future_sample_sizes
x <- sqldf("SELECT COUNT(*) AS N, val_size, z, AVG(EVPI) AS evpi, AVG(EVSI1) AS val1, AVG(EVSI2) AS val2, AVG(EVSI3) AS val3, AVG(EVSI4) AS val4 FROM EVSIs GROUP BY val_size, z")
k <- 1000
par(mar=c(4, 4, 4, 1), xpd=TRUE)

pdf(paste0(settings$output_dir,"EVSI_sim_v2.pdf"), width=7, height=5)
z <- 0.02
y <- x[which(x$z==z),5:8]
max_y <- max(k*y[1,])
plot(c(0,fss), c(0,k*y[1,]), type='l', ylim=c(0,max_y), col='black', xlab="Future sample size", ylab="EVSI (X1000)", lwd=1, xaxt = "n")
axis(1, at=c(0,fss), labels=c(0,fss))
lines(c(0,fss), c(0,k*y[2,]), type='l', ylim=c(0,max(k*y[2,])), col='blue', lwd=1)
lines(c(0,fss), c(0,k*y[3,]), type='l', ylim=c(0,max(k*y[3,])), col='darkgreen', lwd=1)
lines(c(0,fss), c(0,k*y[4,]), type='l', ylim=c(0,max(k*y[4,])), col='orange', lwd=1)
lines(c(0,fss), c(0,k*y[5,]), type='l', ylim=c(0,max(k*y[5,])), col='darkred', lwd=1)
z <- 0.01
y <- x[which(x$z==z),5:8]
lines(c(0,fss), c(0,k*y[1,]), type='l', ylim=c(0,max(k*y[1,])), col='black', lwd=1, lty=2)
lines(c(0,fss), c(0,k*y[2,]), type='l', ylim=c(0,max(k*y[2,])), col='blue', lwd=1, lty=2)
lines(c(0,fss), c(0,k*y[3,]), type='l', ylim=c(0,max(k*y[3,])), col='darkgreen', lwd=1, lty=2)
lines(c(0,fss), c(0,k*y[4,]), type='l', ylim=c(0,max(k*y[4,])), col='orange', lwd=1, lty=2)
lines(c(0,fss), c(0,k*y[5,]), type='l', ylim=c(0,max(k*y[5,])), col='darkred', lwd=1, lty=2)
legend("topleft", inset=c(0,0), horiz=F,  legend=c("Current sample size (n)", "500","1000","2000","4000","8000"), lty=c(0,1,1,1,1,1), col=c('white','black','blue','darkgreen','orange','darkred'), lwd=1, cex=0.6, bty="n")
legend(x=1000, y=max_y*1.04, horiz=T,  legend=c("Threshold (z)", "0.01","0.02"), lty=c(0,2,1), col=c('white','black','black'), lwd=1, cex=0.6, bty="n")
dev.off()





#separate panels
fss <- settings$future_sample_sizes[1:4] #AV doesn't like anything above 4000
x <- sqldf("SELECT COUNT(*) AS N, val_size, z, AVG(EVPI) AS evpi, AVG(EVSI1) AS val1, AVG(EVSI2) AS val2, AVG(EVSI3) AS val3, AVG(EVSI4) AS val4, AVG(EVSI5) AS val5, AVG(EVSI6) AS val6 FROM EVSIs GROUP BY val_size, z")
k <- 1000

par(mfrow=c(2,2))
par(mar=c(4, 4, 4, 1), xpd=TRUE)

dev_n <- 500
y <- x[which(x$val_size==dev_n & x$z %in% c(0.01,0.02)),c(3,6:9)]
max_y <- max(k*y[,-1])
plot(c(0,fss), c(0,k*y[1,-1]), type='l', ylim=c(0,max_y), col='blue', lty=2, xlab="Future sample size", ylab="EVSI (X1000)", lwd=2, xaxt = "n")
axis(1, at=c(0,fss), labels=c(0,fss))
lines(c(0,fss), c(0,k*y[2,-1]), type='l', ylim=c(0,max(k*y[1,])), col='black', lwd=2, lty=1)
legend("topleft", inset=c(0,-0.2), legend=c("threshold (z)\n", "0.01","0.02"), lty=c(1,2), col=c('white','blue','black'), lwd=2, cex=0.6, bty="n", horiz=T)
title(paste("Development sample size:",dev_n))


dev_n <- 1000
y <- x[which(x$val_size==dev_n & x$z %in% c(0.01,0.02)),c(3,6:9)]
max_y <- max(k*y[,-1])
plot(c(0,fss), c(0,k*y[1,-1]), type='l', ylim=c(0,max_y), col='blue', lty=2, xlab="Future sample size", ylab="EVSI (X1000)", lwd=2, xaxt = "n")
axis(1, at=c(0,fss), labels=c(0,fss))
lines(c(0,fss), c(0,k*y[2,-1]), type='l', ylim=c(0,max(k*y[1,])), col='black', lwd=2, lty=1)
legend("topleft", inset=c(0,-0.2), legend=c("threshold (z)\n", "0.01","0.02"), lty=c(1,2), col=c('white','blue','black'), lwd=2, cex=0.6, bty="n", horiz=T)
title(paste("Development sample size:",dev_n))


dev_n <- 2000
y <- x[which(x$val_size==dev_n & x$z %in% c(0.01,0.02)),c(3,6:9)]
max_y <- max(k*y[,-1])
plot(c(0,fss), c(0,k*y[1,-1]), type='l', ylim=c(0,max_y), col='blue', lty=2, xlab="Future sample size", ylab="EVSI (X1000)", lwd=2, xaxt = "n")
axis(1, at=c(0,fss), labels=c(0,fss))
lines(c(0,fss), c(0,k*y[2,-1]), type='l', ylim=c(0,max(k*y[1,])), col='black', lwd=2, lty=1)
legend("topleft", inset=c(0,-0.2), legend=c("threshold (z)\n", "0.01","0.02"), lty=c(1,2), col=c('white','blue','black'), lwd=2, cex=0.6, bty="n", horiz=T)
title(paste("Development sample size:",dev_n))



dev_n <- 4000
y <- x[which(x$val_size==dev_n & x$z %in% c(0.01,0.02)),c(3,6:9)]
max_y <- max(k*y[,-1])
plot(c(0,fss), c(0,k*y[1,-1]), type='l', ylim=c(0,max_y), col='blue', lty=2, xlab="Future sample size", ylab="EVSI (X1000)", lwd=2, xaxt = "n")
axis(1, at=c(0,fss), labels=c(0,fss))
lines(c(0,fss), c(0,k*y[2,-1]), type='l', ylim=c(0,max(k*y[1,])), col='black', lwd=2, lty=1)
legend("topleft", inset=c(0,-0.2), legend=c("threshold (z)\n", "0.01","0.02"), lty=c(1,2), col=c('white','blue','black'), lwd=2, cex=0.6, bty="n", horiz=T)
title(paste("Development sample size:",dev_n))




#Older graph
x <- sqldf("SELECT COUNT(*) AS N, val_size, z, AVG(EVPI) AS evpi, AVG(EVSI1) AS val1, AVG(EVSI2) AS val2, AVG(EVSI3) AS val3, AVG(EVSI4) AS val4, AVG(EVSI5) AS val5, AVG(EVSI6) AS val6 FROM y GROUP BY val_size, z")

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
z <- 0.02
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

