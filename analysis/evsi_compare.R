future_sample_sizes <- c(500,1000,2000,4000,8000,16000)

set.seed(1)
z=0.01
evidence <- list(prev = c(43L, 457L), se = c(43L, 0L), sp = c(65, 392))
res_a <- EVSI_ag(evidence, z, future_sample_sizes, 10^6)
x <- data.frame(prev=rbeta(1000,evidence$prev[1]+1,evidence$prev[2]+1),se=rbeta(1000,evidence$se[1]+1,evidence$se[2]+1),sp=rbeta(1000,evidence$sp[1]+1,evidence$sp[2]+1))
#res_g <- EVSI_g(x[1:100,],z,n_sim=1000,future_sample_sizes)
res_gc <- EVSI_gf(x,z,future_sample_sizes,1000)
#res_gc2 <- CEVSI(as.matrix(x[101:200,]),z,future_sample_sizes,10000,F)

max_y <- max(c(res_a$EVSI, res_gc$EVSI))
plot(c(0,future_sample_sizes),c(0,res_a$EVSI), ylim=c(0,max_y))
#lines(c(0,future_sample_sizes),c(0,res_g$EVSI))
lines(c(0,future_sample_sizes),c(0,res_gc$EVSI),col='red')



set.seed(1)
z=0.02
evidence <- list(prev = c(43L, 457L), se = c(41L, 2L), sp = c(147, 310))
res_a <- EVSI_ag(evidence, z, future_sample_sizes, 10^6)
x <- data.frame(prev=rbeta(1000,evidence$prev[1]+1,evidence$prev[2]+1),se=rbeta(1000,evidence$se[1]+1,evidence$se[2]+1),sp=rbeta(1000,evidence$sp[1]+1,evidence$sp[2]+1))
res_g <- EVSI_generic(x[1:100,],z,n_sim=1000,future_sample_sizes)
res_gc <- EVSI_gf(x,z,future_sample_sizes,1000)
#res_gc2 <- CEVSI(as.matrix(x[101:200,]),z,future_sample_sizes,10000,F)

max_y <- max(c(res_a$EVSI, res_g$EVSI, res_gc$EVSI))
plot(c(0,future_sample_sizes),c(0,res_a$EVSI), ylim=c(0,max_y))
lines(c(0,future_sample_sizes),c(0,res_g$EVSI))
lines(c(0,future_sample_sizes),c(0,res_gc$EVSI),col='red')


