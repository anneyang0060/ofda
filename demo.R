source('FNS_helper.R')
source('sparse.R')

### basic functions
{  
  fun_mu <- function(t){ 2 * sin(2*pi*t) }
  fun_phi <- function(t){
    phi <- matrix(0,length(t),Mpc)
    phi[,1] <- rep(1,length(t))
    for(i in 2:Mpc){
      phi[,i] <- sqrt(2) * cos((i-1)*pi*t)
    }
    return(phi)
  }
  
  gene_data <- function(mk, njk){
    
    t <- lapply(1:mk, function(i) runif(njk[i],a,b))
    e <- lapply(1:mk, function(i) rnorm(njk[i],0,0.5))
    mu <- lapply(t, fun_mu)
    phi <- lapply(t, fun_phi)
    kesi <- c()
    for(i in 1:Mpc){
      kesi <- cbind(kesi,c(rnorm(mk,0,sqrt(lam[i]))))
    }
    y <- lapply(1:mk, function(i) mu[[i]] + 
                  rowSums(matrix(rep(kesi[i,],each=njk[i]),ncol=Mpc)
                          * phi[[i]]) + e[[i]])
    mylist <- list(t,y)
    names(mylist) <- c('t','y')
    return(mylist)
  }
}

### common parameters
{
  Mpc <- 10
  lam <- 0.4*(1:Mpc)^(-2)
  
  a <- 0; b <- 1
  EV1 <- 100; EV2 <- 50
  eval_mu <- seq(a,b,length.out = EV1)
  mu_true <- fun_mu(eval_mu)
  eval_gam_vec <- seq(a,b,length.out = EV2)
  eval_gam_mat <- cbind(rep(eval_gam_vec,each=EV2), rep(eval_gam_vec,EV2))
  gam_true <- c()
  for(i in 1:EV2)
    for(j in 1:EV2){
      gam_true[(i-1)*EV2+j] <- sum(lam * fun_phi(eval_gam_vec[i])
                                   * fun_phi(eval_gam_vec[j]))
    }
  
  G <- 0.7
  Kmax <- 10; K1 <- 4
  mk_mean <- 20; mk_std <- 3
  mk_ <- ceiling(rnorm(Kmax, mk_mean, mk_std))
  mk_[1] <- 40
  njk_mean <- 6; njk_std <- 2
  njk <- sapply(1:(2*mk_mean*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
  njk <- njk[which(njk<=9 & njk>=3)]
  njk <- njk[1:sum(mk_)]
}

### main process
mfull_ <- 0
for(K in 1:Kmax){
  
  set.seed(K)
  # generate data
  {
    njk1 <- njk[(mfull_+1):(mfull_+mk_[K])]
    data <- gene_data(mk_[K], njk1)
    mfull_ <- mfull_ + mk_[K]
  }
  
  res <- ofda_sparse(data,K,K1=5)
  plot(res$mean,type='l',main=paste0('K=',K))
  lines(mu_true,col='gray')
}
