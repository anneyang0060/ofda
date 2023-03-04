ofda_sparse <- function(data, K, K1, L1=10, L2=5, EV1=100, EV2=50, G = 0.7)
{
  a <- 0; b <- 1
  eval_mu <- seq(a,b,length.out = EV1)
  eval_gam_vec <- seq(a,b,length.out = EV2)
  eval_gam_mat <- cbind(rep(eval_gam_vec,each=EV2), rep(eval_gam_vec,EV2))
  
  # initialize
  if(K==1)
  {  
    N <<- 0; mfull <<- 0; N_gam <<- 0
    theta_mu <<- 0; sigma_mu <<- 0
    h_old<<-1
    sigma_gam3 <<- 0; sigma_gam5 <<-0
    
    # lists for bandwidth selection of mean
    res_theta_mu <<- list()
    res_theta_mu$centroids <<- rep(0, L1)
    res_theta_mu$P <<- array(0, dim = c(4,4,EV1,L1))
    res_theta_mu$q <<- array(0, dim = c(4,EV1,L1))
    res_sigma_mu1 <<- list()
    res_sigma_mu1$centroids <<- rep(0, L1)
    res_sigma_mu1$P <<- array(0, dim = c(2,2,EV1,L1))
    res_sigma_mu1$q <<- array(0, dim = c(2,EV1,L1))
    res_sigma_mu2 <<- res_sigma_mu1
    # list of estimating mean
    res_mu <<- res_sigma_mu1
    
    # lists for bandwidth selection of cov
    res_sigma_gam1 <<- list()
    res_sigma_gam1$centroids <<- rep(0, 3)
    res_sigma_gam1$P <<- array(0, dim = c(3,3,EV2^2,3))
    res_sigma_gam1$q <<- array(0, dim = c(3,EV2^2,3))
    res_sigma_gam2 <<- res_sigma_gam1
    res_sigma_gam3 <<- res_sigma_gam1
    res_theta_gam <<- list()
    res_theta_gam$centroids <<- rep(0, 3)
    res_theta_gam$P <<- array(0, dim = c(10,10,EV2^2,3))
    res_theta_gam$q <<- array(0, dim = c(10,EV2^2,3))
    
    # lists for estimating  cov
    res_gam <<- list()
    res_gam$centroids <<- rep(0, L2)
    res_gam$P <<- array(0, dim = c(3,3,EV2^2,L2))
    res_gam$q <<- array(0, dim = c(3,EV2^2,L2))
    
  }
  
  # mean function estimation
  {
    # observe data
    {
      x <- unlist(data$t); y <- unlist(data$y)
      mk <- length(data$t)
      njk1 <- sapply(1:mk, function(i){length(data$t[[i]])})
      NK <<- length(y); N <<- N + NK; mfull <<- mfull + mk
    }
    
    if(K <= K1)
    {
      # theta
      {
        h_theta_mu <<- G * N^(-1/7)
        res_theta_mu <<- online_LCub(x, y, eval_mu, h_theta_mu, L1, 
                                    res_theta_mu, N, NK, 1)
        mu_sec_deri <<- sapply(1:EV1, function(i){
          2*(solve(res_theta_mu$P[,,i,1]+diag(1e-12,4)) %*% matrix(res_theta_mu$q[,i,1],4,1))[3]
        })
        theta_mu <<-  mean(res_theta_mu$P[1,1,6:95 ,1] * mu_sec_deri[6:95]^2)
      }
      
      # sigma
      {
        h_sigma_mu <<- G * N^(-1/5)
        res_sigma_mu1 <<- online_LL(x, y, eval_mu, h_sigma_mu, L1, res_sigma_mu1, N, NK, 1)
        mu <- sapply(1:EV1, function(i){
          (solve(res_sigma_mu1$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu1$q[,i,1],2,1))[1]
        })
        mu_est<<-approx(eval_mu, mu, xout=x, method = 'linear')$y
        r <<- (y-mu_est)^2
        res_sigma_mu2 <<- online_LL(x, r, eval_mu, h_sigma_mu, L1, res_sigma_mu2, N, NK, 1)
        r <<- sapply(1:EV1, function(i){
          (solve(res_sigma_mu2$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu2$q[,i,1],2,1))[1]
        })
        sigma_mu <<- (N-NK)/N * sigma_mu + NK/N * mean(r)
      }
    }
    
    # main regression
    {
      h_mu <<- min((15 * sigma_mu / theta_mu)^(1/5) * N^(-1/5),h_old)
      res_mu <<- online_LL(x, y, eval_mu, h_mu, L1, res_mu, N, NK, 1)
      mu <- sapply(1:EV1, function(i){
        (solve(res_mu$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_mu$q[,i,1],2,1))[1]
      })
      h_old <<- h_mu
    }
  }
  
  # covariance function estimation
  {
    #  observe data
    {
      mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
      res_gam_data <- gene_gam_data(NK, njk1, mk, x, y, mu_est)
      u <- res_gam_data$u; v <- res_gam_data$v
      NK_gam <- sum(njk1 * (njk1 - 1))
      N_gam <<- N_gam + NK_gam
    }
    
    if(K <= K1)
    {
      # theta
      {
        h_theta_gam <- G * N_gam^(-1/8)
        res_theta_gam <<- online_LCub(u, v, eval_gam_mat, h_theta_gam, 3, 
                                     res_theta_gam, N_gam, NK_gam, 2)
        gam_sec_deri <<- sapply(1:EV2^2, function(i){
          2*(sum((solve(res_theta_gam$P[,,i,1]+
                          diag(1e-12,10)) %*% 
                    matrix(res_theta_gam$q[,i,1],10,1))[c(4,6)]))
        })
        gam_sec_deri <<- matrix(gam_sec_deri, EV2, EV2)
        den <<- matrix(res_theta_gam$P[1,1, ,1], EV2, EV2)
        theta_gam <<-mean(den[3:47,3:47] * gam_sec_deri[3:47,3:47]^2)#242
      }
      # sigma
      {
        h_sigma_gam <- G * N_gam^(-1/6)
        res_sigma_gam1 <<- online_LL(u, v, eval_gam_mat, h_sigma_gam,
                                    3, res_sigma_gam1, N_gam, NK_gam, 2)
        gam <<- sapply(1:EV2^2, function(i){
          (solve(res_sigma_gam1$P[,,i,1]+diag(1e-12,3)) %*%
             matrix(res_sigma_gam1$q[,i,1],3,1))[1]
        })
        r <<- (v - sapply(1:NK_gam, function(i){
          gam[which.min(abs(u[i,1]-eval_gam_mat[,1])+
                          abs(u[i,2]-eval_gam_mat[,2]))]
        }))^2
        r[r>(mean(r)+5*sqrt(var(r)))] <<- mean(r)
        r[r<(mean(r)-5*sqrt(var(r)))] <<- mean(r)
        res_sigma_gam2 <<- online_LL(u, r, eval_gam_mat, h_sigma_gam,
                                    3, res_sigma_gam2, N_gam, NK_gam, 2)
        r <<- sapply(1:EV2^2, function(i){
          (solve(res_sigma_gam2$P[,,i,1]+diag(1e-12,3)) %*%
             matrix(res_sigma_gam2$q[,i,1],3,1))[1]
        })
        r <<- matrix(r,EV2,EV2)
        gam <<- matrix(gam, EV2, EV2)
        sigma_eps <<- sigma_mu - mean(diag(gam)) 
        V1 <<- mean(r) + 2*mean(diag(gam))*sigma_eps + sigma_eps^2 - mean(gam^2)
        V1 <<- V1 + mean(diag(r)) + 2*mean(diag(gam))*sigma_eps +
          sigma_eps^2 - mean(diag(gam)^2)
        sigma_gam5 <<- (N_gam-NK_gam)/N_gam * sigma_gam5 + NK_gam/N_gam * V1
      }
      
      # C
      {
        Cgam <<- (2 / 0.04 * 0.6^2 * sigma_gam5 / theta_gam)^(1/6)
      }
    }
    
    # main regression
    {
      h_gam <<- Cgam*N_gam^(-1/6)
      res_gam <<- online_LL(u, v, eval_gam_mat, h_gam, L2, res_gam, N_gam, NK_gam, 2)
      gam <- sapply(1:EV2^2, function(i){
        (solve(res_gam$P[,,i,1]+diag(1e-12,3)) %*% matrix(res_gam$q[,i,1],3,1))[1]
      })
    }
  }
  
  res <- list(mu,gam)
  names(res) <- c('mean','cov')
  return(res)
}