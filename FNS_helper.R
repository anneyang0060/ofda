Epan <- function(z){
  return( 3/4 * (1-z^2) * (abs(z)<1) )
} 

gene_gam_data<-function(NK, njk, mk, x, y, mu_est){
  
  coor1 <- unlist(sapply(1:NK, function(i) 
    rep((1:NK)[i],njk[min(which(i<=cumsum(njk)))]-1)))
  II <- c(0,cumsum(njk))
  coor2 <- unlist(sapply(1:mk,function(i){
    a <- (1:NK)[(II[i]+1):II[i+1]]
    b <- as.vector(sapply(1:length(a), function(j) a[-j]))
  }))
  rm(II)
  u <- t(sapply(1:length(coor1), function(i) c(x[coor1[i]],x[coor2[i]])))
  v <- sapply(1:length(coor1), function(i) 
    (y[coor1[i]]-mu_est[coor1[i]])*(y[coor2[i]]-mu_est[coor2[i]]))
  rm(coor1,coor2)
  
  res <- list(u,v)
  names(res) <- c('u','v')
  return(res)
}

online_LCub <- function(x, y, eval, h, L, res_list, N, n, d){
  
  eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(6+d)) * h}) 
  
  if(K>1){
    idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})
  }else{
    idx <- 1:L
  }
  
  res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
  
  if(d==1){ 
    
    EV <- length(eval)
    
    for(l in 1:L){
      
      Pnew <- array(0, dim = c(4,4,EV)); qnew <- matrix(0,4,EV)
      
      for(i in 1:EV){
        
        side <- cbind(1, x - eval[i], (x - eval[i])^2, (x - eval[i])^3)
        K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
        for(nr in 1:4){
          for(nc in 1:4){
            Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n
          }
          qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
        }
      }
      
      res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
      res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
    }
  }else{
    
    EV <- nrow(eval)
    
    for(l in 1:L){
      Pnew <- array(0, dim = c(10,10,EV)); qnew <- matrix(0,10,EV)
      for(i in 1:EV){
        side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2],
                      (x[,1] - eval[i,1])^2,
                      (x[,1] - eval[i,1])*(x[,2] - eval[i,2]), 
                      (x[,2] - eval[i,2])^2, 
                      (x[,1] - eval[i,1])^3, 
                      (x[,1] - eval[i,1])^2*(x[,2] - eval[i,2]),
                      (x[,1] - eval[i,1])*(x[,2] - eval[i,2])^2,
                      (x[,2] - eval[i,2])^3)# dim: N*10
        K_vec <- (Epan((x[,1] - eval[i,1])/eta[l])/eta[l]
                  *Epan((x[,2] - eval[i,2])/eta[l])/eta[l])
        for(nr in 1:10){
          for(nc in 1:10){
            Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n
          }
          qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
        }
      }
      
      res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
      res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
    }
  }
  
  return(res_list)
}

online_LL <- function(x, y, eval, h, L, res_list, N, n, d){
  
  eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(d+4)) * h})
  
  if(K>1){
    idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})
  }else{
    idx <- 1:L
  }
  
  res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
  
  if(d==1){
    
    EV <- length(eval)
    
    for(l in 1:L){
      
      Pnew <- array(0, dim = c(2,2,EV)); qnew <- matrix(0,2,EV)
      
      for(i in 1:EV){
        side <- cbind(1, x - eval[i])
        K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
        Pnew[,,i] <- matrix(c(
          sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]),
          sum(K_vec*side[,2]*side[,1]),  sum(K_vec*side[,2]^2)
        ),2,2) / n
        qnew[,i] <- matrix(c(
          sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y)
        ),2,1) / n
      }
      
      res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
      res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
      
    }
  }else{
    
    EV <- nrow(eval)
    
    for(l in 1:L){
      
      Pnew <- array(0, dim = c(3,3,EV)); qnew <- matrix(0,3,EV)
      
      for(i in 1:EV){
        side <- cbind(1, x[,1] - eval[i,1], x[,2] - eval[i,2])
        K_vec <- (Epan((x[,1] - eval[i,1])/eta[l])/eta[l]
                  *Epan((x[,2] - eval[i,2])/eta[l])/eta[l])
        Pnew[,,i] <- matrix(c(
          sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]), sum(K_vec*side[,1]*side[,3]),
          sum(K_vec*side[,2]*side[,1]), sum(K_vec*side[,2]^2), sum(K_vec*side[,2]*side[,3]),
          sum(K_vec*side[,3]*side[,1]), sum(K_vec*side[,3]*side[,2]), sum(K_vec*side[,3]^2)
        ),3,3) / n
        qnew[,i] <- matrix(c(
          sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y),
          sum(K_vec*side[,3]*y)
        ),3,1) / n
      }
      
      res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
      res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
    }
  }
  
  return(res_list)
}
