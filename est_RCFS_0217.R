require(tensorA)

### transform the matrix to tilde version
mat_trans <- function(Y, alpha)
{
  Y <- to.tensor(Y)
  alpha_tilde <- sqrt(alpha+1) -1;
  Y_avg <- mean(Y,along = 3);
  Y_bar <- to.tensor(array(Y_avg,c(dim(Y)[1],dim(Y)[2],dim(Y)[3])))
  Y_tilde <- Y + alpha_tilde * Y_bar;
  
  return(Y_tilde)
}

MRC_hat <- function(Y)
{
  p = dim(Y)[1]; q = dim(Y)[2];
  MR = matrix(0, p, p); MC = matrix(0, q, q);
  
  for(i in 1:dim(Y)[3])
  {
    MR = MR+Y[,,i]%*%t(Y[,,i]);
    MC = MC+t(Y[,,i])%*%Y[,,i];
  }
  
  MR <- MR/(prod(dim(Y))); MC <- MC/prod(dim(Y))
  
  return(list(MR = MR, MC = MC))
}

### get the k largest eigenvalues and corresponding eigenvector based on the matrix M
k_eigen_dec <- function(M, k = NA)
{
  res <- eigen(M);
  vals <- res$values;
  vecs <- res$vectors;
  
  ### estimate k by eigenvalue ratio test
  eigv <- sort(vals,decreasing = T);
  k_max <- ceiling(length(eigv)/2);
  ratios <- eigv[1:k_max]/eigv[2:(k_max+1)];
  khat = which.max(ratios)
  
  if(is.na(k))
    k <- khat
  
  # idx_khat <- order(vals,decreasing = T)[1:khat];
  # vals_khat <- vals[idx_khat];
  # Vecs_khat <- vecs[,idx_khat];
  # Vals_khat <- diag(vals_khat);
  
  idx_k <- order(vals,decreasing = T)[1:k];
  vals_k <- vals[idx_k];
  Vecs_k <- vecs[,idx_k];
  Vals_k <- diag(vals_k);
  
  return(list(khat = khat, val_k = Vals_k, vec_k = Vecs_k))
}

H_matrix <- function(F_tilde, R, C, V_R, V_C, Rhat, Chat)
{
  t <- dim(F_tilde)[3]; p = dim(R)[1]; q = dim(C)[1];
  
  V_R_inv <- solve(V_R); V_C_inv <- solve(V_C)
  
  FCCF <- matrix(0, dim(R)[2], dim(R)[2])
  FRRF <- matrix(0, dim(C)[2], dim(C)[2])
  
  RRV <- t(R)%*%Rhat%*%V_R_inv;
  CCV <- t(C)%*%Chat%*%V_C_inv
  
  for(i in 1:t)
  {
    FCCF <- FCCF + F_tilde[,,i]%*%t(C)%*%C%*%t(F_tilde[,,i]);
    
    FRRF <- FRRF + t(F_tilde[,,i])%*%t(R)%*%R%*%F_tilde[,,i];
  }
  
  H_R <- FCCF%*%RRV / (p*q*t);
  H_C <- FRRF%*%CCV / (p*q*t)
  
  return(list(H_R = H_R, H_C = H_C))
}

estimation_RCFS <- function(Y, alpha, k = NA, r = NA)
{
  Y_tilde <- mat_trans(Y,alpha);
  
  res <- MRC_hat(Y_tilde);
  MR_hat <- res$MR; MC_hat <- res$MC
  
  kR_hat <- k_eigen_dec(MR_hat, k);
  R_hat <- sqrt(dim(Y)[1]) * kR_hat$vec_k;
  V_R_hat <- kR_hat$val_k;
  khat <- kR_hat$khat;
  if(is.na(k)) 
    k <- khat;
  
  kC_hat <- k_eigen_dec(MC_hat, r);
  C_hat <- sqrt(dim(Y)[2])* kC_hat$vec_k;
  V_C_hat <- kC_hat$val_k;
  rhat <- kC_hat$khat
  if(is.na(r))
    r <- rhat;
  
  F_hat <- to.tensor(array(0,c(dim(R_hat)[2], dim(C_hat)[2], dim(Y)[3])));
  
  S_hat <- to.tensor(array(0, dim(Y)));
  
  for(i in 1:dim(Y)[3])
  {
    F_hat[,,i] <- t(R_hat)%*%Y[,,i]%*%C_hat / (dim(Y)[1]*dim(Y[2]));
    S_hat[,,i] <- R_hat%*%F_hat[,,i]%*%t(C_hat);
  }
  
  return(list(khat = khat, rhat = rhat, 
              Rhat = R_hat, Chat = C_hat, 
              Fhat = F_hat, Shat = S_hat, 
              VRhat = V_R_hat, VChat = V_C_hat))
}
