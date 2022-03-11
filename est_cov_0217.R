require(tensorA)

Omega_v <- function(Rhat, Chat, Ehat, Fhat, v, i, j)
{
  p = dim(Rhat)[1]; q = dim(Chat)[1]; k = dim(Rhat)[2]; r = dim(Chat)[2]
  T = dim(Ehat)[3]
  
  Omega_v_R11 <- matrix(0,k,k);
  Omega_v_R12 <- matrix(0,k,r);
  Omega_v_R21 <- matrix(0,r,k);
  Omega_v_R22 <- matrix(0,r,r);
  
  Omega_v_C11 <- matrix(0,r,r);
  Omega_v_C12 <- matrix(0,r,r);
  Omega_v_C21 <- matrix(0,k,r);
  Omega_v_C22 <- matrix(0,k,k);
  
  for(t in (1+v):T)
  {
    F_t <- Fhat[,,t]; F_tv <- Fhat[,,t-v]
    E_ti <- Ehat[i,,t]; E_tv_i <- Ehat[i,,t-v]
    E_tj <- Ehat[,j,t]; E_tv_j <- Ehat[,j,t-v]
    
    temp_R <- t(C_hat) %*% E_t_i %*% t(E_tv_i) %*% C_hat
    temp_C <- t(R_hat) %*% E_t_j %*% t(E_tv_j) %*% R_hat
    
    Omega_v_R11 = Omega_v_R11 + F_t %*% temp_R %*% t(F_tv)
    Omega_v_C11 = Omega_v_C11 + t(F_t) %*% temp_C %*% F_tv
    
    Omega_v_R12 = Omega_v_R12 + F_t %*% temp_R
    Omega_v_C12 = Omega_v_C12 + t(F_t) %*% temp_C
    
    Omega_v_R21 = Omega_v_R21 + temp_R %*% t(F_t)
    Omega_v_C21 = Omega_v_C21 + temp_C %*% F_t
    
    Omega_v_R22 = Omega_v_R22 + temp_R
    Omega_v_C22 = Omega_v_C22 + temp_C
  }
  
  return(list(R11 = Omega_v_R11/(q*T), R12 = Omega_v_R12/(q*T), R21 = Omega_v_R21/(q*T), R22 = Omega_v_R22/(q*T), C11 = Omega_v_C11/(p*T), C12 = Omega_v_C12/(p*T), C21 = Omega_v_C21/(p*T), C22 = Omega_v_C22/(p*T)))
}

cov_loading_HAC <- function(Y, R_hat, C_hat, F_hat, VR, VC, i, j, alpha)
{
  p = dim(R_hat)[1]; q = dim(C_hat)[1]; k = dim(F_hat)[1]; r = dim(F_hat)[2]; T = dim(F_hat)[3]
  
  E_hat <- to.tensor(array(0,dim(Y)));
  for(t in 1:T)
  {
    E_hat[,,t] <- Y[,,t] - to.tensor(R_hat %*% F_hat[,,t] %*% t(C_hat))
  }
  
  F_mean <- mean(F_hat, along = 3)
  F_bar <- to.tensor(array(F_mean, dim(F_hat)))
  F_hat0 <- F_hat
  F_hat <- F_hat0 - F_bar
  
  ### when v = 0
  res_0 <- Omega_v(R_hat, C_hat, E_hat, F_hat, 0, i, j, alpha)
  Omega_v_R11 <- res_0$R11;
  Omega_v_R12 <- res_0$R12;
  Omega_v_R21 <- res_0$R21;
  Omega_v_R22 <- res_0$R22;
  Omega_v_C11 <- res_0$C11;
  Omega_v_C12 <- res_0$C12;
  Omega_v_C21 <- res_0$C21;
  Omega_v_C22 <- res_0$C22;
  
  ### when v = 1,2,....
  mr = round(log(q*T)); mc = round(log(p*T))
  mm = min(mr,mc)
  
  for(v in 1:mm)
  {
    temp_res <- Omega_v(R_hat, C_hat, E_hat, F_hat, v, i, j, alpha)
    Omega_v_R11 <- Omega_v_R11 + (temp_res$R11 + t(temp_res$R11)) * (1 - v/(1+mr))
    Omega_v_R12 <- Omega_v_R12 + (temp_res$R12 + t(temp_res$R12)) * (1 - v/(1+mr))
    Omega_v_R21 <- Omega_v_R21 + (temp_res$R21 + t(temp_res$R21)) * (1 - v/(1+mr))
    Omega_v_R22 <- Omega_v_R22 + (temp_res$R22 + t(temp_res$R22)) * (1 - v/(1+mr))
    Omega_v_C11 <- Omega_v_C11 + (temp_res$C11 + t(temp_res$C11)) * (1 - v/(1+mc))
    Omega_v_C12 <- Omega_v_C12 + (temp_res$C12 + t(temp_res$C12)) * (1 - v/(1+mc))
    Omega_v_C21 <- Omega_v_C21 + (temp_res$C21 + t(temp_res$C21)) * (1 - v/(1+mc))
    Omega_v_C22 <- Omega_v_C22 + (temp_res$C22 + t(temp_res$C22)) * (1 - v/(1+mc))
  }
  
  if(mr > mc)
  {
    for(v in (mm+1):mr)
    {
      temp_res <- Omega_v(R_hat, C_hat, E_hat, F_hat, v, i, j, alpha)
      Omega_v_R11 <- Omega_v_R11 + (temp_res$R11 + t(temp_res$R11)) * (1 - v/(1+mr))
      Omega_v_R12 <- Omega_v_R12 + (temp_res$R12 + t(temp_res$R12)) * (1 - v/(1+mr))
      Omega_v_R21 <- Omega_v_R21 + (temp_res$R21 + t(temp_res$R21)) * (1 - v/(1+mr))
      Omega_v_R22 <- Omega_v_R22 + (temp_res$R22 + t(temp_res$R22)) * (1 - v/(1+mr))
    }
  }
  
  else if(mr < mc)
  {
    for(v in (mm+1):mc)
    {
      temp_res <- Omega_v(R_hat, C_hat, E_hat, F_hat, v, i, j, alpha)
      Omega_v_C11 <- Omega_v_C11 + (temp_res$C11 + t(temp_res$C11)) * (1 - v/(1+mc))
      Omega_v_C12 <- Omega_v_C12 + (temp_res$C12 + t(temp_res$C12)) * (1 - v/(1+mc))
      Omega_v_C21 <- Omega_v_C21 + (temp_res$C21 + t(temp_res$C21)) * (1 - v/(1+mc))
      Omega_v_C22 <- Omega_v_C22 + (temp_res$C22 + t(temp_res$C22)) * (1 - v/(1+mc))
    }
  }
  
  ### choice of best beta
  alpha_1 <- -0.5 * (sum(diag(F_mean %*% Omega_v_R22 %*% t(F_mean))))^(-1) * sum(diag(Omega_v_R12 %*% t(F_mean) + F_mean %*% Omega_v_R21));
  alpha_2 <- -0.5 * (sum(diag(t(F_mean) %*% Omega_v_C22 %*% F_mean)))^(-1) * sum(diag(Omega_v_C12 %*% F_mean + t(F_mean) %*% Omega_v_C21));
  
  VR_inv <- solve(VR); VC_inv <- solve(VC)
  
  if(alpha == 0)
  {
    Sigma_Ri <- VR_inv %*% Omega_v_R11 %*% VR_inv
    Sigma_Cj <- VC_inv %*% Omega_v_C11 %*% VC_inv
  }
  else
  {
    Sigma_Ri = VR_inv %*% (Omega_v_R11 + alpha*Omega_v_R12%*%t(F_mean) + alpha*F_mean%*%Omega_v_R21 + alpha^2*F_mean%*%Omega_v_R22%*%t(F_mean)) %*% VRinv
    Sigma_Cj = VC_inv %*% (Omega_v_C11 + alpha*Omega_v_C12%*%F_mean + alpha*t(F_mean)%*%Omega_v_C21 + alpha^2*t(F_mean)%*%mega_v_C22%*%F_mean) %*% VCinv
  }
  
  return(list(Ri = Sigma_Ri, Cj = Sigma_Cj, alpha_1 = alpha_1, alpha_2 = alpha_2))
}