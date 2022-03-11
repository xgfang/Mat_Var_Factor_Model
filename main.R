source("data_gen_0217.R")
source("est_cov_0217.R")
source("est_RCFS_0217.R")

space_distance <- function(A, B)
{
  if(!all(dim(A) == dim(B)))
    stop('space_distance(): Matrices do have the same shape!')
  
  M <- A%*%solve(crossprod(A))%*%t(A) - B%*%solve(crossprod(B))%*%t(B);
  dis <- base::norm(M,type = "2")
  
  return(dis)
}

one_mc_sim <- function(R, C, alpha, p, q, k, r, size, setting, i, seed)
{
  cat("Process MC #: ", toString(i), "\n")
  set.seed(seed)
  
  dt_gen <- mat_fac_model(setting, p, q, k, r, size)
  Y <- dt_gen$Y; Fac <- dt_gen$Fac
  Y_tilde <- mat_trans(Y, alpha); F_tilde <- mat_trans(Fac, alpha)
  
  est <- estimation_RCFS(Y_tilde, alpha, k, r)
  k_hat <- est$khat; r_hat <- est$rhat; VR_hat <- est$VRhat; VC_hat <- est$VChat
  R_hat <- est$Rhat; C_hat <- est$Chat; F_hat <- est$Fhat
  
  spdistR <- space_distance(R_hat, R);
  spdistC <- space_distance(C_hat, C);
  
  H_RC <- H_matrix(F_tilde, R, C, VR_hat, VC_hat, R_hat, C_hat)
  HR <- H_RC$H_R; HC <- H_RC$H_C
  HR_inv <- solve(HR); HC_inv <- solve(HC)
  
  errorF <- to.tensor(array(NA,dim = dim(F_tilde)))
  for(t in 1:size)
  {
    errorF[,,t] <- F_hat[,,t] - to.tensor(HR_inv%*%Fac[,,t]%*%t(HC_inv));
  }
  errorF_bar <- mean.tensor(errorF, along = 3)
  errorF_bar_norm <- base::norm(errorF_bar, type = "2")
  
  ### may modify here
  row = col = 1;
  
  errorR <- R_hat - R%*%HR
  errorR_row <- errorR[row,]
  
  # Sigma_res <- cov_loading_HAC(Y, R_hat, C_hat, F_hat, VR_hat, VC_hat, row, col, alpha)
  # SigmaR_row <- Sigma_res$Ri
  # SigmaR_row_flat <- matrix(SigmaR_row, nrow = 1)
  
  res <- c(k_hat, r_hat, spdistR, spdistC, errorF_bar_norm, errorR_row)
  return(res)
}

model_serial <- function(alpha, p, q, k, r, size, setting, n_iter)
{
  set.seed(9876)
  R <- matrix(runif(p*k,-1,1), nrow = p, ncol = k)
  C <- matrix(runif(q*r,-1,1), nrow = q, ncol = r)
  
  res_model <- matrix(NA, nrow = n_iter, ncol = 5+k)
  
  for(i_iter in 1:n_iter)
  {
    res_model[i_iter, ] <- one_mc_sim(R, C, alpha, p, q, k, r, size, setting, i = i_iter, seed = i_iter )
  }
  
  return(res_model)
}

#################
## Main Function 
#################

n_iter = 100
alpha = -1
# pqset <- rbind(c(20,20),c(20,100),c(100,100))
p = q = 20
k = r = 3
# tfset = c(0.5,1,1.5,2)
size <- 0.5*p*q
setting = 1

res_model <- array(NA, dim = c(n_iter, 5+k, dim(pqset)[1], length(tfset)))

res_model[,,1,1] <- model_serial(alpha, p, q, k, r, size, setting, n_iter)

# for(pq in 1:dim(pqset)[1])
# {
#   p = pqset[pq, 1]; q = pqset[pq, 2]
#   for(tf in 1:length(tfset))
#   {
#     size = tfset[tf]*p*q
#     res_model[ , , pq, tf] <- model_serial(alpha, p, q, k, r, size, setting, n_iter)
#   }
# }