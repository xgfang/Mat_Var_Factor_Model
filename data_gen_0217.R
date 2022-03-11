require(tensorA)

matvar_normal <- function(p, q, size, mean = NA, cov1 = NA, cov2 = NA)
{
  if(is.na(mean))
  {
    mean = array(0,c(p,q));
  }
  
  if(is.na(cov1))
  {
    cov1 = diag(p);
    A = cov1;
  }
  else
  {
    A = chol(cov1);
  }
  
  if(is.na(cov2))
  {
    cov2 = diag(q);
    B = cov2;
  }
  else
  {
    B = chol(cov2);
  }
  
  Z = array(rnorm(p*q*size),c(p,q,size));
  
  Y =  to.tensor(array(0, c(p,q,size)));
  
  for(i in 1:size)
  {
    Y[,,i] <- mean + A %*% Z[,,i] %*% t(B);
  }
  
  return(Y)
}

VAR_1 <- function(p,q,t,coef)
{
  Y <- to.tensor(array(0,c(p,q,t)));
  
  for(i in 1:p){
    for(j in 1:q){
      Y[i,j,] <- arima.sim(list(order = c(1,0,0), ar = 0.5), n = t);
    }
  }
  
  return(Y)
}

mat_fac_model <- function(setting, p, q, k, r, size)
{
  R <- matrix(runif(p*k,-1,1), nrow = p, ncol = k)
  C <- matrix(runif(q*r,-1,1), nrow = q, ncol = r)
  
  if(setting == 1)
  {
    E = matvar_normal(p, q, size)
    Fac = matvar_normal(k, r, size)
  }
  
  else if(setting == 2)
  {
    E = VAR_1(p,q,size,coef = 0.5)
    Fac = VAR_1(k,r,size,coef = 0.1)
  }
  else if(setting == 3)
  {
    Ecov_R <- matrix(1/p,p,p);
    Ecov_C <- matrix(1/q,q,q);
    
    diag(Ecov_R) <- 1; diag(Ecov_C) <- 1
    
    E <- matvar_normal(p,q,size,cov1 = Ecov_R,cov2 = Ecov_C)
    Fac <- matvar_normal(k,r,size)
  }
  
  Y = to.tensor(array(0,dim = c(p,q,size)))
  for(i in 1:size){
    Y[,,i] <- to.tensor(R %*% Fac[,,i] %*% t(C))+ E[,,i];
  }
  
  return(list(Y = Y, Fac = Fac))
  
}

