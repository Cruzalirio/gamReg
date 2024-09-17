#' @title 
#' Estimation of spatial models with 
#' heteroskedastic innovations with GAMLSS
#' 
#' @description
#' A set 
#' 

gamsarar <- function(x,z,y,wt1,wt2,threshold = 1e-4, max_iter = 200){
  diff_lambda <- 10000
  diff_rho <- 10000 
  iter_count <- 0
  lambda_fin <- 0
  n <- nrow(x)
  rho_fin <- 0
  B_f <- (diag(length(y))-lambda_fin*wt2)
  A_f <- (diag(length(y))-rho_fin*wt1)
  while(diff_lambda > threshold & diff_rho> threshold){
    
    m1 <- gamlss::gamlss(B_f%*%A_f%*%y~B_f%*%x-1,~z-1)
    
    var2 <- diag(c(predict(m1,what="sigma",type="response")^2))
    
    vero_ini <- function(rhos,x, z, y,wt1, wt2, var2,n){
      
      beta_ini <- solve(t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%x)%*%t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%(diag(length(y))-rhos[1]*wt1)%*%y #paso 1
      vv <- t((diag(length(y))-rhos[1]*wt1)%*%y-x%*%beta_ini)%*%
        t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%
        (diag(length(y))-rhos[2]*wt2)%*%((diag(length(y))-rhos[1]*wt1)%*%y-x%*%beta_ini)
      ve <- var2;L <- chol(ve)
      ldetW <- 2*sum(diag(log(L))) 
      lc <- -n/2*log(pi)-1/2*ldetW+log(det(diag(length(y))-rhos[2]*wt2))+
        log(det(diag(length(y))-rhos[1]*wt1))-1/2*vv
      return(-lc)
    }
    
    nC <- parallel::detectCores()
    cl <- parallel::makeCluster(nC-1) 
    parallel::setDefaultCluster(cl=cl)
    p_ini <- optimParallel(par=c(0,0),vero_ini,
                         x=x, y=y,z=z, wt1=wt1, wt2=wt2,var2=var2,n=n,
                         lower=rep(-0.99999, 2), upper = rep(0.99999, 2),
                         parallel=list(loginfo=TRUE))
    
    parallel::setDefaultCluster(cl=NULL)
    parallel::stopCluster(cl)
    
    lambda_ini <- p_ini$par[2]
    rho_ini <- p_ini$par[1]
    
    B <- (diag(length(y))-lambda_ini*wt2)
    A <- (diag(length(y))-rho_ini*wt1)
    
    m2 <- gamlss(B%*%A%*%y~B%*%x-1,~z-1)
    
    var2 <- diag(c(predict(m2,what="sigma",type="response")^2))
    
    vero_fin <- function(rhos, x, z, y,wt1, wt2, var2, n){
      
      beta_ini <- solve(t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%x)%*%t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%(diag(length(y))-rhos[1]*wt1)%*%y 
      vv <- t((diag(length(y))-rhos[1]*wt1)%*%y-x%*%beta_ini)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%((diag(length(y))-rhos[1]*wt1)%*%y-x%*%beta_ini)
      ve <- var2
      L <- chol(ve)
      ldetW <- 2*sum(diag(log(L))) 
      lc <- -n/2*log(pi)-1/2*ldetW+log(det(diag(length(y))-rhos[2]*wt2))+
        log(det(diag(length(y))-rhos[1]*wt1))-1/2*vv
      return(-lc)
    }
    
    nC <- parallel::detectCores()
    cl <- parallel::makeCluster(nC-1) 
    parallel::setDefaultCluster(cl=cl)
    p_fin <- optimParallel(par=c(p_ini$par[1],p_ini$par[2]),vero_ini,
                         x=x, y=y,z=z, wt1=wt1, wt2=wt2,var2=var2,n=n,
                         lower=rep(-0.99999, 2), upper = rep(0.99999, 2),
                         parallel=list(loginfo=TRUE))
    
    parallel::setDefaultCluster(cl=NULL)
    parallel::stopCluster(cl)
    
    lambda_fin <- p_fin$par[2]
    rho_fin <- p_fin$par[1]
    
    diff_lambda <- lambda_fin-lambda_ini
    diff_rho <- rho_fin-rho_ini
    B_f <- (diag(length(y))-lambda_fin*wt2)
    A_f <- (diag(length(y))-rho_fin*wt1)
    
    m_f <- gamlss::gamlss(B_f%*%A_f%*%y~B_f%*%x-1,~z-1)
    
    iter_count <- iter_count + 1
    if(iter_count > max_iter) {
      stop("No converge")
    }
    
  }
  
  coef <- c("(Intercept)" = m_f$mu.coefficients[1], x1 = m_f$mu.coefficients[2], m_f$mu.coefficients[3], 
           "(rho)"=rho_fin[1],"(lambda)"=lambda_fin[1],
           "alpha_0"=m_f$sigma.coefficients[1], "alpha_1"=m_f$sigma.coefficients[2], "alpha_2"=m_f$sigma.coefficients[3])
  return(coef)
}
