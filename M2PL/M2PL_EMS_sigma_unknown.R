#' -----------------------------------------------------------------------------
#' title:  Variable selection for M2PL models by EMS with unknown Sigma
#' author: Laixu3107
#' date:   2021.04.15
#' -----------------------------------------------------------------------------


if(sys.nframe() == 0L){rm(list=ls()); gc()}
# ---- Necessary packages ----
library(glmnet)     # logistic regression
library(mvtnorm)    # dmvnorm
library(progress)   # progress bar

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("./esti_sigma.cpp")


# ---- Some useful functions ----
Grid_pts <- function(K, lb=-4, ub=4, np=11){
  # K     : no. of latent variables.
  # lb, ub: lower bound, upper bound.
  # np    : no. of grid points for each dimension.
  # ouput : grid points.
  
  p_list <- list()
  for(k in 1:K){
    p_list[[k]] <- seq(from=lb, to=ub, length.out=np)
  }
  grid_pts <- as.matrix(expand.grid(p_list, stringsAsFactors=FALSE))
  colnames(grid_pts) <- NULL
  return(grid_pts)
}


Compute_p_tilde <- function(y, x_grid, A, b, sigma, iter){
  # Data:   2019.09.28.   Author: Laixu3107.
  # y:      responses of all examinees, a matrix. 
  # x_grid: grid points of x.
  # return: probabilities or weights for all grid points.
  
  G <- nrow(x_grid)
  N <- nrow(y)
  
  pb_estep <- progress_bar$new(format=sprintf("EMS:%03d E-step [:bar] :percent eta::eta", iter),
                               total=N+1, clear=TRUE, width=60, show_after=0)
  pb_estep$tick(0)  # progress bar
  
  log_phi_x <- dmvnorm(x_grid, mean=rep(0,K), sigma=sigma, log=T)
  p2y1      <- plogis(x_grid%*%A + rep(1, G)%*%t(b))
  log_p2y1  <- log(p2y1)
  log_p2y0  <- log(1-p2y1)
  
  log_p2y <- matrix(data=0,nrow=G,ncol=J)
  tp <- rep(0, N*G) # tilde p(xi| yi, A, b, sigma)
  wg <- rep(0, G)   # weights for sigma estimation
  
  pb_estep$tick(1)
  
  for(i in 1:N){
    
    log_p2y[, y[i,]==1] <- log_p2y1[, y[i,]==1]
    log_p2y[, y[i,]==0] <- log_p2y0[, y[i,]==0]
    
    log_p2yx <- rowSums(log_p2y) + log_phi_x
    p2yx <- exp(log_p2yx)
    
    tp[1:G+(i-1)*G] <- p2yx/sum(p2yx)
    wg <- wg + tp[1:G+(i-1)*G]
    
    pb_estep$tick(1)
  }
  
  return(list(tp=tp, wg=wg))
}


List_Candidate_Models <- function(K){
  # Laixu3107, Dec 16, 2020.
  # List all candidate sub-models.
  # K:      number of latent traits.
  # Output: A matrix, each row denotes a sub-model.
  
  models <- matrix(data=0, nrow=2^K, ncol=K)
  s <- 0
  for(k in 0:K){
    com_k <- t(combn(K, k))
    for(l in 1:nrow(com_k)){
      s <- s + 1
      models[s, com_k[l,]] <- 1 
    }
  }
  return(models)
}


Calcu_wscov <- function(x_grid, wg, N){
  # Laixu3107, Apr 22, 2021
  # Calculate weighted sample covariance and correlation matrices.
  G <- nrow(x_grid)
  K <- ncol(x_grid)
  
  scov <- matrix(data=0,nrow=K,ncol=K)
  for(g in 1:G){
    scov <- scov + x_grid[g,]%*%t(x_grid[g,]) *wg[g]
  }
  scov <- scov/N
  
  scor <- matrix(data=0,nrow=K,ncol=K)
  for(k in 1:K){scor[k,] <- scov[k,]/scov[k,k]}
  for(k in 1:K){scor[,k] <- scov[,k]/scov[k,k]}
  
  scov <- as.matrix(forceSymmetric(scov))
  scor <- as.matrix(forceSymmetric(scor))
  
  return(list(wscov=scov, wscor=scor))
}


Culcu_neg2loglik <- function(S, sigma, N){
  N*log(det(sigma)) + N*sum(diag(solve(sigma)%*%S))
}


# -------- EMS for M2PL with unknown Sigma --------
EMS_MIRT <- function(y,            # data set, all responses of all subjects
                     A_init,       # initial value of A     for EMS algorithm
                     b_init,       # initial value of b     for EMS algorithm
                     sigma_init,   # initial value of sigma for EMS algorithm
                     fixed         # which item not do model selection
){
  
  # ---- Dimension of A & y ----
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent traits
  N <- nrow(y)       # number of subjects
  
  # ---- Grid points ----
  x_grid <- Grid_pts(K=K, lb=-4, ub=4, np=11)
  G      <- nrow(x_grid)
  
  y_samp <- apply(y, MARGIN=2, FUN=function(x){rep(x,each=G,times=1)})
  storage.mode(y_samp) <- "integer"
  
  p_samp <- rep(0, times=N*G)
  x_samp <- matrix(data=0,nrow=N*G,ncol=K)
  for(i in 1:N){
    x_samp[(1:G)+(i-1)*G,] <- x_grid
  }
  
  # ---- Sub-models ----
  Models <- List_Candidate_Models(K)
  NM     <- nrow(Models)
  
  # ---- Settings ----
  coef_mat  <- matrix(data=0, nrow=NM, ncol=K+1)
  bic_j     <- rep(Inf, NM) # bic of all models for each item
  opt_vec   <- rep(0, J)    # opt model for each item
  bic_vec   <- rep(Inf, J)  # opt bic for each item
  bic_min   <- Inf          # min bic in MS-step
  
  BICs_list <- list()
  Coef_list <- list()
  for(j in 1:J){
    Coef_list[[j]] <- coef_mat
    BICs_list[[j]] <- bic_j
  }
  
  # ---- EMS initialization ----
  A_c     <- A_init;      A_new     <- A_c
  Mod_c   <- A_c!=0;      Mod_new   <- Mod_c
  b_c     <- b_init;      b_new     <- b_c
  sigma_c <- sigma_init;  sigma_new <- sigma_c
  
  # ---- EMS iterations ----
  iter <- 0
  bic_seq <- rep(0,100)  # expected bic sequence generated by EMS
  time_total <- proc.time()
  while(iter < 100){
    
    iter <- iter + 1
    
    # ---- E-step: ----
    estep  <- Compute_p_tilde(y, x_grid, A_c, b_c, sigma_c, iter=iter)
    p_samp <- estep$tp
    wg     <- estep$wg
    
    # ---- MS-step: ----
    pb_ems <- progress_bar$new(format=sprintf("EMS:%03d MS-step [:bar] :percent eta::eta", iter),
                               total=(J-K)*NM+1+K, clear=TRUE, width=60, show_after=0)
    pb_ems$tick(0)  # progress bar
    
    # -- Estimate the correlation matrix --
    wscov_wscor <- Calcu_wscov(x_grid, wg, N)
    wscov <- wscov_wscor$wscov
    wscor <- wscov_wscor$wscor
    
    sigma_new  <- calcu_sigma_cmle_cpp(scov=wscov, scor=wscor, tol=1e-6)
    neg2loglik <- Culcu_neg2loglik(S=wscov, sigma=sigma_new, N)
    
    pb_ems$tick(1)
    
    # -- Estimate A and b for each sub-model by using glmnet --
    # -- for fixed items:
    for(j in fixed){
      
      exclude <- which(A_init[, j]==0)
      
      fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                    alpha=1, lambda=0, standardize=FALSE, exclude=exclude)
      
      A_new[, j] <- as.vector(coef(fit))[-1]
      b_new[j]   <- as.vector(coef(fit))[1]
      bic_vec[j] <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
      
      pb_ems$tick(1)
      
    }
    
    # -- for non-fixed items:
    for(j in (1:J)[-fixed]){
      
      for(iM in 1:NM){
        
        exclude <- which(Models[iM, ]==0)
        if(length(exclude) == K){
          fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                        alpha=1, lambda=1e4, standardize=FALSE)
        }
        else{
          fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                        alpha=1, lambda=0, standardize=FALSE, exclude=exclude)
        }
        
        coef_mat[iM,] <- as.vector(coef(fit))
        bic_j[iM]     <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
        
        pb_ems$tick(1)
      }
      
      Coef_list[[j]] <- coef_mat
      BICs_list[[j]] <- bic_j
      
      opt_j      <- which.min(bic_j)
      opt_vec[j] <- opt_j
      bic_vec[j] <- bic_j[opt_j]
      
      A_new[, j] <- coef_mat[opt_j, -1]
      b_new[j]   <- coef_mat[opt_j, 1]
      
    }
    
    Mod_new <- A_new!=0
    bic_min <- sum(bic_vec) + neg2loglik
    bic_seq[iter] <- bic_min
    
    # ---- Display the new parameters ----
    # cat("A_new:\n");      print(A_new);     cat("\n");
    # cat("b_new:\n");      print(b_new);     cat("\n")
    # cat("sigma_new:\n");  print(sigma_new); cat("\n")
    # cat("bic.min:", bic_min, "\t")
    # cat("neg2loglik:", neg2loglik, "\t")
    # cat("bic_Ab:", sum(bic_vec), "\n")
    
    # ---- Stop criterion ----
    if(all(Mod_new==Mod_c)){
      err_A <- max(abs((A_new[Mod_new]-A_c[Mod_c])/A_c[Mod_c])) # maximum relative difference
      err_sigma <- max(abs((sigma_new-sigma_c)/sigma_c))
      # cat("err_A:", err_A, "\n")
      # cat("err_sigma:", err_sigma, "\n")
      if(all(c(err_A, err_sigma) < 0.1)) {break}
      
    }
    # else{
    #   cat("Models wasn't same.\n")
    # }
    
    # ---- Replace the current model & parameters ----
    A_c     <- A_new
    b_c     <- b_new
    Mod_c   <- Mod_new
    sigma_c <- sigma_new
    
  } # end while
  
  bic_seq <- bic_seq[1:iter]
  
  time_total <- proc.time() - time_total
  time_total <- as.numeric(time_total[3])
  # cat("use.time:", time_total, "\n")
  
  # ---- Return output ----
  result <- list(A_opt=A_new,
                 b_opt=b_new,
                 sigma_opt=sigma_new,
                 bic_seq=bic_seq,
                 iter_ems=iter,
                 time_ems=time_total
  )
  return(result)
  
}


# ------------------------------------------------------------------------------
if(sys.nframe() == 0L){
  
  library(magrittr)
  library(MASS)
  
  # ---- True model & parameters setting ----
  A_t <- c(2,1.5,1,0,0,0,0,0,0,
           0,0,0,2,1.5,1,0,0,0,
           0,0,0,0,0,0,2,1.5,1)
  A_t <- matrix(A_t, nrow=3, ncol=9, byrow=T)
  fixed <- c(1,4,7)
  
  J <- ncol(A_t)  # no. of items
  K <- nrow(A_t)  # no. of latent traits
  N <- 500        # no. of subjects
  
  b_t     <- rep(0,J)
  sigma_t <- matrix(.1,K,K); diag(sigma_t) <- 1

  # ---- Generate random samples ----
  set.seed(1)
  x <- mvrnorm(n=N, mu=rep(0,K), Sigma=sigma_t) # latent traits
  y <- x %>%
    `%*%` (A_t) %>%
    `+` (matrix(data=b_t,nrow=N,ncol=J,byrow=T)) %>%
    plogis(q=.) %>%
    rbinom(n=N*J, size=1, prob=.) %>%
    matrix(data=., nrow=N, ncol=J, byrow=F)
  
  # ---- EMS initial model & parameters ----
  A_init <- matrix(data=1/J, nrow=K, ncol=J, byrow=TRUE); A_init[,fixed] <- A_t[,fixed]
  b_init <- rep(0,J)
  sigma_init <- matrix(0.5,K,K);  diag(sigma_init) <- 1
  
  # ---- EMS ----
  EMS_MIRT(y=y,                    # data set, all responses of all subjects
           A_init=A_init,          # initial value of A     for EMS algorithm
           b_init=b_init,          # initial value of b     for EMS algorithm
           sigma_init=sigma_init,  # initial value of sigma for EMS algorithm
           fixed=fixed             # which item not do model selection
  )
}
