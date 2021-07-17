#' -----------------------------------------------------------------------------
#' title:  Variable selection for M2PL models by EMSSP with unknown Sigma 
#' author: Laixu3107
#' date:   2021.04.16
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


List_Candidate_Models <- function(K){
  # Laixu3107, Dec 16, 2020
  # List all candidate sub-models.
  # K:      number of latent traits.
  # Output: A matrix, each row denotes a model.
  
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



E_STEP4 <- function(y, x_grid, A, b, sigma, cumu, pct, iter){
  # Laixu3107, Apr 16, 2021
  # y     : responses of all subjects, a matrix.
  # x_grid: grid points of x, note that all subjects have same grid points.
  # A, b  : current parameters of A, b.
  # sigma : current parameters of sigma.
  # cumu  : cumulative percent of p tilde for each subject.
  # pct   : ratio to number of p tilde and total.
  # return: y_samp_sub, x_samp_sub, p_samp_sub. (and wg)
  
  N <- nrow(y)
  J <- ncol(y)
  G <- nrow(x_grid)
  K <- ncol(x_grid)
  
  pb_estep <- progress_bar$new(format=sprintf("EMS:%03d E-step [:bar] :percent eta::eta", iter),
                               total=2*N+1, clear=TRUE, width=60, show_after=0)
  pb_estep$tick(0)  # progress bar
  
  
  log_phi_x <- dmvnorm(x_grid, mean=rep(0,K), sigma=sigma, log=T)
  p2y1      <- plogis(x_grid%*%A + rep(1,G)%*%t(b))
  log_p2y1  <- log(p2y1)
  log_p2y0  <- log(1-p2y1)
  
  log_p2y <- matrix(data=0,nrow=G,ncol=J)
  tp <- rep(0,G)   # tilde p(xg| yi, A, b, sigma)
  wg <- rep(0,G)   # weights for sigma estimation
  
  # -- work matrix and vector --
  order_mat <- matrix(data=0, nrow=N, ncol=G);storage.mode(order_mat) <- "integer"
  p_ord_mat <- matrix(data=0, nrow=N, ncol=G)
  g_vec <- rep(0,N)
  s_vec <- rep(0,N)
  for(i in 1:N){
    
    log_p2y[, y[i,]==1] <- log_p2y1[, y[i,]==1]
    log_p2y[, y[i,]==0] <- log_p2y0[, y[i,]==0]
    
    log_p2yx <- rowSums(log_p2y) + log_phi_x
    p2yx     <- exp(log_p2yx)
    tp       <- p2yx/sum(p2yx)
    wg       <- wg + tp
    
    order_mat[i, ] <- order(tp, decreasing=TRUE)
    p_ord_mat[i, ] <- tp[order_mat[i, ]]
    
    
    s <- 0
    for(g in 1:G){
      s <- s + p_ord_mat[i, g]
      # print(s);print(s > cumu); print(g/G>pct)
      if(s > cumu & g/G>pct){break}
    }
    g_vec[i] <- g
    s_vec[i] <- s
    pb_estep$tick(1)
  }
  
  # -- y_samp_sub, x_samp_sub and p_samp_sub --
  y_samp_sub <- matrix(data=0, nrow=sum(g_vec), ncol=J); storage.mode(y_samp_sub) <- "integer"
  x_samp_sub <- matrix(data=0, nrow=sum(g_vec), ncol=K)
  p_samp_sub <- rep(0, sum(g_vec))
  
  s <- 0
  for(i in 1:N){
    y_samp_sub[(1:g_vec[i])+s, ] <- matrix(y[i, ], nrow=g_vec[i], ncol=J, byrow=T)
    x_samp_sub[(1:g_vec[i])+s, ] <- x_grid[order_mat[i, 1:g_vec[i]], ]
    p_samp_sub[(1:g_vec[i])+s]   <- p_ord_mat[i, 1:g_vec[i]]/s_vec[i]
    
    s <- s + g_vec[i]
    pb_estep$tick(1)
  }
  
  output <- list(y_samp_sub=y_samp_sub,
                 x_samp_sub=x_samp_sub,
                 p_samp_sub=p_samp_sub,
                 wg=wg
  )
  
  return(output)
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



# -------- EMSSP for M2PL with unknown Sigma --------
EMS_MIRT_sp <- function(y,            # data set, all responses of all subjects
                        A_init,       # initial value of A     for EMS algorithm
                        b_init,       # initial value of b     for EMS algorithm
                        sigma_init,   # initial value of sigma for EMS algorithm
                        fixed,        # which item not do model selection
                        pts_cumu,     # how much points will be use for each subject
                        pts_pct       # how much points will be use for each subject
){
  
  # ---- Dimension of A & y ----
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent traits
  N <- nrow(y)       # number of subjects
  
  # ---- Grid points & sub-models ----
  x_grid    <- Grid_pts(K=K, lb=-4, ub=4, np=11)
  
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
  bic_seq <- c()  # expected bic sequence generated by EMS
  time_total <- proc.time()
  while(iter < 100){
    
    iter <- iter + 1
    
    # ---- E-step: ----
    estep <- E_STEP4(y, x_grid, A_c, b_c, sigma_c, cumu=pts_cumu, pct=pts_pct, iter)
    y_samp_sub <- estep$y_samp_sub
    x_samp_sub <- estep$x_samp_sub
    p_samp_sub <- estep$p_samp_sub
    wg         <- estep$wg
    
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
    
    # -- Estimation without penalty for each sub-model using glmnet --
    # -- for fixed items:
    for(j in fixed){
      
      exclude <- which(A_init[, j]==0)
      
      fit <- glmnet(x=x_samp_sub, y=y_samp_sub[,j], weights=p_samp_sub, family="binomial",
                    alpha=1, lambda=0, standardize=FALSE, exclude=exclude)
      
      A_new[, j] <- as.vector(coef(fit))[-1]
      b_new[j]   <- as.vector(coef(fit))[1]
      bic_vec[j] <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
      
      pb_ems$tick(1)
      
    }
    
    # -- for non-fixed items:
    for(j in (1:J)[-fixed]){
      
      for(iM in 1:NM){
        
        exclude <- which(Models[iM,]==0)
        if(length(exclude) == K){
          fit <- glmnet(x=x_samp_sub, y=y_samp_sub[,j], weights=p_samp_sub, family="binomial",
                        alpha=1, lambda=1e4, standardize=FALSE)
        }
        else{
          fit <- glmnet(x=x_samp_sub, y=y_samp_sub[,j], weights=p_samp_sub, family="binomial",
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
    bic_seq <- c(bic_seq, bic_min)
    
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
    A_c   <- A_new
    b_c   <- b_new
    Mod_c <- Mod_new
    sigma_c <- sigma_new
    
  } # end while
  
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
# ---- A simple test ----
if(sys.nframe() == 0L){
  
  library(magrittr)   # %>% 
  library(MASS)       # mvrnorm
  
  # ---- True model & parameters setting ----
  A_t <- c(2,1.5,1,0,0,0,0,0,0,
           0,0,0,2,1.5,1,0,0,0,
           0,0,0,0,0,0,2,1.5,1)
  A_t <- matrix(A_t, nrow=3, ncol=9, byrow=T)
  fixed <- c(1,4,7)
  
  J <- ncol(A_t)  # no. of items
  K <- nrow(A_t)  # no. of latent traits
  N <- 1000       # no. of subjects
  
  b_t     <- rep(0,J)
  sigma_t <- matrix(.1,K,K); diag(sigma_t) <- 1
  
  # ---- Generate random samples ----
  set.seed(4)
  x <- mvrnorm(n=N, mu=rep(0,K), Sigma=sigma_t) # latent traits
  y <- x %>%
    `%*%` (A_t) %>%
    `+` (matrix(data=b_t,nrow=N,ncol=J,byrow=T)) %>%
    plogis(q=.) %>%
    rbinom(n=N*J, size=1, prob=.) %>%
    matrix(data=., nrow=N, ncol=J, byrow=F)
  
  A_init <- matrix(data=1/J, nrow=K, ncol=J); A_init[,fixed] <- A_t[,fixed]
  b_init <- rep(0,J)
  sigma_init <- matrix(0.5,K,K); diag(sigma_init) <- 1
  
  pts_pct  <- 0.05
  pts_cumu <- 0.90
  
  EMS_MIRT_sp(y=y,                   # data set, all responses of all subjects
              A_init=A_init,         # initial value of A     for EMS algorithm
              b_init=b_init,         # initial value of b     for EMS algorithm
              sigma_init=sigma_init, # initial value of sigma for EMS algorithm
              fixed=fixed,           # which item not do model selection
              pts_cumu=pts_cumu,     # number of grid points each dimension
              pts_pct=pts_pct        # number of grid points each dimension
  )
}
