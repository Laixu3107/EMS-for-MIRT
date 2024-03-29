#' -----------------------------------------------------------------------------
#' title :  Variable selection for M2PL models by EMS with known Sigma
#' author:  Laixu3107
#' date  :  2020.12.20
#' -----------------------------------------------------------------------------


if(sys.nframe() == 0L){rm(list=ls()); gc()}
# ---- Necessary packages ----
library(glmnet)     # logistic regression
library(mvtnorm)    # dmvnorm
library(progress)   # progress bar


# ---- Some useful functions ----
Grid_pts <- function(K, lb=-4, ub=4, np=11){
  # K     : no. of latent variables.
  # lb, ub: lower bound, upper bound.
  # np    : no. of grid points for each variable.
  # output: grid points.
  
  p_list <- list()
  for(k in 1:K){
    p_list[[k]] <- seq(from=lb, to=ub, length.out=np)
  }
  grid_pts <- as.matrix(expand.grid(p_list, stringsAsFactors=FALSE))
  colnames(grid_pts) <- NULL
  return(grid_pts)
}


Compute_p_tilde <- function(y, x_grid, A, b, mu, sigma){
  # y        : responses of all subjects, a matrix with N*J.
  # x_grid   : grid points of x, note that all subjects have same grid points.
  # A, b     : current item parameters.
  # mu, sigma: parameters of p.d.f. of x grid.
  # return   : p tilde.
  
  M <- nrow(x_grid)
  N <- nrow(y)
  
  log_px  <- dmvnorm(x_grid, mean=mu, sigma=sigma, log=T)
  logit_p <- x_grid%*%A + rep(1, M)%*%t(b)
  p       <- plogis(logit_p)
  
  prob <- rep(0, N*M)
  for(i in 1:N){
    log_pyx <- rowSums(dbinom(rep(1, M)%*%t(y[i, ]), 1, p, log=T))
    loglik  <- log_px + log_pyx
    lik     <- exp(loglik)
    prob[1:M+(i-1)*M] <- lik/sum(lik)
  } 
  return(prob)
}


List_Candidate_Models <- function(K){
  # List all candidate sub-models.
  # K     : number of latent traits.
  # output: A matrix, each row denotes a model.
  
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


# -------- EMS M2PL function with known Sigma --------
EMS_MIRT <- function(y,            # data set, responses of all subjects
                     mu_t,         # true mu (mean vector of latent traits)
                     sigma_t,      # true sigma (covariance of latent traits)
                     fixed,        # which item not do model selection
                     A_init,       # initial value of A for EMS algorithm
                     b_init,       # initial value of b for EMS algorithm
                     grid_num      # number of grid points for each latent trait
){
  
  # ---- Dimension of A & y ----
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent traits
  N <- nrow(y)       # number of subjects
  
  # ---- Grid points & sub-models----
  x_grid <- Grid_pts(K=K, lb=-4, ub=4, np=grid_num)
  M      <- nrow(x_grid)
  
  y_samp <- apply(y, MARGIN=2, FUN=function(x){rep(x,each=M,times=1)})
  p_samp <- rep(NA, times=N*M)
  x_samp <- matrix(0,N*M,K)    
  for(i in 1:N){
    x_samp[(1:M)+(i-1)*M,] <- x_grid
  }
  
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
  A_c   <- A_init;  A_new   <- A_c
  Mod_c <- A_c!=0;  Mod_new <- Mod_c
  b_c   <- b_init;  b_new   <- b_c
  
  # ---- EMS iterations ----
  iter <- 0
  bic_seq <- c()  # expected bic sequence generated by EMS algorithm
  time_total <- proc.time()
  while(iter < 100){
    
    iter <- iter + 1
    
    pb_ems <- progress_bar$new(format=sprintf("EMS:%03d [:bar] :percent eta::eta", iter),
                               total=(J-K)*NM+1+K, clear=TRUE, width=60, show_after=0)
    pb_ems$tick(0)  # progress bar
    
    # ---- E-step: ----
    p_samp <- Compute_p_tilde(y, x_grid, A_c, b_c, mu_t, sigma_t)
    pb_ems$tick(1)
    
    # ---- MS-step: ----
    # -- Estimation without penalty for each sub-model using glmnet --
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
    bic_min <- sum(bic_vec)
    bic_seq <- c(bic_seq, bic_min)
    
    # ---- Display the new parameters ----
    cat("A_new:\n");    print(A_new)
    cat("bic.min:", bic_min, "\t")
    
    # ---- Stop criterion ----
    if(all(Mod_new==Mod_c)){
      err <- max(abs(A_new[Mod_new]-A_c[Mod_c])/abs(A_c[Mod_c])) # maximum relative difference
      cat("err:", err, "\n")
      if(err < 0.1){break}
    }
    else{
      cat("Models wasn't same.\n")
    }
    
    # ---- Replace the current model & parameters ----
    A_c   <- A_new
    b_c   <- b_new
    Mod_c <- Mod_new
    
  } # end while
  
  time_total <- proc.time() - time_total
  time_total <- as.numeric(time_total[3])
  cat("use.time:", time_total, "\n")
  
  # ---- Return output ----
  result <- list(A_opt=A_new,
                 b_opt=b_new,
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
  A_t <- c(1,1,1,0,0,0,
           0,0,0,1,1,1)
  A_t <- matrix(A_t, nrow=2, ncol=6, byrow=T)
  fixed <- c(1,4)
  
  J <- ncol(A_t)  # no. of items
  K <- nrow(A_t)  # no. of latent traits
  N <- 1000       # no. of subjects
  
  b_t     <- rep(0,J)
  mu_t    <- rep(0,K)
  sigma_t <- matrix(.1,K,K); diag(sigma_t) <- 1
  
  # ---- Generate random samples ----
  set.seed(1)
  x <- mvrnorm(n=N, mu=mu_t, Sigma=sigma_t)   # latent variable, here just only to generate y
  y <- x %>%
    `%*%` (A_t) %>%
    `+` (matrix(data=b_t,nrow=N,ncol=J,byrow=T)) %>%
    plogis(q=.) %>%
    rbinom(n=N*J, size=1, prob=.) %>%
    matrix(data=., nrow=N, ncol=J, byrow=F)
  
  # ---- EMS initial model & parameters ----
  A_init <- matrix(data=1/J, nrow=K, ncol=J, byrow=TRUE); A_init[,fixed] <- A_t[,fixed]
  b_init <- rep(0,J)
  grid_num <- 11
  
  # ---- EMS ----
  EMS_MIRT(y=y,               # data set, responses of all subjects
           mu_t=mu_t,         # true mu (mean vector of latent traits)
           sigma_t=sigma_t,   # true sigma (covariance of latent traits)
           fixed=fixed,       # which item not do model selection
           A_init=A_init,     # initial value of A for EMS algorithm
           b_init=b_init,     # initial value of b for EMS algorithm
           grid_num=grid_num  # number of grid points each dimension
  )
}
