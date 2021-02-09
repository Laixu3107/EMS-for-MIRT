#' ----------------------------------
#' title : EML1 for MIRT
#' author: Laixu3107
#' date  : 2020.12.17
#' ----------------------------------


if(sys.nframe()==0L){rm(list=ls()); gc()}
# ---- Needed Packages ----
library(glmnet)     # L1-logistic regression
library(mvtnorm)    # dmvnorm
library(progress)   # progress bar

# ---- Some useful function ----
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


# -------- EML1 MIRT function --------
EM_MIRT <- function(y,            # data set, all responses of all subjects
                    mu_t,         # true mu (mean vector of latent traits)
                    sigma_t,      # true sigma (covariance of latent traits)
                    fixed,        # which item not do sub-model selection
                    lambda_list,  # penalized parameters list
                    A_init,       # initial value of A for EM algorithm
                    b_init,       # initial value of b for EM algorithm
                    grid_num      # number of grid points for each latent trait
){
  
  # ---- Matrix dimension ----
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent traits
  N <- nrow(y)       # number of subjects
  
  # ---- Grid points ----
  x_grid <- Grid_pts(K=K, lb=-4, ub=4, np=grid_num)
  M      <- nrow(x_grid)
  
  y_samp <- apply(y, MARGIN=2, FUN=function(x){rep(x,each=M,times=1)})
  p_samp <- rep(NA, times=N*M)
  x_samp <- matrix(0,N*M,K)    
  for(i in 1:N){
    x_samp[(1:M)+(i-1)*M,] <- x_grid
  }
  
  # ---- Simulation settings ----
  nla    <- length(lambda_list)
  b_list <- list()
  A_list <- list()
  for(ila in 1:nla){
    b_list[[ila]] <- rep(0, J)
    A_list[[ila]] <- matrix(data=0, nrow=K, ncol=J)
  }
  bic_vec      <- rep(Inf, nla)
  cr_vec       <- rep(0, nla)
  time_em1_vec <- rep(0, nla)
  time_em2_vec <- rep(0, nla)
  iter_em1_vec <- rep(0, nla)
  iter_em2_vec <- rep(0, nla)
  
  bics <- rep(Inf, J)
  
  # ---- Simulation ----
  # count <- 0
  time_total <- proc.time()
  
  for(ila in 1:nla){
    
    cat(sprintf("---- ila: %02d ----\n", ila))

    # ---- EM L1 initial ----
    if(ila==1){
      A_c   <- A_init
      Mod_c <- A_c!=0
      b_c   <- b_init
      
      A_new   <- A_c
      Mod_new <- Mod_c
      b_new   <- b_c
    }
    
    # ---- EM L1 iteration ----
    t_em1 <- proc.time()
    iter_em1 <- 0
    
    while(iter_em1 < 100){
      
      iter_em1 <- iter_em1 + 1
      pb_em <- progress_bar$new(format=sprintf("la:%02d EM-L1:%03d [:bar] :percent eta::eta", ila, iter_em1),
                                total=J+1, clear=TRUE, width=60, show_after=0)
      pb_em$tick(0)  # progress bar
      
      # ---- E-step: ----
      p_samp <- Compute_p_tilde(y, x_grid, A_c, b_c, mu_t, sigma_t)
      pb_em$tick(1)
      
      # ---- M-step: ----
      for(j in fixed){
        
        excl <- which(A_init[,j]==0)
        
        fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                      alpha=1, lambda=0, standardize=FALSE, exclude=excl)
        
        bics_j    <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
        A_new[,j] <- as.vector(coef(fit)[-1])
        b_new[j]  <- as.numeric(coef(fit)[1])
        bics[j]   <- bics_j
        
        pb_em$tick(1)
      }
      
      for(j in (1:J)[-fixed]){
        
        fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                      alpha=1, lambda=lambda_list[ila], standardize=FALSE)
        
        bics_j    <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
        A_new[,j] <- as.vector(coef(fit)[-1])
        b_new[j]  <- as.numeric(coef(fit)[1])
        bics[j]   <- bics_j
        pb_em$tick(1)
      }
      
      Mod_new <- A_new!=0
      # ---- Display the new parameters ----
      cat("A_new:\n");   print(A_new)
      
      # ---- Stop criterion ----
      if(all(Mod_new==Mod_c)){
        err <- max(abs(A_new[Mod_new]-A_c[Mod_c])/abs(A_c[Mod_c])) # maximum relative difference
        cat("err:", err, "\n")
        if(err < 0.1){
          break
        }
      }
      else{
        cat("Models wasn't same.\n")
      }
      
      # ---- Replace current parameter ----
      A_c   <- A_new
      b_c   <- b_new
      Mod_c <- Mod_new
      
    } # end while EM-L1
    
    t_em1 <- proc.time() - t_em1
    cat("t_em1:", as.numeric(t_em1[3]), "seconds.\n")
    
    time_em1_vec[ila] <- as.numeric(t_em1[3])
    iter_em1_vec[ila] <- iter_em1
    
    
    # ---- EM parameters estimate ----
    # ---- Note that: if the model is same to previous lambda, skip this one. 
    if(ila > 1){
      if(all(Mod_new==(A_list[[ila-1]]!=0))){
        
        A_list[[ila]] <- A_list[[ila-1]]
        b_list[[ila]] <- b_list[[ila-1]]
        bic_vec[ila]  <- bic_vec[ila-1]
        next
      }
    }
    
    A_c <- A_new
    b_c <- b_new
    
    t_em2 <- proc.time()
    iter_em2 <- 0
    while(iter_em2 < 100){
      
      iter_em2 <- iter_em2 + 1
      pb_em <- progress_bar$new(format=sprintf("la:%02d EM:%03d [:bar] :percent eta::eta", ila, iter_em2),
                                total=J+1, clear=TRUE, width=60, show_after=0)  
      pb_em$tick(0)  # progress bar
      
      # ---- E-step: ----
      p_samp <- Compute_p_tilde(y=y, x=x_grid, A=A_c, b=b_c, mu=mu_t, sigma=sigma_t)
      pb_em$tick(1)
      
      # ---- M-step: ----
      for(j in 1:J){
        
        if(sum(A_c[,j]==0)==K){
          # -- if jth column are all zeros, skip this estimation --
          A_new[,j] <- A_c[,j]
          b_new[j]  <- b_c[j]
        }
        else{
          excl <- which(A_c[,j]==0)
          fit <- glmnet(x=x_samp, y=y_samp[,j], weights=p_samp, family="binomial",
                        alpha=1, lambda=0,  standardize=FALSE, exclude=excl)
          
          bics_j    <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
          A_new[,j] <- as.vector(coef(fit)[-1])
          b_new[j]  <- as.numeric(coef(fit)[1])
          bics[j]   <- bics_j
        }
        
        pb_em$tick(1)
      }
      
      # ---- Display the new parameters ----
      cat("A_new:\n");   print(A_new)
      
      # ---- Stop criterion ----
      # -- Note: Why we check whether the models same?
      #          One time, the non-zero parameter turn to zero through estimate.
      
      if(all((A_new!=0)==(A_c!=0))){
        err <- max(abs(A_new[A_new!=0]-A_c[A_c!=0])/abs(A_c[A_c!=0])) # maximum relative difference
        cat("err:", err, "\n")
        if(err < 0.1){
          break
        }
      }
      else{
        cat("Models wasn't same.\n")
      }
      
      # ---- Replace current parameter ----
      A_c <- A_new
      b_c <- b_new
      
    } # end while EM
    
    t_em2 <- proc.time() - t_em2
    cat("em2_use:", as.numeric(t_em2[3]), "seconds.\n")
    
    time_em2_vec[ila] <- as.numeric(t_em2[3])
    iter_em2_vec[ila] <- iter_em2
    
    
    # ---- Calculate bic ----
    # bic <- Compute_Bic(y, x_grid, A_new, b_new, prob_samp)
    bic <- sum(bics)
    
    # ---- Collect result ----
    A_list[[ila]] <- A_new
    b_list[[ila]] <- b_new
    bic_vec[ila]  <- bic
    
    A_t_sub       <- A_t[,-fixed]
    A_opt_sub     <- A_list[[ila]][,-fixed]
    cr_vec[ila]   <- sum((A_t_sub!=0)==(A_opt_sub!=0))/(J-K)/K
    
  } # end for(ila in 1:nla)


  # ---- Model selection ----
  time_total <- proc.time() - time_total
  time_total <- as.numeric(time_total[3])
  cat("use_total:", time_total, "\n")
  
  bic_min <- which.min(bic_vec)
  A_opt <- A_list[[bic_min]]
  b_opt <- b_list[[bic_min]]
  
  # ---- Return output ----
  result <- list(A_list=A_list,
                 b_list=b_list,
                 bic_vec=bic_vec,
                 bic_min=bic_min,
                 time_em1_vec=time_em1_vec,
                 time_em2_vec=time_em2_vec,
                 iter_em1_vec=iter_em1_vec,
                 iter_em2_vec=iter_em2_vec,
                 time_total=time_total,
                 A_opt=A_opt,
                 b_opt=b_opt
  )
  return(result)
  
} # end EML1 MIRT function


# ------------------------------------------------------------------------------
# ---- A simple test ----
if(sys.nframe() == 0L){
  library(magrittr)   # %>% 
  library(MASS)       # mvrnorm
  
  # ---- load true model ----
  A_t <- c(1,1,1,0,0,0,
           0,0,0,1,1,1)
  A_t <- matrix(A_t, nrow=2, ncol=6, byrow=T)
  fixed <- c(1,4)
  
  J <- ncol(A_t)  # no. of items
  K <- nrow(A_t)  # no. of latent traits
  N <- 1000       # no. of subjects
  
  # ---- true parameter setting ----
  b_t     <- rep(0,J)
  mu_t    <- rep(0,K)
  sigma_t <- matrix(.1,K,K); diag(sigma_t) <- 1
  
  # ---- generate random sample ----
  set.seed(1)
  x <- mvrnorm(n=N, mu=mu_t, Sigma=sigma_t)   # latent variable, here just only to generate y
  y <- x %>%
    `%*%` (A_t) %>%
    `+` (matrix(data=b_t,nrow=N,ncol=J,byrow=T)) %>%
    plogis(q=.) %>%
    rbinom(n=N*J, size=1, prob=.) %>%
    matrix(data=., nrow=N, ncol=J, byrow=F)
  
  # ---- EM initial parameters ----
  A_init <- matrix(data=1/J, nrow=K, ncol=J, byrow=TRUE); A_init[,fixed] <- A_t[,fixed]
  b_init <- rep(0,J)
  grid_num <- 11
  lambda_list <- c(0.01, 0.05, 0.1, 0.15)
  # ---- EM L1 ----
  EM_MIRT(y=y,               
          mu_t=mu_t, 
          sigma_t=sigma_t,
          fixed=fixed,
          lambda_list,
          A_init=A_init,
          b_init=b_init,
          grid_num=grid_num
  )
}
