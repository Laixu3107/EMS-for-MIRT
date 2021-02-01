#' --------------------------------
#' title :  E-MS sp for M2PL model
#' author:  Laixu3107
#' date  :  2020.01.17
#' --------------------------------


if(sys.nframe() == 0L){rm(list=ls()); gc()}
# ---- Needed Packages ----
library(glmnet)     # logistic regression
library(mvtnorm)    # dmvnorm
library(progress)   # progress bar

# ---- Some Useful Function ----
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

E_STEP4 <- function(y, x_grid, A, b, mu, sigma, cumu=0.9, pct=0.05, iter=0){
  # Compute sorted p tilde,
  # y        : responses of all subjects, a matrix with N*J.
  # x_grid   : grid points of x, note that all subjects have same grid points.
  # A, b     : current item parameters.
  # mu, sigma: parameters of p.d.f. of x grid.
  # cumu     : cumulative percent of p tilde for each subject.
  # pct      : percentage of the least grid points used for calculation.
  # iter     : EMS iteration number.
  # return   : y_samp_sub, x_samp_sub, p_samp_sub.
  
  M <- nrow(x_grid)
  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(x_grid)
  
  pb_estep <- progress_bar$new(format=sprintf("EMS:%03d E-step [:bar] :percent eta::eta", iter),
                               total=2*N+1, clear=TRUE, width=60, show_after=0)
  pb_estep$tick(0)  # progress bar
  
  log_px  <- dmvnorm(x_grid, mean=mu, sigma=sigma, log=T)
  logit_p <- x_grid%*%A + rep(1, M)%*%t(b)
  p       <- plogis(logit_p)
  pb_estep$tick(1)
  
  # -- work matrix and vector --
  order_mat <- matrix(data=0, nrow=N, ncol=M)
  p_ord_mat <- matrix(data=0, nrow=N, ncol=M)
  storage.mode(order_mat) <- "integer"

  m_vec <- rep(0, N)
  s_vec <- rep(1, N)
  for(i in 1:N){
    
    pyx <- p
    pyx[, y[i,]==0] <- 1 - p[, y[i,]==0]
    log_pyx <- rowSums(log(pyx)) # 12.37
    
    loglik  <- log_px + log_pyx
    lik     <- exp(loglik)
    prob    <- lik/sum(lik)
    
    order_mat[i, ] <- order(prob, decreasing=TRUE)
    p_ord_mat[i, ] <- prob[order_mat[i, ]]
    
    s <- 0
    for(j in 1:M){
      s <- s + p_ord_mat[i, j]
      if(s > cumu & j/M>pct){break}
    }
    m_vec[i] <- j
    s_vec[i] <- s
    pb_estep$tick(1)
  }
  
  # -- y_samp_sub, x_samp_sub --
  y_samp_sub <- matrix(data=0, nrow=sum(m_vec), ncol=J); storage.mode(y_samp_sub) <- "integer"
  x_samp_sub <- matrix(data=0, nrow=sum(m_vec), ncol=K)
  p_samp_sub <- rep(0, sum(m_vec))
  
  s <- 0
  for(i in 1:N){
    y_samp_sub[(1:m_vec[i])+s, ] <- matrix(y[i, ], nrow=m_vec[i], ncol=J, byrow=T)
    x_samp_sub[(1:m_vec[i])+s, ] <- x_grid[order_mat[i, 1:m_vec[i]], ]
    p_samp_sub[(1:m_vec[i])+s]   <- p_ord_mat[i, 1:m_vec[i]]/s_vec[i]

    s <- s + m_vec[i]
    pb_estep$tick(1)
  }
  
  output <- list(y_samp_sub=y_samp_sub,
                 x_samp_sub=x_samp_sub,
                 p_samp_sub=p_samp_sub
  )
  
  return(output)
}


# -------- main function --------
EMS_MIRT_sp <- function(y,         # data set, responses of all subjects
                        mu_t,      # true mu (mean vector of latent traits)
                        sigma_t,   # true sigma (covariance of latent traits)
                        fixed,     # which item not do model selection
                        A_init,    # initial value of A for EMS algorithm
                        b_init,    # initial value of b for EMS algorithm
                        grid_num,  # number of grid points for each latent variable
                        pts_cumu,  # cumulative percent of p tilde for each subject
                        pts_pct    # percentage of the least grid points used for calculation
){
  
  # ---- Matrix dimension ----
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent traits
  N <- nrow(y)       # number of subjects
  
  # ---- Grid points & sub-models ----
  x_grid <- Grid_pts(K=K, lb=-4, ub=4, np=grid_num)
  
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
  
  # ---- E-MS initialization ----
  A_c   <- A_init;  A_new   <- A_c
  Mod_c <- A_c!=0;  Mod_new <- Mod_c
  b_c   <- b_init;  b_new   <- b_c
  
  # ---- E-MS iteration ----
  iter <- 0
  bic_seq <- c()  # expected bic sequence generated by E-MS algorithm
  time_total <- proc.time()
  while(iter < 100){
    
    iter <- iter + 1
    
    # ---- E-step: ----
    estep <- E_STEP4(y, x_grid, A_c, b_c, mu_t, sigma_t, cumu=pts_cumu, pct=pts_pct, iter)
    y_samp_sub <- estep$y_samp_sub
    x_samp_sub <- estep$x_samp_sub
    p_samp_sub <- estep$p_samp_sub
    # print(length(p_samp_sub))
    # print(length(p_samp_sub)/(N*M))
    
    # ---- MS-step: ----
    # ---- Estimation without penalty for each sub-model using glmnet.
    pb_ems <- progress_bar$new(format=sprintf("EMS:%03d MS-step [:bar] :percent eta::eta", iter),
                               total=(J-K)*NM+K, clear=TRUE, width=60, show_after=0)
    pb_ems$tick(0)  # progress bar
    
    # -- for fixed item:
    for(j in fixed){
      
      exclude <- which(A_init[, j]==0)
      
      fit <- glmnet(x=x_samp_sub, y=y_samp_sub[,j], weights=p_samp_sub, family="binomial",
                    alpha=1, lambda=0, standardize=FALSE, exclude=exclude)
      
      A_new[, j] <- as.vector(coef(fit))[-1]
      b_new[j]   <- as.vector(coef(fit))[1]
      bic_vec[j] <- (1-fit$dev.ratio)*fit$nulldev + log(N)*fit$df
      
      pb_ems$tick(1)
      
    }
    
    # -- for non-fixed item:
    for(j in (1:J)[-fixed]){
      
      for(iM in 1:NM){
        
        exclude <- which(Models[iM,]==0)
        if(length(exclude) == K){  # null model
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
    
    # ---- Replace current parameter ----
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
  
  # ---- load true model ----
  A_t <- c(1,1,1,0,0,0,0,0,0,
           0,0,0,1,1,1,0,0,0,
           0,0,0,0,0,0,1,1,1)
  A_t <- matrix(A_t, nrow=3, ncol=9, byrow=TRUE)
  fixed <- c(1,4,7)
  
  J <- ncol(A_t)  # no. of items
  K <- nrow(A_t)  # no. of latent traits
  N <- 1000       # no. of subjects
  
  # ---- true parameter setting ----
  b_t     <- rep(0,J)
  mu_t    <- rep(0,K)
  sigma_t <- matrix(.1,K,K); diag(sigma_t) <- 1
  
  # ---- generate random sample ----
  set.seed(4)
  x <- mvrnorm(n=N, mu=mu_t, Sigma=sigma_t)   # latent variable, here just only to generate y
  y <- x %>%
    `%*%` (A_t) %>%
    `+` (matrix(data=b_t,nrow=N,ncol=J,byrow=T)) %>%
    plogis(q=.) %>%
    rbinom(n=N*J, size=1, prob=.) %>%
    matrix(data=., nrow=N, ncol=J, byrow=F)
  
  # --- E-MS initial parameters & settings ----
  A_init <- matrix(data=1/J, nrow=K, ncol=J, byrow=TRUE); A_init[,fixed] <- A_t[,fixed]
  b_init <- rep(0,J)
  grid_num <- 11
  pts_pct  <- 0.05
  pts_cumu <- 0.90
  
  # ---- E-MS sp ----
  EMS_MIRT_sp(y=y,               # data set, responses of all subjects
              mu_t=mu_t,         # true mu (mean vector of latent traits)
              sigma_t=sigma_t,   # true sigma (covariance of latent traits)
              fixed=fixed,       # which item not do model selection
              A_init=A_init,     # initial value of A for EMS algorithm
              b_init=b_init,     # initial value of b for EMS algorithm
              grid_num=grid_num, # number of grid points for each latent variable
              pts_cumu=pts_cumu, # cumulative percent of p tilde for each subject
              pts_pct=pts_pct    # percentage of the least grid points used for calculation
              
  )
}
