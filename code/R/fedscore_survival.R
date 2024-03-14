library(mle.tools)
library(rjson)
library(survival)

# ODAC for homogeneous survival data
# ODACH for heterogeneous survival data

############################################################################################################
############################################################################################################
############################################################################################################
# Algorithms Based on
# Duan, Rui, et al. "Learning from local to global: An efficient distributed algorithm for modeling time-to-event data." Journal of the American Medical Informatics Association 27.7 (2020): 1028-1036.
# &
# Luo, Chongliang, et al. "ODACH: a one-shot distributed algorithm for Cox model with heterogeneous multi-center data." Scientific reports 12.1 (2022): 6627.

# code adapted from https://github.com/Penncil/pda

# Steps (ODAC):
# 1. initial \beta and \V from local site (ODAC & ODACH)

# 2. initial communication: transfer result in 1 to other sites (ODAC & ODACH)

# 3.a  for other sites (ODAC)
#     1) required: U, W, Z
#     1) compute first and second derivative
#     2) transfer 1) to local site

# 3.b  for other sites (ODACH)
#     1) required: U, W, Z
#     2) compute first and second derivative
#     3) transfer 1) to local site

# 4. get surrogate likelihood

# 5. get estimators

# 6. get variance of the two estimators


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# Basic functions for ODAL, ODAL2, robust ODAL


##### Get beta (optimal) estimates of Cox regression
##### input:
##### data_matrix, matrix for regression
##### Note that the response in the data_matrix should be change to "y" for simplicity and standardization
get_beta <- function(data_matrix) {
  model <- coxph(Surv(label_time, label_status) ~ ., data = data_matrix)
  beta <- model$coefficients
  return(beta)
}

##### Get variance estimates of Cox regression
##### input:
##### data_matrix, matrix for regression
##### Note that the response in the data_matrix should be change to "y" for simplicity and standardization
get_V <- function(data_matrix) {
  model <- coxph(Surv(label_time, label_status) ~ ., data = data_matrix)
  Var <- model$var
  return(diag(Var))
}

# Calculate the initial value beta bar for surrogate likelihood function
get_bbar <- function(modelIndex, site_list, p, i_fold) {
  bhat_wt_sum <- rep(0, p)
  wt_sum <- rep(0, p)
  for (i in length(site_list)) {
    bhat <- read_beta(site_list[i], modelIndex, i_fold) %>% unlist
    Vhat <- read_Sigma(site_list[i], modelIndex, i_fold) %>% unlist
    bhat_wt_sum <- bhat_wt_sum + bhat/Vhat
    wt_sum <- 1 / Vhat
  }
  T_all <- get_unique_time(site_list, modelIndex, i_fold)
  nt <- length(T_all)
  bbar <- bhat_wt_sum/wt_sum # initial value of beta
  res <- list(T_all = T_all, nt = nt, bbar = bbar)
  return(res)
}

# Get unique time points for the given vector of site names
get_unique_time <- function(site_list, modelIndex, i_fold) {
  t_all <- c()
  for (i in 1:length(site_list)) {
    t_temp <- read_tunique(site_list[i], modelIndex, i_fold)
    t_all <- c(t_all, t_temp)
  }
  return(sort(unique(t_all)))
}

# Calculate U, W, Z for and write into intermediate.json
# Notes: It is common in Cox regression to add fake data points at each time point to calculate the summary statistics in the risk sets. 
# In this code, pfdata is created to hold these fake data points.
get_UWZ <- function(data_matrix, modelIndex, site_nameR, site_list, i_fold) {
  # get design matrix first:
  design <- get_designX(data_matrix)
  dat <- cbind(data_matrix[, c("label_time", "label_status")], design)
  p <- ncol(design)
  # get results from site_beta.json:
  res <- get_bbar(modelIndex, site_list, p, i_fold)
  bbar <- res$bbar
  T_all <- res$T_all
  nt <- res$nt
  t_max <- max(res$T_all) + 1
  # reform matrix:
  pfdata <- cbind(T_all, 0, matrix(0, nt, p)) %>% as.data.frame()
  colnames(pfdata) <- colnames(dat)
  pfdata <- rbind(dat, pfdata) %>% as.data.table() 
  pfdata <- pfdata[, 'interval':=cut(pfdata$label_time, breaks = c(T_all, t_max), labels = 1:nt, right = F)][order(pfdata$interval),]
  pfdata$interval[is.na(pfdata$interval)] <- nt
  X <- as.matrix(pfdata[, -c(1, 2)][,-'interval'])
  # Get summary statistics U, Z and W
  eXb <- c(exp(X %*% bbar))
  X2 <- X[, 1] * X
  for (i in 2:ncol(X)) {
    X2 <- cbind(X2, X[,i] * X)
  }
  UWZ <- eXb * cbind(1, X, X2)
  # functions for calculating column-wise (reverse) cumsum
  # see details in 'pda' package
  UWZ <- rcpp_aggregate(x = UWZ, indices = pfdata$interval, cumulative = T, reversely = T)
  # functions for calculating column-wise (reverse) cumsum
  # see details in 'pda' package and helpers.R
  U <- UWZ[,1] - c(nt:1)
  W <- UWZ[, 2:(p + 1)]
  Z <- array(UWZ[, -c(1:(p + 1))], c(nt, p, p))
  return(list(U = U, W = W, Z = Z, n = nrow(data_matrix)))
}


##### Obtain derivatives (by ODAC) at remote site and write into intermediate.json
##### input: 
##### data_matrix: data matrix at remote site
get_dev_ODAC <- function(data_matrix, modelIndex, site_nameR, site_list, i_fold) {
  U <- 0
  W <- 0
  Z <- 0
  # Sum up U, W, Z:
  for (i in 1:length(site_list)) {
    dev <- get_UWZ(data_matrix, modelIndex, site_nameR, site_list, i_fold)
    U <- U + dev$U
    W <- W + dev$W
    Z <- Z + dev$Z
  }
  
  # get design matrix
  design <- get_designX(data_matrix)
  data_matrix <- cbind(data_matrix[, c("label_time", "label_status")], design)
  p <- ncol(data_matrix) - 2
  
  bbar <- get_bbar(modelIndex, site_list, p, i_fold)$bbar
  T_all <- get_unique_time(site_list, modelIndex, i_fold)
  
  d <- c(table(c(data_matrix[data_matrix$label_status == T, "label_time"], T_all)) - 1)
  
  X <- as.matrix(data_matrix[data_matrix$label_status == T, -c(1, 2)])
  
  logL_D1 <- apply(X, 2, sum) - apply(d * W / U, 2, sum, na.rm = T)
  W2 <- array(NA, c(dim(W), p))  
  
  for (i in 1:p) {
    W2[,,i] <- W[,i] * W
  }
  logL_D2 <- apply(d * (W2 - U * Z) / U^2, c(2, 3), sum, na.rm = T)
  
  derivatives <- list(bbar = bbar,
                     T_all = T_all,
                     U = U,
                     Z = Z,
                     site = site_nameR,
                     site_size = nrow(data_matrix),
                     logL_D1 = logL_D1,
                     logL_D2 = logL_D2)
  return(derivatives)
}


##### obtain derivatives (by ODACH) at remote site, using beta from local site
##### input: 
##### data_matrix: data matrix at remote site
get_dev_ODACH <- function(data_matrix, modelIndex, site_nameR, site_list, i_fold) {
  dev <- get_UWZ(data_matrix,  modelIndex, site_nameR, site_list, i_fold)
  time <- data_matrix$label_time
  status <- data_matrix$label_status
  n <- length(time)
  
  # get design matix
  design <- get_designX(data_matrix)
  p <- ncol(design)
  data_matrix <- cbind(data_matrix[, c("label_time", "label_status")], design)
  X <- design
  bbar <- get_bbar(modelIndex, site_list, p, i_fold)$bbar
  
  
  hasTies <- any(duplicated(data_matrix$time))
  logL_D2 <- -matrix(rcpp_coxph_logL_hessian(bbar, time = time, event = status, z = X), p, p)
  
  if (hasTies) {
    logL_D1 <- -rcpp_coxph_logL_gradient_efron(bbar, time = time, event = status, z = X)
  } else {
    logL_D1 <- -rcpp_coxph_logL_gradient(bbar, time = time, event = status, z = X)
  }
  
  derivatives <- list(bbar = bbar, site = site_nameR, n = nrow(data_matrix), logL_D1 = logL_D1, logL_D2 = logL_D2)
  return(derivatives)
}

# get single surrogate likelihood estimates for one site
get_btilde <- function(data_matrix, sitename, site_list, modelIndex, heterogeneity, i_fold) {
  
  time <- data_matrix$label_time
  status <- data_matrix$label_status
  n <- length(time)
  X <- get_designX(data_matrix)
  n <- length(time)
  p <- ncol(X)
  data_matrix <- cbind(data_matrix[, c("label_time", "label_status")], X)
  hasTies <- any(duplicated(data_matrix$label_time))
  
  D1_all <- rep(0, p)
  D2_all <- matrix(0, p, p)
  N <- 0
  
  for (i in 1:length(site_list)) {
    Dev <- read_intermediate(sitename, modelIndex, heterogeneity, i_fold)
    D1_all <- D1_all + Dev$D1
    D2_all <- D2_all + Dev$D2
    N <- N + Dev$N_R
  }
  bbar <- get_bbar(modelIndex, site_list, p, i_fold)$bbar

  if (hasTies) {
    logL_local <- function(beta) -rcpp_coxph_logL_efron(beta, time, status, X)
    logL_local_D1 <- function(beta) -rcpp_coxph_logL_gradient_efron(beta, time, status, X)
    logL_local_D2 <- function(beta) -matrix(-rcpp_coxph_logL_hessian(beta, time, status, X), p, p)
  } else {
    logL_local <- function(beta) -rcpp_coxph_logL(beta, time, status, X)
    logL_local_D1 <- function(beta) -rcpp_coxph_logL_gradient(beta, time, status, X)
    logL_local_D2 <- function(beta) -matrix(-rcpp_coxph_logL_hessian(beta, time, status, X), p, p)
  }
  
  logL_diff_D1 <- D1_all/N - logL_local_D1(bbar)/n
  logL_diff_D2 <- D2_all/N - logL_local_D2(bbar)/n
  logL_tilde <- function(b) - (logL_local(b)/n + sum(b*logL_diff_D1) + 0.5 * t(b-bbar) %*% logL_diff_D2 %*% (b-bbar))
  
  print("Start to solve optim")
  startTime <- Sys.time()
  result_opt <- optim(par = bbar,
                     fn = logL_tilde,
                     hessian = T,
                     control = list(trace = T, reltol = 1e-3))
  endTime <- Sys.time()
  print(endTime - startTime)
  if (result_opt$convergence != 0) {
    print(result_opt$convergence)
    print("Warning! Not successful completion")
  } else {
    print("Successfuly finished optimization of surrogate likelihood")
  }
  result <- list(converge = result_opt$convergence, 
                 beta_global = result_opt$par, 
                 H = result_opt$hessian)
  return(result)
}

# Get final estimated ~beta by weighting local estimated ~beta
get_federatedCox_i <- function(data_matrix, site_list, modelIndex, heterogeneity, i_fold) {
  p <- ncol(get_designX(data_matrix)) 
  names <- colnames(get_designX(data_matrix))
  K <- length(site_list)
  btilde_wt_sum <- rep(0, p)
  wt_sum <- rep(0, p)
  
  for (i in 1:length(site_list)) {
    surr <- read_tbeta(site_list[i], modelIndex, heterogeneity, i_fold)
    btilde_wt_sum <- btilde_wt_sum %>% as.vector
    btilde_wt_sum <- btilde_wt_sum + matrix(surr$hessian, nrow = p, ncol = p) %*% unlist(surr$beta_global)
    wt_sum <- wt_sum + matrix(surr$hessian, nrow = p, ncol = p)
  }
  # inv_var weighted estimation
  return(list(beta_global = solve(wt_sum, btilde_wt_sum),
              V_global= solve(wt_sum) * K,
              coef_names = names))
}

# Get final estimated ~beta by weighting local estimated ~beta
get_federatedCox_final <- function(data_matrix, site_list, var, modelIndex, heterogeneity, i_fold = "NA", csv_name) {
  p <- ncol(get_designX(data_matrix)) 
  names <- colnames(get_designX(data_matrix))
  K <- length(site_list)
  btilde_wt_sum <- rep(0, p)
  wt_sum <- rep(0, p)
  
  for (i in 1:length(site_list)) {
    surr <- read_tbeta(site_list[i], modelIndex, heterogeneity, i_fold)
    btilde_wt_sum <- btilde_wt_sum %>% as.vector
    btilde_wt_sum <- btilde_wt_sum + matrix(surr$hessian, nrow = p, ncol = p) %*% unlist(surr$beta_global)
    wt_sum <- wt_sum + matrix(surr$hessian, nrow = p, ncol = p)
  }
  # inv_var weighted estimation
  coef <- solve(wt_sum, btilde_wt_sum) 
  coef <-  as.data.frame(coef)
  rownames(coef) <- names
  file_name <- paste(csv_name, sprintf("var%d_model%s_fold%s.csv", var, modelIndex, i_fold),sep = "")
  write.csv(coef, file_name)
}


get_federatedCox <- function(data_train_all, site_list, var_list, heterogeneity, uni_cut, fold = "no_cv", 
                             csv_name = "../../output/survival/fedsurv/coef_global_homoS1_trainCV") {
  if (fold != "no_cv") {
    for (fold in 1:fold) {
      cat("starting fold:")
      cat(fold, "\n")
      data_train <- get_fold(data_train_all, fold)
      data_train <- transform_df_fixed(data_train, uni_cut)
      
      for (model_index in 1:length(var_list)) {
        cat("starting model:")
        cat(model_index, "\n")
        var_m <- var_list[1:model_index]
        data_train_temp <- data_train[, c("label_status", "label_time", var_m)]
        model <- get_federatedCox_i(data_train_temp, site_list, model_index, heterogeneity, fold)
        coef <- model$beta_global %>% as.data.frame()
        rownames(coef) <- model$coef_names
        file_name <- paste(csv_name, sprintf("_model%d_fold%s.csv", model_index, fold), sep = "")
        write.csv(coef, file_name)
      }
    }
  } else {
    data_train <- transform_df_fixed(data_train_all, uni_cut)
    for (model_index in 1:length(var_list)) {
      cat("starting model:")
      cat(model_index, "\n")
      var_m <- var_list[1:model_index]
      data_train_temp <- data_train[, c("label_status", "label_time", var_m)]
      model <- get_federatedCox_i(data_train_temp, site_list, model_index, heterogeneity, fold)
      coef <- model$beta_global %>% as.data.frame()
      rownames(coef) <- model$coef_names
      file_name <- paste(csv_name, sprintf("_model%d_fold%s.csv", model_index, 1), sep="")
      write.csv(coef, file_name)
    }
  }
}

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# User level function

##### function for remote site to calculate intermediate values and write them to intermediate.json
##### file: ranking_site_beta.json, ranking_intermediate.json
##### input:
##### local and remote sitename
##### data_matrix, remote site data, col must be variable ranked by importance.
write_Ds_ranking <- function(sitenameL, sitenameR, data_matrix) {
  COl <- colnames(data_matrix)
  V2 <- COL[-1]
  list_read <- read_site_beta_ranking (sitename = sitenameL)
  # get info from ranking_site_beta.json
  Lbeta <- list_read$beta_local
  V1 <- list_read$variable
  # calculate D
  Dev <- get_dev(Lbeta, data_matrix)
  # check if models are correctly matched:
  if (V1 != V2) {
    print("Models don't match! Recheck if input matrix is formatted correctly.")
    return()
  }
  # write in intermediate.json
  write_intermediate_ranking(sitenameL = sitenameL, sitenameR = sitenameR, Lbeta = Lbeta, 
                             vairableList = V1, D1 = Dev$logL_D1, D2 = Dev$logL_D2)
}




######################################################
# get federated ranking
# weights: either default setting using 1/K or user specified standardized weights (sum of weights equals to 1)
# path: path of the "rank" folder
AutoScore_Frank <- function(weight = "default", K = 2, path = "rank/") {
  files <- list.files(path = path, pattern = "*.csv")
  rank <- read.csv(paste(path, files[1], sep = ""))
  char <- gsub("\\..*", "", files[1])
  colnames(rank) <- c("var", paste("imp", char, sep = ""), char)
  tbl <- rank
  for (i in 2:length(files)) {
    rank <- read.csv(paste(path, files[i], sep = ""))
    char <- gsub("\\..*","", files[i])
    colnames(rank) <- c("var", paste("imp", char,sep = ""), char)
    tbl <- left_join(tbl, rank)
  }
  dat <- select(tbl, matches("^rank")) %>% as.matrix
  if (typeof(weight) == "character") {
    rank_fed <- dat %*% rep(1/K, K)
  } else {
    if (sum(weight) != 1) {
      print("Error! Input for 'weight' should be either a vector of standardized weights or use default option")
    }
    rank_fed <- dat %*% weight
  }
  tbl$rank_fed <- rank_fed
  tbl <- tbl %>% arrange(rank_fed)
  return(tbl)
}
