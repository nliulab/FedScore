setwd('~/DukeNUS/fedscoreH')
repo <- '~/DukeNUS/fedscoreH/'

dyn.load(paste0(repo,"C/mylik_odal.so"))  
dyn.load(paste0(repo,"C/mylik_gradient.so"))
dyn.load(paste0(repo,"C/mylik_gradient_second.so"))
dyn.load(paste0(repo,"C/cov_indep_odal.so"))
Sys.setenv(RGL_USE_NULL=TRUE)
library(rjson)
library(meta)
library(logistf)
library(matlib)
library(glmmML)
library(MASS)

var_func <- function(par, K, n, length_par, x_all, y) {
  #######
  # A
  #######
  dev_score_tmp <- array(dim = c(K,length_par^2))
  for (k in 1:K) {
    temp1 <- .C("cal_dev_score_all",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
               as.double(par),as.double(length_par),result=double(length_par^2))
    dev_score_tmp[k,] <- temp1[["result"]]
  }
  
  # diagnal of the A matrix
  dev_score_not_full <- array(dim = c(K,length_par^2))
  for (i in 1:length_par) {
    dev_score_not_full[,1 + (length_par+1)*(i-1)] <- dev_score_tmp[,i]
  }
  
  # funtion to get a symmatric matrix
  makeSymm <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
  
  # rearrange the matrix
  dev_score <- array(dim = c(K,length_par^2))
  for (row in 1:K) {
    tmp_matrix <- matrix(dev_score_not_full[row,], nrow = length_par, byrow = TRUE)
    tmp_matrix[upper.tri(tmp_matrix, diag=FALSE)] <- dev_score_tmp[row, ((length_par+1):((((length_par^2)+length_par))/2))]
    tmp_matrix2 = makeSymm(tmp_matrix)
    dev_score[row,] = as.vector(tmp_matrix2)
  }
  
  # A in sandwich
  A.tmp = apply(dev_score, 2, sum)
  A = matrix(A.tmp, nrow = length_par, byrow = TRUE)
  
  #######
  # B
  #######
  score = array(dim = c(K,length_par))
  for(k in 1:K){
    temp2<- .C("cal_score_all",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
               as.double(par),as.double(length_par),result=double(length_par))
    score[k,] = temp2[["result"]]
  }
  B = t(score) %*% score 
  
  # get the result
  tmp = solve(A) %*% B %*% solve(A)
  esvar = tmp
  
  return(esvar)
}

run_dCLR <- function(site1, site2, uni_cut, var_list) {
  K <- 2
  df_temp <- rbind(site1, site2)
  df2 = df_temp[,!names(df_temp) %in% c("label")]
  df3 = as.data.frame(scale(df2))

  df <- cbind(df_temp$label, df3)
  colnames(df)[colnames(df) == "df_temp$label"] <- "label"
  df <- df[, c(var_list, "label")]
  df$site <- c(rep(1, nrow(site1)), rep(2, nrow(site2)))
  n <- c(nrow(site1), nrow(site2))
  length_par <- length(var_list)
  outcm_name <- 'label'
  fit2 <- glmmboot(label ~ .-site, cluster = site, data = df)
  espar_esm_gold_1 =  c(fit2$coefficients)
  espar_esm_gold_var_1 = c(diag(fit2$variance))
  gold_MLE = espar_esm_gold_1
  print(gold_MLE)
  length_par <- length(gold_MLE)
  y <- matrix(NA, nrow = K, ncol = max(n))
  for (k in 1:K) {
    y[k, 1:n[k]] <- df[df$site == k, outcm_name]
  }
  
  x_tmp <- df[, -which(names(df) == outcm_name)]
  x_all <- matrix(NA, nrow = K, ncol = length(var_list) * max(n))
  for (k in 1:K) {
    x_all[k, 1:(n[k]*length(var_list))] <- unlist(x_tmp[x_tmp$site == k, -which(names(x_tmp) == "site")])
  }

  lik_local_list = matrix(NA, nrow = K, ncol = length_par)
  local_var_list = matrix(NA, nrow = K, ncol = length_par^2)
  for (local_num in 1:K) {
    # local likelihood
    lik_local <- function(par) {
      temp <- .C("mylik_all",as.integer(n[local_num]),as.double(y[local_num,][!is.na(y[local_num,])]),
                as.double(x_all[local_num,][!is.na(x_all[local_num,])]),
                as.double(par),as.double(length_par),result=double(1))
      each_n <- n[local_num]*(n[local_num]-1)/2
      lik <- temp[["result"]]/each_n
      return(lik)
    }
    # get the result for local sites
        espar_esm_local=NA
        op_local=optim(espar_esm_gold_1,
                       lik_local,
                       control = list(fnscale=-1,maxit=1000),
                       method = "Nelder-Mead")
        espar_esm_local=op_local$par
        lik_local_list[local_num,] = espar_esm_local
        local_var_list[local_num,] = var_func(espar_esm_local, K, n, length_par, x_all, y)
  }
  print(lik_local_list)
  print(local_var_list)

  #### weighting of the broadcast value
  est_tmp = rep(0, length_par)
  est_tmp_2 = rep(0, length_par^2)
  for (index in 1:K){
    tmp1 = local_var_list[index,]
    tmp2 = lik_local_list[index,]
    print(tmp1)
    print(tmp2)
    if (!any(is.na(tmp1)) & !any(is.na(tmp2))){
          est_tmp = est_tmp + t(Ginv(matrix(tmp1, length_par, length_par))%*% tmp2)
          est_tmp_2 = est_tmp_2 + as.vector(Ginv(matrix(tmp1, length_par, length_par)))
    } else {
      est_tmp = est_tmp
      est_tmp_2 = est_tmp_2
    }
  }
  esm_init_bc = c(t(Ginv(matrix(est_tmp_2,length_par,length_par)) %*%t(est_tmp)))
  print(esm_init_bc)

  ## calculate gradients
  ## with the beta-bar obtain the gradients from all the sites with the initial values
  # first order
  ourlik_all=function(par){
    grad=matrix(0,nrow = K,ncol=length_par)
    for(k in 1:K){
      temp<- .C("mylik_gradient",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
                as.double(par),as.double(length_par),result=double(length_par))
      each_n = n[k]*(n[k]-1)/2
      grad[k,]=temp[["result"]]/each_n
    }
    return(grad) # grad divided by the size of each site
  }
  # get the gradients
  grads = ourlik_all(esm_init_bc)
  mean_grads = apply(grads, 2, mean)

  print(mean_grads)

  # second order
  ourlik_all_second=function(par)
  {
    grad=matrix(0,nrow = K,ncol=length_par^2)
    for(k in 1:K){
      temp<- .C("mylik_gradient_second",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
                as.double(par),as.double(length_par), result=double(length_par^2))
      each_n = n[k]*(n[k]-1)/2
      grad[k,]=temp[["result"]]/each_n
    }
    return(grad) # grad divided by the size of each site
  }
  # get the gradients
  second_grads = ourlik_all_second(esm_init_bc)
  mean_2_grads = matrix(c(apply(second_grads, 2, mean)), length_par, length_par)
  
  print(mean_2_grads)

  ## construct surrogate likelihood function

  # every site its own data and the
  est_matrix_second = matrix(NA, nrow = K, ncol = length_par)
  first_grd_var = second_grd_var = matrix(NA, nrow = K, ncol = length_par^2)
  for (local_num in 1:K){
    print(local_num)
    # local likelihood
    lik_local <- function(par){
      temp<- .C("mylik_all",as.integer(n[local_num]),as.double(y[local_num,][!is.na(y[local_num,])]),as.double(x_all[local_num,][!is.na(x_all[local_num,])]),
                as.double(par),as.double(length_par),result=double(1))
      each_n = n[local_num]*(n[local_num]-1)/2
      lik=temp[["result"]]/each_n
      return(lik)
    }

    ############ step 4: #########################
    ## construct surrogate likelihood function
    # l̃1(β)=l1(β)+{∇l(β¯)−∇l1(β¯)}(β−β¯)+(β−β¯)T{∇2l(β¯)−∇2l1(β¯)}(β−β¯)/2,
    surr_lik_second <- function(par) {
      lik_local(par) +
        (mean_grads - grads[1,]) %*% (par - esm_init_bc) +
        (0.5*(t(par - esm_init_bc)
              %*% (mean_2_grads - matrix(second_grads[1,], length_par, length_par))
              %*% (par - esm_init_bc)))
    }

        esm_final_tmp_second = NA
        op_final_tmp_second=optim(esm_init_bc,
                                  surr_lik_second,
                                  control = list(fnscale = -1,maxit = 1000),
                                  method = "Nelder-Mead")
        esm_final_tmp_second = op_final_tmp_second$par

        # assign the final answer
        est_matrix_second[local_num, ] = esm_final_tmp_second
        second_grd_var[local_num,] = c(var_func(esm_final_tmp_second, K, n, length_par, x_all, y))

  }

  ## Final Esitimation of dCLR(distributed Logistic Regression)
  ### A. Each site has its own surrogate pairwise likelihood estimate
  dclr_est_by_site = data.frame(est_matrix_second,row.names = paste0("dCLR_by_site",1:K))
  dclr_est_by_site


  ### B. An overall weighted estimate
  #### weighting of the final value
  # weighted average ### first gradient
  #### weighting of the broadcast value
  # weighted average ### second gradient
  est_tmp = rep(0, length_par)
  est_tmp_2 = rep(0, length_par^2)
  for (index in 1:K){
    tmp1 = second_grd_var[index,]
    tmp2 = est_matrix_second[index,]
    if (!any(is.na(tmp1)) & !any(is.na(tmp2))){
      est_tmp = est_tmp + t(ginv(matrix(tmp1, length_par, length_par) + diag(1e-6, length_par)) %*% est_matrix_second[index,])
      est_tmp_2 = est_tmp_2 + as.vector(ginv(matrix(tmp1, length_par, length_par) + diag(1e-6, length_par)))
    } else {
      est_tmp = est_tmp
      est_tmp_2 = est_tmp_2
    }
  }
  print(est_tmp_2)

  esm_final_second = c(t(solve(matrix(est_tmp_2,length_par,length_par)) %*% t(est_tmp)))
  esm_final_second_var = var_func(esm_final_second, K, n, length_par, x_all, y)

  dclr_est_weighted = esm_final_second
  print(dclr_est_weighted)

  return(dclr_est_weighted)
}

run_dCLR_cat <- function(site1, site2, uni_cut, var_list) {
  K <- 2
  df <- rbind(site1, site2)
  df <- df[, c(var_list, "label")]
  df$site <- c(rep(1, nrow(site1)), rep(2, nrow(site2)))
  n <- c(nrow(site1), nrow(site2))
  df[var_list] <- lapply(df[var_list], as.factor)
  print(df)

  df_encoded <- model.matrix(~ ., data = df[, c(var_list, "site")])
  
  print(colnames(df_encoded))
  length_par <- ncol(df_encoded) - 2
  print(length_par)
  outcm_name <- 'label'
  fit2 <- glmmboot(label ~ .-site, cluster = site, data = df)
  espar_esm_gold_1 =  c(fit2$coefficients)
  espar_esm_gold_var_1 = c(diag(fit2$variance))
  gold_MLE = espar_esm_gold_1
  print(gold_MLE)
  
  y <- matrix(NA, nrow = K, ncol = max(n))
  for (k in 1:K) {
    y[k, 1:n[k]] <- df[df$site == k, outcm_name]
  }
  
  x_tmp <- df_encoded[, -which(colnames(df_encoded) %in% c("site", "(Intercept)"))]  # Remove site column
  x_all <- matrix(NA, nrow = K, ncol = length_par * max(n))
  for (k in 1:K) {
    x_all[k, 1:(n[k] * length_par)] <- unlist(x_tmp[df$site == k, ])
  }
  lik_local_list = matrix(NA, nrow = K, ncol = length_par)
  local_var_list = matrix(NA, nrow = K, ncol = length_par^2)
  print(dim(lik_local_list))
  
  for (local_num in 1:K) {
    print(local_num)
    # local likelihood
    lik_local <- function(par) {
      temp <- .C("mylik_all",as.integer(n[local_num]),as.double(y[local_num,][!is.na(y[local_num,])]),
                 as.double(x_all[local_num,][!is.na(x_all[local_num,])]),
                 as.double(par),as.double(length_par),result=double(1))
      each_n <- n[local_num]*(n[local_num]-1)/2
      lik <- temp[["result"]]/each_n
      return(lik)
    }
    # get the result for local sites
    espar_esm_local=NA
    op_local=optim(espar_esm_gold_1,
                   lik_local,
                   control = list(fnscale=-1,maxit=1000),
                   method = "Nelder-Mead")
    espar_esm_local <- op_local$par
    print(espar_esm_local)
    print(var_func(espar_esm_local, K, n, length_par, x_all, y))
    lik_local_list[local_num,] = espar_esm_local
    local_var_list[local_num,] = var_func(espar_esm_local, K, n, length_par, x_all, y)
  }
  print(lik_local_list)
  print(local_var_list)

  # #### weighting of the broadcast value
  est_tmp = rep(0, length_par - 1)
  est_tmp_2 = rep(0, length_par^2)
  for (index in 1:K){
    tmp1 = local_var_list[index,]
    tmp2 = lik_local_list[index,]
    print(tmp1)
    print(tmp2)
    if (!any(is.na(tmp1)) & !any(is.na(tmp2))){
      est_tmp = est_tmp + t(Ginv(matrix(tmp1, length_par, length_par))%*% tmp2)
      est_tmp_2 = est_tmp_2 + as.vector(Ginv(matrix(tmp1, length_par, length_par)))
    } else {
      est_tmp = est_tmp
      est_tmp_2 = est_tmp_2
    }
  }
  esm_init_bc = c(t(Ginv(matrix(est_tmp_2,length_par,length_par)) %*%t(est_tmp)))
  print(esm_init_bc)

  ## calculate gradients
  ## with the beta-bar obtain the gradients from all the sites with the initial values
  # first order
  ourlik_all=function(par){
    grad <- matrix(0, nrow = K, ncol = length_par)
    for (k in 1:K) {
      temp <- .C("mylik_gradient",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
                as.double(par),as.double(length_par),result=double(length_par))
      each_n = n[k]*(n[k]-1)/2
      grad[k,]=temp[["result"]]/each_n
    }
    return(grad) # grad divided by the size of each site
  }
  # get the gradients
  grads = ourlik_all(esm_init_bc)
  mean_grads = apply(grads, 2, mean)

  print(mean_grads)

  # second order
  ourlik_all_second=function(par)
  {
    grad=matrix(0,nrow = K,ncol=length_par^2)
    for(k in 1:K){
      temp<- .C("mylik_gradient_second",as.integer(n[k]),as.double(y[k,][!is.na(y[k,])]),as.double(x_all[k,][!is.na(x_all[k,])]),
                as.double(par),as.double(length_par), result=double(length_par^2))
      each_n = n[k]*(n[k]-1)/2
      grad[k,]=temp[["result"]]/each_n
    }
    return(grad) # grad divided by the size of each site
  }
  # get the gradients
  second_grads = ourlik_all_second(esm_init_bc)
  mean_2_grads = matrix(c(apply(second_grads, 2, mean)), length_par, length_par)

  print(mean_2_grads)

  ## construct surrogate likelihood function

  # every site its own data and the
  est_matrix_second = matrix(NA, nrow = K, ncol = length_par)
  first_grd_var = second_grd_var = matrix(NA, nrow = K, ncol = length_par^2)
  for (local_num in 1:K){
    print(local_num)
    # local likelihood
    lik_local <- function(par){
      temp<- .C("mylik_all",as.integer(n[local_num]),as.double(y[local_num,][!is.na(y[local_num,])]),as.double(x_all[local_num,][!is.na(x_all[local_num,])]),
                as.double(par),as.double(length_par),result=double(1))
      each_n = n[local_num]*(n[local_num]-1)/2
      lik=temp[["result"]]/each_n
      return(lik)
    }

    ############ step 4: #########################
    ## construct surrogate likelihood function
    # l̃1(β)=l1(β)+{∇l(β¯)−∇l1(β¯)}(β−β¯)+(β−β¯)T{∇2l(β¯)−∇2l1(β¯)}(β−β¯)/2,
    surr_lik_second <- function(par) {
      lik_local(par) +
        (mean_grads - grads[1,]) %*% (par - esm_init_bc) +
        (0.5*(t(par - esm_init_bc)
              %*% (mean_2_grads - matrix(second_grads[1,], length_par, length_par))
              %*% (par - esm_init_bc)))
    }

    esm_final_tmp_second = NA
    op_final_tmp_second=optim(esm_init_bc,
                              surr_lik_second,
                              control = list(fnscale = -1,maxit = 1000),
                              method = "Nelder-Mead")
    esm_final_tmp_second = op_final_tmp_second$par

    # assign the final answer
    est_matrix_second[local_num, ] = esm_final_tmp_second
    second_grd_var[local_num,] = c(var_func(esm_final_tmp_second, K, n, length_par, x_all, y))

  }
   
  ## Final Esitimation of dCLR(distributed Logistic Regression)
  ### A. Each site has its own surrogate pairwise likelihood estimate
  dclr_est_by_site = data.frame(est_matrix_second,row.names = paste0("dCLR_by_site",1:K))
  dclr_est_by_site


  ### B. An overall weighted estimate
  #### weighting of the final value
  # weighted average ### first gradient
  #### weighting of the broadcast value
  # weighted average ### second gradient
  est_tmp = rep(0, length_par)
  est_tmp_2 = rep(0, length_par^2)
  for (index in 1:K){
    tmp1 = second_grd_var[index,]
    tmp2 = est_matrix_second[index,]
    if (!any(is.na(tmp1)) & !any(is.na(tmp2))){
      est_tmp = est_tmp + t(ginv(matrix(tmp1, length_par, length_par) + diag(1e-6, length_par)) %*% est_matrix_second[index,])
      est_tmp_2 = est_tmp_2 + as.vector(ginv(matrix(tmp1, length_par, length_par) + diag(1e-6, length_par)))
    } else {
      est_tmp = est_tmp
      est_tmp_2 = est_tmp_2
    }
  }
  print(est_tmp_2)

  esm_final_second = c(t(solve(matrix(est_tmp_2,length_par,length_par)) %*% t(est_tmp)))
  esm_final_second_var = var_func(esm_final_second, K, n, length_par, x_all, y)

  dclr_est_weighted = esm_final_second
  print(dclr_est_weighted)

  return(dclr_est_weighted)
}
