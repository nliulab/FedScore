library(mle.tools)
library(rjson)


# ODAL2 for homogeneous data

############################################################################################################
############################################################################################################
############################################################################################################
# Algorithms Based on Rui Duan, et al. "Learning from electronic health records across multiple sites: 
#                 A communication-efficient and privacy-preserving distributed algorithm". 
#'                Journal of the American Medical Informatics Association, 2020, https://doi.org/10.1093/jamia/ocz199
# 


# Steps (ODAL2):
# 1. initial value obtained from local site

# 2. initial communication: transfer result in 1 to other sites

# 3.  for other sites
#     1) compute first and second derivative
#     2) transfer 1) to local site

# 4. get surrogate likelihood

# 5. get ODAL estimators

# 6. get variance of the two estimators

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
#Standardize data input


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# functions writing and reading json file
# each time write / read for one model only

##### write cutoff.json
##### input: 
##### sitename
##### train data
write_cutoff <- function(df, quantiles = c(0, 0.05, 0.2, 0.8, 0.95, 1), max_cluster = 5, categorize = "quantile", sitename, path = "cutoff.json"){
  cutoff <- get_cut_vec(df, quantiles, max_cluster, categorize)
  cutoff_new <- list("site" = sitename, "cutoff" = cutoff)
  all_cutoff <- rjson::fromJSON(file = path)
  # print(all_cutoff)
  
  K = length(all_cutoff)
  found = FALSE
  if(K>0){
  for (i in 1:K) {
    if (all_cutoff[[i]]$site == sitename) {
      all_cutoff[[i]] = cutoff_new
      found = TRUE
    }
  }
  }
  if (!found) {
    all_cutoff[[K+1]] = cutoff_new
  }
  json_data = rjson::toJSON(all_cutoff, 2)
  #print(json_data)
  write(json_data, path)
}


##### read cutoff.json
##### input:
##### sitename
read_cutoff <- function(sitename,path = "cutoff.json") {
  all_cutoff <- rjson::fromJSON(file = path)
  K = length(all_cutoff)
  for (i in 1:K) {
    if (all_cutoff[[i]]$site == sitename) {
      return(all_cutoff[[i]]$cutoff)
    }
  }
  return()
}


##### write score.json
##### input:
##### sitename
# write score table obtained based on each training data
# \sitename is the site name where  training data used comes from
# write_score <- function(score, sitename){
#   all_score <- rjson::fromJSON(file = "score.json")
#   score_new <- list(site = sitename, score_list = score)
#   
#   K = length(all_score)
#   found = FALSE
#   if(K>0){
#     for (i in 1:K) {
#       print(sitename)
#       print(all_score[[i]])
#       if (all_score[[i]]$site == sitename) {
#         all_score[[i]] = score_new
#         found = TRUE
#       }
#     }
#   }
#   if (!found) {
#     all_score[[K+1]] = score_new
#   }
#   json_data = rjson::toJSON(all_score, 2)
#   #print(json_data)
#   write(json_data, "score.json")
# }

##### read score.json
##### input:
##### sitename
read_score <- function(sitename){
  all_score <- rjson::fromJSON(file = "score.json")
  K = length(all_score)
  if(K>0){
    for (i in 1:K) {
      if (all_score[[i]]$site == sitename) {
        return(all_score[[i]]$score_list)
      }
    }
    return()
  }
}



##### write site_beta.json
##### input: 
##### sitename
##### modelIndex, range from 1 to 20
##### L_beta, optim beta calculated when choosing the site as local site
##### variableList: predictors included in the model

write_site_beta <- function(sitename, modelIndex, Lbeta, variableList) {
  mylist = list(site = sitename, model_index = modelIndex, beta_local = Lbeta, variable = variableList)
  all_site_beta <- rjson::fromJSON(file = "site_beta.json")
  K = length(all_site_beta)
  found = FALSE
  for (i in 1:K) {
    if (all_site_beta[[i]]$site == sitename & 
        all_site_beta[[i]]$model_index == modelIndex) {
      all_site_beta[[i]] = mylist
      found = TRUE
    }
  }
  if (!found) {
    all_site_beta[[K+1]] = mylist
  }
  json_data = rjson::toJSON(all_site_beta, 2)
  #print(json_data)
  write(json_data, "site_beta.json")
}

# test
# write_site_beta(sitename = "HA", modelIndex = 40, Lbeta = 233, vairableList = c("a", "b"))


##### read site_beta.json
##### input: sitename
##### return the optim beta by the site that is selected as local site
read_site_beta <- function(sitename = "SG", modelIndex) {
  all_site_beta <- rjson::fromJSON(file = "site_beta.json")
  K = length(all_site_beta)
  for (i in 1:K) {
    if (all_site_beta[[i]]$site == sitename & 
        all_site_beta[[i]]$model_index == modelIndex) {
      return(all_site_beta[[i]]$beta_local)
    }
  }
  return()
}

##### write intermediate.json
##### input: 
##### sitenameL: local site name
##### sitenameR: remote site name
##### modelIndex, range from 1 to 20
##### L_beta, optim beta calculated when choosing the site as local site
##### variableList: predictors included in the model
##### D1, D2: gradients calculated using local_site beta and remote site data
##### note that D2 need to be store by rows in json file.

write_intermediate <- function(sitenameL, sitenameR, modelIndex, Lbeta, variableList, D1, D2,n) {
  mylist = list(local_site = sitenameL, remote_site = sitenameR, model_index = modelIndex, beta_local = Lbeta, variable = variableList, D1= D1, D2 = D2, N_R = n)
  all_site_beta <- rjson::fromJSON(file = "intermediate.json")
  # print(all_site_beta)
  K = length(all_site_beta)
  found = FALSE
  for (i in 1:K) {
    if (all_site_beta[[i]]$local_site == sitenameL & 
        all_site_beta[[i]]$remote_site == sitenameR & 
        all_site_beta[[i]]$model_index == modelIndex) {
      all_site_beta[[i]] = mylist
      found = TRUE
    }
  }
  if (!found) {
    all_site_beta[[K+1]] = mylist
  }
  json_data = rjson::toJSON(all_site_beta, 2)
  #print(json_data)
  write(json_data, "intermediate.json")
}

# test
# write_intermediate (sitenameL = "C", sitenameR = "B",modelIndex = 2, Lbeta = 233, variableList = c("a", "b"),
#                     D1 = 0.5, D2 = matrix(5:8, nrow = 2, ncol = 2), n = 100)

##### read intermediate.json
##### input:
##### sitenameL: local site name
##### sitenameR: remote site name
##### modelIndex, range from 1 to 20
read_intermediate <- function(sitenameL, sitenameR, modelIndex) {
  all_site_beta <- rjson::fromJSON(file = "intermediate.json")
  K = length(all_site_beta)
  for (i in 1:K) {
    if (all_site_beta[[i]]$local_site == sitenameL & 
        all_site_beta[[i]]$remote_site == sitenameR & 
        all_site_beta[[i]]$model_index == modelIndex) {
      return(all_site_beta[[i]])
    }
  }
  return()
}


##### write ranking_site_beta.json
##### input: 
##### sitename
##### L_beta, optim beta calculated when choosing the site as local site
##### variableList: predictors included in the model

# write_site_beta_ranking <- function(sitename, Lbeta, vairableList) {
#   mylist = list(site = sitename,  beta_local = Lbeta, variable = vairableList)
#   all_site_beta <- rjson::fromJSON(file = "ranking_site_beta.json")
#   K = length(all_site_beta)
#   found = FALSE
#   for (i in 1:K) {
#     if (all_site_beta[[i]]$site == sitename) {
#       all_site_beta[[i]] = mylist
#       found = TRUE
#     }
#   }
#   if (!found) {
#     all_site_beta[[K+1]] = mylist
#   }
#   json_data = toJSON(all_site_beta, 2)
#   #print(json_data)
#   write(json_data, "ranking_site_beta.json")
# }

# test
# write_site_beta_ranking(sitename = "HA", Lbeta = 233, vairableList = c("a", "b"))


##### read ranking_site_beta.json
##### input: sitename
##### return the optim beta by the site that is selected as local site
# read_site_beta_ranking <- function(sitename = "SG") {
#   all_site_beta <- rjson::fromJSON(file = "site_beta_ranking.json")
#   K = length(all_site_beta)
#   for (i in 1:K) {
#     if (all_site_beta[[i]]$site == sitename) {
#       return(all_site_beta[[i]])
#     }
#   }
#   return()
# }

##### write ranking_intermediate.json
##### input: 
##### sitenameL: local site name
##### sitenameR: remote site name
##### modelIndex, range from 1 to 20
##### L_beta, optim beta calculated when choosing the site as local site
##### variableList: predictors included in the model
##### D1, D2: gradients calculated using local_site beta and remote site data
##### n: sample size of remote size
##### note that D2 need to be store by rows in json file.

# write_intermediate_ranking <- function(sitenameL, sitenameR, Lbeta, vairableList, D1, D2, n) {
#   mylist = list(local_site = sitenameL, remote_site = sitenameR,  beta_local = Lbeta, 
#                 variable = vairableList, D1= D1, D2 = 0, N_R = n)
#   mylist$D2 = list(D2[1,], D2[2,])
#   all_site_beta <- rjson::fromJSON(file = "ranking_intermediate.json")
#   K = length(all_site_beta)
#   found = FALSE
#   for (i in 1:K) {
#     if (all_site_beta[[i]]$local_site == sitenameL & 
#         all_site_beta[[i]]$remote_site == sitenameR) {
#       all_site_beta[[i]] = mylist
#       found = TRUE
#     }
#   }
#   if (!found) {
#     all_site_beta[[K+1]] = mylist
#   }
#   json_data = rjson::toJSON(all_site_beta, 2)
#   #print(json_data)
#   write(json_data, "ranking_intermediate.json")
# }

# test
# write_intermediate_ranking (sitenameL = "C", sitenameR = "B", Lbeta = 233, vairableList = c("a", "b"),
#                     D1 = 0.5, D2 = matrix(5:8, nrow = 2, ncol = 2), n = 233)

##### read intermediate.json
##### input:
##### sitenameL: local site name
##### sitenameR: remote site name
##### modelIndex, range from 1 to 20
# read_intermediate_ranking <- function(sitenameL, sitenameR) {
#   all_site_beta <- rjson::fromJSON(file = "ranking_intermediate.json")
#   K = length(all_site_beta)
#   for (i in 1:K) {
#     if (all_site_beta[[i]]$local_site == sitenameL & 
#         all_site_beta[[i]]$remote_site == sitenameR) {
#       return(all_site_beta[[i]]$site)
#     }
#   }
#   return()
# }


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# basic functions for ODAL, ODAL2, robust ODAL


##### Get beta (optimal) estimates of a logistic regression
##### input:
##### data_matrix, matrix for regression
##### Note that the response in the data_matrix should be change to "y" for simplicity and standardization
get_beta <- function(data_matrix){
  model <- glm(as.formula("label ~ ."), data=data_matrix,family="binomial", control = list(maxit = 50))
  beta <- model$coefficients
  return(beta)
}


##### Get log likelihood of a logistic regression
##### input: model_matrix, optimal beta
get_LLH <- function(beta, data_matrix){
  model <- glm(as.formula("label ~ ."),data=data_matrix,family="binomial")
  design = model.matrix(model)
  Y = data_matrix$label
  llh <- sum(Y*(as.matrix(design)%*%t(t(beta)))-log(1+exp(as.matrix(design)%*%t(t(beta)))))/length(Y)
  return(llh)
}


##### Sigmoid function
##### map input by sigmoid function
expit <- function(x){
  1/(1+exp(-x))
}


##### first order intermediate term
##### input:
##### beta: local estimation
##### data_matrix: data from remote site
##### return a vector of length p + 1
LLHD1 <- function(beta, data_matrix){
  model <- glm(as.formula("label ~ ."),data=data_matrix,family="binomial")
  design = model.matrix(model)
  Y = data_matrix$label 
  t(Y-expit(design%*%t(t(beta)))) %*% design / length(Y)
}


##### second order intermediate term
##### beta: local estimation
##### data_matrix: data from remote site
##### return a (p +1 ) by (p + 1) matrix
LLHD2 = function(beta,data_matrix){
  model <- glm(as.formula("label ~ ."),data=data_matrix,family="binomial")
  design = model.matrix(model)
  Z=expit(design%*%beta)
  t(c(-Z*(1-Z))*design)%*%design / nrow(design)
}


##### obtain derivatives at remote site, using beta from local site
##### input: 
##### beta_local: beta calculated at local site
##### data_matrix: data matrix at remote site
get_dev <- function(beta_local, data_matrix){
  # get first and second order derivatives and put in a list: 
  derivatives <- list(
    logL_D1 = LLHD1(beta_local,data_matrix),
    logL_D2 = LLHD2(beta_local,data_matrix)
  )
  
  return(derivatives)
}

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# User level function

##### function to input Lbeta calculated by each site
##### file: site_beta.json
##### input:
##### sitename
##### data_matrix, col must be variable ranked by importance.
write_Lbeta <- function(sitename, site_list, data_train, var, cutoff){
  N = length(var)    # total number of models
  # get unified cutoff calculated from cutoff.json
  # cutoff <- get_uni_cut(site_list)
  # transform cts variable to categorical:
  data_transformed = transform_df_fixed(data_train, cutoff)
  for(i in 0:(N-1)){
    k = i+1
    var_model = var[1:k]
    mat= data_transformed[,c("label",var_model)]
    beta = get_beta(mat)
    write_site_beta(sitename = sitename, modelIndex = k, Lbeta = beta, variableList = var_model)
    print(sprintf("Local beta for model %d has been written into file", k))
  }
}

# calculate parameters for the final selected model
write_Lbeta_final <- function(sitename, site_list, data_train, var, cutoff, modelIndex){
  N = length(var)
  data_transformed = transform_df_fixed(data_train, cutoff)
  beta = get_beta(data_transformed[,c("label",var)])
  write_site_beta(sitename = sitename, modelIndex = modelIndex, Lbeta = beta, variableList = var)
  print(sprintf("Local beta for model %s has been written into file", modelIndex))
}


##### function for remote site to calculate intermediate values and write them to intermediate.json
##### file: site_beta.json, intermediate.json
##### input:
##### sitename
##### data_matrix, remote site data, col must be variable ranked by importance.
write_Ds <- function(sitenameL, sitenameR, site_list, data_train, var, cutoff){
  N = length(var) # total number of models
  # get unified cutoff calculated from cutoff.json
  # cutoff <- get_uni_cut(site_list)
  # transform cts variable to categorical:
  data_transformed = transform_df_fixed(data_train, cutoff)
  for(i in 0:(N-1)){
    k = i+1
    var_model = var[1:k]
    mat= data_transformed[,c("label",var_model)]
    Lbeta  = read_site_beta(sitename = sitenameL, modelIndex = k) %>% unlist # get info from site_beta.json
    # calculate D
    Dev = get_dev(Lbeta, mat)
    D1 = as.vector(Dev$logL_D1)
    D2 = as.matrix(Dev$logL_D2)
    dimnames(D1) <- NULL
    dimnames(D2) <- NULL
    # write in intermediate.json
    write_intermediate(sitenameL = sitenameL, sitenameR = sitenameR, modelIndex = k, Lbeta = Lbeta, 
                       variableList = var_model, D1, D2, n = nrow(mat))
    print(sprintf("D1 and D2 for model %d has been written into file", k))
  }

}

# prepare parameters for final model
write_Ds_final <- function(sitenameL, site_list, data_train, var, cutoff, modelIndex){
  N = length(site_list) -1 
  site_other = site_list[site_list != sitenameL]
  data_transformed = transform_df_fixed(data_train, cutoff)
  mat= data_transformed[,c("label",var)]
  for(i in 1:N){
    sitenameR = site_other[i]
    Lbeta  = read_site_beta(sitename = sitenameL, modelIndex = modelIndex) %>% unlist # get info from site_beta.json
    print(Lbeta)
    print(mat)
    # calculate D
    Dev = get_dev(Lbeta, mat)
    D1 = as.vector(Dev$logL_D1)
    D2 = as.matrix(Dev$logL_D2)
    dimnames(D1) <- NULL
    dimnames(D2) <- NULL
    # write in intermediate.json
    write_intermediate(sitenameL = sitenameL, sitenameR = sitenameR, modelIndex = modelIndex, Lbeta = Lbeta, variableList = var, D1, D2, n = nrow(mat))
  }
  
}


##### function to input Lbeta calculated by each site
##### file: ranking_site_beta.json
##### input:
##### sitename
##### data_matrix, col must be variable ranked by importance.
write_Lbeta_ranking <- function(sitename, data_matrix){
  COl = colnames(data_matrix)
  COL = COL[-1]
  L_beta = get_beta(data_matrix[,k])
  write_site_beta_ranking(sitename = sitename, Lbeta = L_beta, vairableList = COL)
  print("Local beta for model %d has been written into file")
}


##### function for remote site to calculate intermediate values and write them to intermediate.json
##### file: ranking_site_beta.json, ranking_intermediate.json
##### input:
##### local and remote sitename
##### data_matrix, remote site data, col must be variable ranked by importance.
write_Ds_ranking <- function(sitenameL, sitenameR, data_matrix){
  COl = colnames(data_matrix)
  V2 = COL[-1]
  list_read = read_site_beta_ranking (sitename = sitenameL)
  # get info from ranking_site_beta.json
  Lbeta = list_read$beta_local
  V1 = list_read$variable
  # calculate D
  Dev = get_dev(Lbeta, data_matrix)
  # check if models are correctly matched:
  if(!(V1==V2)) {
      print("Models don't match! Recheck if input matrix is formatted correctly")
      return()
    }
  # write in intermediate.json
  write_intermediate_ranking(sitenameL = sitenameL, sitenameR = sitenameR, Lbeta = Lbeta, 
                       vairableList = V1, D1 = Dev$logL_D1, D2 = Dev$logL_D2)
}

##### function to build federated logistic regression
##### including estimated variance
##### return federated beta, and federated vcov
##### input: 
##### method: method for building surrogate likelihood
#####     1) ODAL2 (mean)
#####     2) ODAL  (mean)
#####     3) robust ODAL (median)
#####     currently only available for ODAL2
##### data_matrix: all data at local site
##### variableList: the vector of ranked variables (all)
##### return one global model
get_federatedLR <- function(method = 1, data_matrix, sitenameL, sitenameR_list, modelIndex){
  ####### get data matrix for the required model
  dat <- data_matrix
  
  ####### get info from local site
  beta_local <- get_beta(dat) # beta

  ####### get info from remote sites by reading ranking_intermediate.json
  size = length(sitenameR_list)
  derivatives_all <- vector(mode = "list", length = size)
  for(i in 1:size){
    list_read <- read_intermediate(sitenameL, sitenameR_list[i], modelIndex)
    
    # check if the model required matches with model stored in Json file
    # if(sum(var == unlist(list_read$variable))==modelIndex){
    #   print("Error! model doesn't match Json file! Recheck previous steps!")
    #   return(NULL)
    # }
    # store all first and second order D:
    derivatives_all[[i]]$logL_D1 <- list_read$D1
    derivatives_all[[i]]$logL_D2 <- list_read$D2
    derivatives_all[[i]]$n <- list_read$N_R
  }

  ####### second order surrogate likelihood
  # calculation for first order approx:
  D1_all <- 0
  # D1: vector of length p
  for(j in 1:size){
    D1_all <- D1_all + derivatives_all[[j]]$logL_D1 * derivatives_all[[j]]$n
  }
  ## calculation for second order approx:
  D2_all <- 0
  # D2: p x p matrix
  for(j in 1:size){
    D2_all <- D2_all + derivatives_all[[j]]$logL_D2 * derivatives_all[[j]]$n
  }
  
  N = 0
  for(i in 1:length(derivatives_all)){
    N = N + derivatives_all[[i]]$n
  }
  N = N + nrow(data_matrix) # number of total observations across all sites
  # print(D1_all/N)
  # print(length(D1_all/N))
  # print(beta_local)
  # print(length(beta_local))
  # print("Above is beta_local")
  # print(LLHD1(beta_local, dat))
  # print(length(LLHD1(beta_local, dat)))
  # print("602")
  # print(D1_all/N - LLHD1(beta_local, dat))
  # print("603")
  # print("Above is D1_all/N - LLHD1(beta_local, dat)")
  # print(matrix(D2_all/N, ncol = sqrt(length(D2_all))))
  # print("Above is D2_all/N in matrix format")
  # print(N)
  # print(LLHD2(beta_local, dat))
  # print("Above is LLHD2")
  ###### optimize the surrogate likelihood:
  surr = function(beta){
    -(get_LLH(beta, dat) +
        (D1_all/N - LLHD1(beta_local, dat)) %*% beta +
        t(beta-beta_local) %*% (matrix(D2_all/N, ncol = sqrt(length(D2_all))) - LLHD2(beta_local, dat)) %*% (beta-beta_local) / 2)
   }
  # general-purpose optimization, may need to check if it can always converge...
  # currently choose Newton method
  print("Start to solve optim")
  startTime <- Sys.time()
  # result_opt <- optim(par = beta_local,
  #                     fn = surr,
  #                     control = list(trace=T))
  # result_opt <- optim(par = beta_local,
  #                     fn = surr,
  #                     method = "CG",
  #                     )
  result_opt <- optim(par = beta_local,
                      fn = surr,
                      control = list(trace = T, reltol = 1e-3))
  endTime <- Sys.time()
  print(endTime - startTime)
  if(result_opt$convergence!=0){
    print(result_opt$convergence)
    print("Warning! Not successful completion")}
  else
  print("Successfuly finished optimization of surrogate likelihood")
  result <- list(converge = result_opt$convergence, beta_global = result_opt$par)
  return(result)
}



######################################################
# get federated ranking
# weights: either default setting using 1/K or user specified standardized weights (sum of weights eqauls to 1)
AutoScore_Frank <- function(weight = "default", K = 2){
  # loop through files in the "outputs/ranks/" folder:
  files = list.files(path = "outputs/ranks/", pattern="*.csv")
  rank = read.csv(paste("outputs/ranks/", files[1], sep = ""))
  char = gsub("\\..*","",files[1])
  colnames(rank) = c("var", paste("imp", char,sep = ""), char)
  tbl = rank
  for(i in 2:length(files)){
    rank = read.csv(paste("outputs/ranks/", files[i], sep = ""))
    char = gsub("\\..*","",files[i])
    colnames(rank) = c("var", paste("imp", char,sep = ""), char)
    tbl = left_join(tbl, rank)
  }
  dat = select(tbl, matches("^rank")) %>% as.matrix
  if(typeof(weight) == "character"){
    rank_fed = dat %*% rep(1/K, K)
  }
  else{
    if(!sum(weight)==1){
      print("error! Input for 'weight' should be either a vector of standardized weights or use default option")
    }
    rank_fed = dat %*% weight
  }
  tbl$rank_fed = rank_fed
  return(tbl)
}




