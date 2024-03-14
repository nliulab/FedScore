########################################
########################################
########################################
# Internal functions for federated AutoScore
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# Functions for writing and reading json file

##### write cutoff.json
##### input: 
##### sitename
##### train data
write_cutoff <- function(df, quantiles = c(0, 0.05, 0.2, 0.8, 0.95, 1), max_cluster = 5, 
                         categorize = "quantile", sitename, path = "cutoff.json") {
  cutoff <- get_cut_vec(df, quantiles, max_cluster, categorize)
  cutoff_new <- list("site" = sitename, "cutoff" = cutoff)
  if (file.exists(path)) {
    all_cutoff <- rjson::fromJSON(file = path)
  } else {
    all_cutoff <- list()
  }
  K <- length(all_cutoff)
  found <- FALSE
  if (K > 0) {
    for (i in 1:K) {
      if (all_cutoff[[i]]$site == sitename) {
        all_cutoff[[i]] <- cutoff_new
        found <- TRUE
      }
    }
  }
  if (!found) {
    all_cutoff[[K + 1]] <- cutoff_new
  }
  json_data <- rjson::toJSON(all_cutoff, 2)
  write(json_data, path)
}

##### get_cutoff
##### input: 
##### sitename
##### train data
read_all_cutoff <- function(site_names, path = "cutoff.json") {
  tmp <- c()
  for (i in site_names) {
    tmp <- c(tmp,read_cutoff(i, path))
  }
  return(tmp)
}


##### read cutoff.json
##### input:
##### sitename
read_cutoff <- function(sitename, path = "cutoff.json") {
  all_cutoff <- rjson::fromJSON(file = path)
  K <- length(all_cutoff)
  for (i in 1:K) {
    if (all_cutoff[[i]]$site == sitename) {
      print(all_cutoff[[i]]$cutoff)
      return(all_cutoff[[i]]$cutoff)
    }
  }
  return()
}



##### write site_beta.json
##### input: 
##### sitename
##### modelIndex, range from 1 to 20
##### \hat{\beta_j}
##### \hat{Vj}
##### a vector for unique time points at this site
##### variableList: predictors included in the model
write_site_localest <- function(sitename, modelIndex, Lbeta, LSigma, t_unique, variableList, i_fold) {
  mylist <- list(site = sitename, model_index = modelIndex, beta_local = Lbeta, Sigma_local = LSigma, 
                 t_unique = t_unique, variable = variableList, i_fold = i_fold)
  if (file.exists("../../output/survival/site_beta.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/survival/site_beta.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  found <- FALSE
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$site == sitename & all_site_beta[[i]]$model_index == modelIndex & all_site_beta[[i]]$i_fold == i_fold) {
        all_site_beta[[i]] <- mylist
        found <- TRUE
      }
    }
  }
  if (!found) {
    all_site_beta[[K + 1]] <- mylist
  }
  json_data <- rjson::toJSON(all_site_beta, 2)
  write(json_data, "../../output/survival/site_beta.json")
}

##### function to input local results calculated by each site
##### file: site_beta.json
##### input:
##### sitename
##### data_matrix, col must be variable ranked by importance.
write_Lbeta_i <- function(sitename, data_train, var, cutoff, i_fold) {
  N <- length(var)
  data_transformed <- transform_df_fixed(data_train, cutoff)
  for (i in 0:(N - 1)) {
    k <- i + 1
    var_model <- var[1:k]
    mat <- data_transformed[, c("label_time", "label_status", var_model)]
    beta <- get_beta(mat)
    V <- get_V(mat)
    t_unique <- unique(mat$label_time)
    write_site_localest(sitename, modelIndex = k, beta, V, t_unique, var_model, i_fold)
    print(sprintf("Local estimations for model %d has been written into file", k))
  }
}

write_Lbeta <- function(sitename, data_train_all, var, cutoff, fold = "no_cv") {
  if (fold != "no_cv") {
    for (i in 1:fold) {
      data_train_temp <- get_fold(data_train_all, i)
      write_Lbeta_i(sitename, data_train_temp, var, cutoff, i)
    }
  } else {
    write_Lbeta_i(sitename, data_train_all, var, cutoff, fold)
  }
}

# calculate parameters for the final selected model
write_Lbeta_final <- function(sitename, data_train, var, cutoff, modelIndex) {
  if (modelIndex == "change_ref") {
    data_transformed <- transform_df_fixed(data_train, cutoff)
    mat <- data_transformed[, c("label_time", "label_status", var)]
  } else {
    mat <- data_train
  }
  N <- length(var)
  beta <- get_beta(mat)
  V <- get_V(mat)
  t_unique <- unique(mat$label_time)
  write_site_localest(sitename, modelIndex, beta, V, t_unique, var, i_fold = "NA")
  print(sprintf("Local beta for model %s has been written into file", modelIndex))
}


##### read site_beta.json (local beta)
##### input: sitename
##### return the optim beta by the site that is selected as local site
read_beta <- function(sitename = "SG", modelIndex, i_fold) {
  if (file.exists("../../output/survival/site_beta.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/survival/site_beta.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$site == sitename & all_site_beta[[i]]$model_index == modelIndex & all_site_beta[[i]]$i_fold == i_fold) {
        return(all_site_beta[[i]]$beta_local)
      }
    }
  }
  return()
}

##### read site_beta.json (local Sigma)
##### input: sitename
##### return the optim Sigma by the site that is selected as local site
read_Sigma <- function(sitename = "SG", modelIndex, i_fold) {
  if (file.exists("../../output/survival/site_beta.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/survival/site_beta.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$site == sitename & all_site_beta[[i]]$model_index == modelIndex & all_site_beta[[i]]$i_fold == i_fold) {
        return(all_site_beta[[i]]$Sigma_local)
      }
    }
  }
  return()
}

##### read site_beta.json (unique time points)
##### input: sitename
##### return the unique time points by the site that is selected as local site
read_tunique <- function(sitename = "SG", modelIndex, i_fold) {
  if (file.exists("../../output/survival/site_beta.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/survival/site_beta.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$site == sitename & all_site_beta[[i]]$model_index == modelIndex & all_site_beta[[i]]$i_fold == i_fold) {
        return(all_site_beta[[i]]$t_unique)
      }
    }
  }
  return()
}

##### write intermediate.json
##### input: 
##### sitenameL: local site name
##### sitenameR: remote site name
##### modelIndex, range from 1 to 20
##### D1, D2
##### N_R: the sample size for the inout remote site.
# Previous notes from FedScore:
# Note that D2 need to be store by rows in json file. Must be read by row accordingly
write_intermediate <- function(sitenameR, modelIndex, D1, D2, n, variableList, heterogeneity, i_fold) {
  mylist <- list(remote_site = sitenameR, model_index = modelIndex, D1=D1, D2=D2, N_R = n, 
                 variableList = variableList, heterogeneity = heterogeneity, i_fold = i_fold)
  if (file.exists("../../output/survival/intermediate.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/survival/intermediate.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  found <- FALSE
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$remote_site == sitenameR & all_site_beta[[i]]$model_index == modelIndex 
          & all_site_beta[[i]]$heterogeneity == heterogeneity & all_site_beta[[i]]$i_fold == i_fold) {
        all_site_beta[[i]] <- mylist
        found <- TRUE
      }
    }
  }
  if (!found) {
    all_site_beta[[K + 1]] <- mylist
  }
  json_data <- rjson::toJSON(all_site_beta, 2)
  write(json_data, "../../output/survival/intermediate.json")
}


##### read intermediate.json
##### input:
##### sitenameR: remote site name
##### modelIndex, range from 1 to 20
##### return D1 and D2
read_intermediate <- function(sitename, modelIndex, heterogeneity, i_fold) {
  if (file.exists("../../output/survival/intermediate.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/survival/intermediate.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  for (i in 1:K) {
    if (all_site_beta[[i]]$remote_site == sitename & all_site_beta[[i]]$model_index == modelIndex 
        & all_site_beta[[i]]$heterogeneity == heterogeneity & all_site_beta[[i]]$i_fold == i_fold) {
      return(all_site_beta[[i]])
    }
  }
  return()
}


################
################
##### function which lets remote site calculate intermediate values and write them into intermediate.json
##### file: site_beta.json, intermediate.json
##### input:
##### sitename
##### data_matrix, remote site data, col must be variable ranked by importance.
##### when data_matrix is available, there's no need to use i_fold?
write_Ds_i <- function(sitenameR, site_list, data_matrix, var, cutoff, heterogeneity = T, i_fold) {
  N <- length(var)
  choice <- ifelse(heterogeneity, "hete", "homo")
  data_transformed <- transform_df_fixed(data_matrix, cutoff)
  for (i in 0:(N - 1)) {
    k <- i + 1
    var_model <- var[1:k]
    mat <- data_transformed[, c('label_status','label_time', var_model)]
    Lbeta <- read_beta(sitename = sitenameR, modelIndex = k, i_fold = i_fold) %>% unlist # get info from site_beta.json
    # Calculate D
    if (heterogeneity) {
      Dev <- get_dev_ODACH(mat, k, sitenameR, site_list, i_fold)
    } else {
      Dev <- get_dev_ODAC(mat, k, sitenameR, site_list, i_fold)
    }
    D1 <- as.vector(Dev$logL_D1)
    D2 <- as.matrix(Dev$logL_D2)
    dimnames(D1) <- NULL
    dimnames(D2) <- NULL
    # write in intermediate.json
    write_intermediate(sitenameR = sitenameR, modelIndex = k,  variableList = var_model, D1 = D1, D2 = D2, 
                       n = nrow(mat), heterogeneity = choice, i_fold = i_fold)
    print(sprintf("D1 and D2 for model %d has been written into file", k))
  }
}

write_Ds <- function(sitenameR, site_list, data_train_all, var, cutoff, heterogeneity = T, fold = "no_cv") {
  if (fold != "no_cv") {
    for (i in 1:fold) {
      data_train_temp <- get_fold(data_train_all, i)
      write_Ds_i(sitenameR, site_list, data_train_temp, var, cutoff, heterogeneity, i)
    }
  } else {
    write_Ds_i(sitenameR, site_list, data_train_all, var, cutoff, heterogeneity, fold)
  }
}

# prepare parameters for final model
write_Ds_final <- function(sitenameL, site_list, data_train, var, cutoff, modelIndex, heterogeneity = T) {
  if (modelIndex == "change_ref") {
    data_transformed <- transform_df_fixed(data_train, cutoff)
    mat <- data_transformed[, c("label_time", "label_status", var)]
  } else {
    mat <- data_train
  }
  N <- length(site_list) - 1 
  choice <- ifelse(heterogeneity, "hete", "homo")
  site_other <- site_list[site_list != sitenameL]
  for (i in 1:N) {
    sitenameR <- site_other[i]
    Lbeta <- read_beta(sitename = sitenameL, modelIndex = modelIndex, i_fold = "NA") %>% unlist # get info from site_beta.json
    # Calculate D
    if (heterogeneity) {
      Dev <- get_dev_ODACH(mat, modelIndex, sitenameR, site_list, i_fold = "NA")
    } else {
      Dev <- get_dev_ODAC(mat, modelIndex, sitenameR, site_list, i_fold = "NA")
    }
    D1 <- as.vector(Dev$logL_D1)
    D2 <- as.matrix(Dev$logL_D2)
    dimnames(D1) <- NULL
    dimnames(D2) <- NULL
    # write in intermediate.json
    write_intermediate(sitenameR = sitenameR, modelIndex = modelIndex, variableList = var, D1 = D1, D2 = D2, 
                       n = nrow(mat), heterogeneity = choice, i_fold = "NA")
  }
}

# obtain and write beta_global and Hessian using data from one site.
write_surr_i <- function(data_matrix, sitename, site_list, var_list, heterogeneity, i_fold) {
  # loop through all m variables
  for (i in 1:length(var_list)) {
    print("starting model index:")
    print(i)
    var <- var_list[1:i]
    dat <- data_matrix[, c("label_time", "label_status", var)]
    res <- get_btilde(dat, sitename, site_list, i, heterogeneity, i_fold)
    write_tbeta(sitename, i, heterogeneity, res$converge, res$beta_global, res$H, i_fold)
  }
}

write_surr <- function(data_train_all, sitename, site_list, var_list, heterogeneity, uni_cut, fold = "no_cv") {
  if (fold != "no_cv") {
    for (i in 1:fold) {
      print("starting fold:")
      print(i)
      data_train_temp <- get_fold(data_train_all, i)
      data_train_temp <- transform_df_fixed(data_train_temp, uni_cut)
      write_surr_i(data_train_temp, sitename, site_list, var_list, heterogeneity, i)
    }
  } else {
    write_surr_i(transform_df_fixed(data_train_all, uni_cut), sitename, site_list, var_list, heterogeneity, fold)
  }
}


write_surr_final <- function(data_matrix, sitename, site_list, modelIndex = "final", heterogeneity, i_fold = "NA") {
  res <- get_btilde(data_matrix, sitename, site_list, modelIndex = "final", heterogeneity, i_fold)
  write_tbeta(sitename, "final", heterogeneity, res$converge, res$beta_global, res$H, i_fold)
}

# Write local estimated ~beta into intermediate.json for input of get_federatedCox()
write_tbeta <- function(sitenameL, modelIndex, heterogeneity, converge, beta_global, hessian, i_fold) {
  mylist <- list(site = sitenameL, model_index = modelIndex, heterogeneity = heterogeneity, 
                converge = converge, beta_global = beta_global, hessian = hessian, i_fold = i_fold)
  if (file.exists("../../output/survival/surr.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/survival/surr.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  found <- FALSE
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$site == sitenameL & all_site_beta[[i]]$model_index == modelIndex &
          all_site_beta[[i]]$heterogeneity == heterogeneity & all_site_beta[[i]]$i_fold == i_fold) {
        all_site_beta[[i]] <- mylist
        found <- TRUE
      }
    }
  }
  if (!found) {
    all_site_beta[[K + 1]] <- mylist
  }
  json_data <- rjson::toJSON(all_site_beta, 2)
  write(json_data, "../../output/survival/surr.json")
}

# read local estimated ~beta into intermediate.json for input of get_federatedCox()
read_tbeta <- function(sitename, modelIndex, heterogeneity, i_fold) {
  if (file.exists("../../output/survival/surr.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/survival/surr.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$site == sitename & all_site_beta[[i]]$model_index == modelIndex &
          all_site_beta[[i]]$heterogeneity == heterogeneity & all_site_beta[[i]]$i_fold == i_fold) {
        return(all_site_beta[[i]])
      }
    }
  }
  return()
}


######################################################
# reformat and write output of AutoScore_rank() to csv file
# rank: vector, output by AutoScore_rank()
# name: character, example: "imp1", "imp2",... assigned for each site
# rank: character, example: "rank1", "rank2".... assigned for each site accordingly
write_rank <- function(data, csv_name) {
  var <- names(data) %>% as.vector()
  tbl <- as.data.frame(cbind(var, data))
  tbl$rank <- c(1:length(data))
  colnames(tbl) <- c("var", "imp", "rank")
  write.csv(tbl, file = csv_name, row.names = F)
  print(sprintf("reformatted rank has been written into %s", csv_name))
}

######################################################
# get scores consist of integers from a coefficient 
transform_score <- function(vec, myrange = 100) {
  c1 <- range(vec)
  l <- c1[2] - c1[1]
  c2 <- vec * myrange / l
  return(round(c2))
}

######################################################
# transform score tbl to all positive, by adding a positive integer
transform_scoretbl <- function(vec) {
  return(vec - min(vec))
}

######################################################
# get unified cut off for all sites
# read from cutoff.json 
# return unified cutoff for variables
get_uni_cut <- function(site_list, path = "cutoff.json", weight= NULL, K) {
  cut_all <- vector(mode = "list", length = length(site_list))
  for (i in 1:length(site_list)) {
    cutoff <- read_cutoff(site_list[i], path)
    cut_all[[i]] <- cutoff
    var_name <- names(cutoff)
  }
  if (is.null(weight)) {
    weight <- rep(1/K, K)
  }
  cut_uni <- vector(mode = "list", length = length(var_name))
  # use loop to check if the length of cutoff are the same; add later
  for (i in 1:length(var_name)) {
    N <- length(cut_all[[1]][[i]])
    myvec <- matrix(nrow = length(site_list), ncol = N, byrow = T)
    for (j in 1:length(cut_all)) {
      if (typeof(cut_all[[j]][[i]]) == "character") {
        cut_uni[[i]] <- "let_binary"
        break
      } else {
        if (length(cut_all[[j]][[i]]) != N) {
          cut_all[[j]][[i]] <- cut_all[[1]][[i]]
        }
        myvec[j, ] <- as.vector(cut_all[[j]][[i]])}
    }
    cut_uni[[i]] <- colSums(myvec * weight) %>% round
  }
  names(cut_uni) <- var_name
  return(cut_uni)
}

######################################################
# new function to 
# transform continuous variable in df to categorical
# input: list read from cutoff.json
#        training data at local site (full)
#        variableList
transform_df <- function(dat, cutoff, variableList) {
  var_cts <- names(cutoff)
  var_other <-  variableList[!variableList %in% var_cts]
  df_temp <- dat[, var_cts]
  df <- c()
  for (i in 1:length(var_cts)) {
    if (typeof(cutoff[[i]]) == "character") {
      df <- cbind(df, df_temp[, i])
    } else {
      vec <- cut(df_temp[, i], breaks = cutoff[[i]], include.lowest = T)
      df <- cbind(df, vec)
    }
  }
  colnames(df) <- var_cts
  df <- cbind(df, dat[, var_other])
  head(df)
}


# Internal function in AutoScore:
#' @title Internal function: Calculate cut_vec from the training set (AutoScore Module 2)
#' @param df training set to be used for calculate the cut vector
#' @param categorize  Methods for categorize continuous variables. Options include "quantile" or "kmeans" (Default: "quantile").
#' @param quantiles Predefined quantiles to convert continuous variables to categorical ones. (Default: c(0, 0.05, 0.2, 0.8, 0.95, 1)) Available if \code{categorize = "quantile"}.
#' @param max_cluster The max number of cluster (Default: 5). Available if \code{categorize = "kmeans"}.
#' @return cut_vec for \code{transform_df_fixed}
get_cut_vec <-
  function(df,
           quantiles = c(0, 0.05, 0.2, 0.8, 0.95, 1),
           #by default
           max_cluster = 5,
           categorize = "quantile") {
    # Generate cut_vec for downstream usage
    cut_vec <- list()
    
    for (i in setdiff(names(df), c("label", "label_time", "label_status"))) {
      # for factor variable
      if (is.factor(df[, i])) {
        if (length(levels(df[, i])) < 10)
          (next)()
        else
          warning("WARNING: The number of categories should be less than 10",
                  i)
      }
      
      if (categorize == "quantile") {
        cut_off_tmp <- quantile(df[, i], quantiles, na.rm=T)
        cut_off_tmp <- unique(cut_off_tmp)
        cut_off <- signif(cut_off_tmp, 3)
      } else if (categorize == "k_means") {
        clusters <- kmeans(na.omit(df[, i]), max_cluster)
        cut_off_tmp <- c()
        for (j in unique(clusters$cluster)) {
          cut_off_tmp <- append(cut_off_tmp, min(df[, i][clusters$cluster == j], na.rm=T))
        }
        cut_off_tmp <- append(cut_off_tmp, max(df[, i], na.rm=T))
        cut_off_tmp <- sort(cut_off_tmp)
        cut_off_tmp <- unique(cut_off_tmp)
        cut_off <- signif(cut_off_tmp, 3)
        cut_off <- unique(cut_off)
      } else {
        stop('ERROR: please specify correct method for categorizing:  "quantile" or "k_means".')
      }
      
      l <- list(cut_off)
      names(l)[1] <- i
      cut_vec <- append(cut_vec, l)
    }
    ## delete min and max for each cut-off (min and max will be captured in the new dataset)
    if (length(cut_vec) != 0) {
      for (i in 1:length(cut_vec)) {
        if (length(cut_vec[[i]]) <= 2) {
          cut_vec[[i]] <- c("let_binary")
        } else {
          cut_vec[[i]] <- cut_vec[[i]][2:(length(cut_vec[[i]]) - 1)]
        }
      }
    }
    return(cut_vec)
    
  }



# Internal function in AutoScore
#' @title Internal function: Categorizing continuous variables based on cut_vec (AutoScore Module 2)
#' @param df dataset(training, validation or testing) to be processed
#' @param cut_vec fixed cut vector
#' @return  Processed \code{data.frame} after categorizing based on fixed cut_vec
#' @export
transform_df_fixed <- function(df, cut_vec) {
  j <- 1
  for (i in setdiff(names(df), c("label", "label_time", "label_status"))) {
    if (is.factor(df[, i])) {
      df[, i] <- factor(car::recode(var = df[, i], recodes = "NA = 'Unknown'"))
      if (length(levels(df[, i])) < 10)
        (next)()
      else
        stop("ERROR: The number of categories should be less than 9")
    }
    vec <- df[, i]
    cut_vec_new <- cut_vec[[j]]
    if (cut_vec_new[1] == "let_binary") {
      vec[vec != getmode(vec)] <- paste0("not_", getmode(vec))
      vec <- as.factor(vec)
      df[, i] <- vec
    } else {
      if (min(vec, na.rm=T) < cut_vec[[j]][1]) {
        cut_vec_new <- c(floor(min(df[, i], na.rm=T)) - 100, cut_vec_new)
      }
      if (max(vec, na.rm=T) >= cut_vec[[j]][length(cut_vec[[j]])]) {
        cut_vec_new <- c(cut_vec_new, ceiling(max(df[, i], na.rm=T) + 100))
      }
      
      cut_vec_new_tmp <- signif(cut_vec_new, 3)
      cut_vec_new_tmp <- unique(cut_vec_new_tmp)
      df_i_tmp <-  cut(
        df[, i],
        breaks = cut_vec_new_tmp,
        right = F,
        include.lowest = F,
        dig.lab = 3
      )

      ## delete min and max for the Interval after discretion: validation_set
      if (min(vec, na.rm=T) < cut_vec[[j]][1])
        levels(df_i_tmp)[1] <- gsub(".*,", "(,", levels(df_i_tmp)[1])
      if (max(vec, na.rm=T) >= cut_vec[[j]][length(cut_vec[[j]])])
        levels(df_i_tmp)[length(levels(df_i_tmp))] <-
        gsub(",.*", ",)", levels(df_i_tmp)[length(levels(df_i_tmp))])
      
      df[, i] <- as.factor(ifelse(is.na(df[, i]), "*Unknown", as.character(df_i_tmp)))
    }
    
    j <- j + 1
  }
  return(df)
}



# Internal function of AutoScore
#' @title Internal Function: Add baselines after second-step logistic regression (part of AutoScore Module 3)
#' @param df A \code{data.frame} used for logistic regression
#' @param coef_vec Generated from logistic regression
#' @return Processed \code{vector} for generating the scoring table
add_baseline <- function(df, coef_vec) {
  names(coef_vec) <- gsub("[`]", "", names(coef_vec)) # remove the possible "`" in the names
  df <- subset(df, select = names(df)[!names(df) %in% c("label", "label_time", "label_status")])
  coef_names_all <- unlist(lapply(names(df), function(var_name) {
    paste0(var_name, levels(df[, var_name]))
  }))
  coef_vec_all <- numeric(length(coef_names_all))
  names(coef_vec_all) <- coef_names_all
  # Remove items in coef_vec that are not meant to be in coef_vec_all
  # (i.e., the intercept)
  coef_vec_core <-
    coef_vec[which(names(coef_vec) %in% names(coef_vec_all))]
  i_coef <-
    match(x = names(coef_vec_core),
          table = names(coef_vec_all))
  coef_vec_all[i_coef] <- coef_vec_core
  coef_vec_all
}




# from pda package:
# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_coxph_logL <- function(beta, time, event, z) {
  .Call('_pda_rcpp_coxph_logL', PACKAGE = 'pda', beta, time, event, z)
}

rcpp_coxph_logL_gradient <- function(beta, time, event, z) {
  .Call('_pda_rcpp_coxph_logL_gradient', PACKAGE = 'pda', beta, time, event, z)
}

rcpp_coxph_logL_hessian <- function(beta, time, event, z) {
  .Call('_pda_rcpp_coxph_logL_hessian', PACKAGE = 'pda', beta, time, event, z)
}

rcpp_coxph_logL_efron <- function(beta, time, event, z) {
  .Call('_pda_rcpp_coxph_logL_efron', PACKAGE = 'pda', beta, time, event, z)
}

rcpp_coxph_logL_gradient_efron <- function(beta, time, event, z) {
  .Call('_pda_rcpp_coxph_logL_gradient_efron', PACKAGE = 'pda', beta, time, event, z)
}

rcpp_coxph_logL_gradient_efron_dist <- function(beta, ind_machine_, useLocal, dj_cutoff, time, event, z) {
  .Call('_pda_rcpp_coxph_logL_gradient_efron_dist', PACKAGE = 'pda', beta, ind_machine_, useLocal, dj_cutoff, time, event, z)
}

rcpp_aggregate <- function(x, indices, simplify = TRUE, cumulative = FALSE, reversely = FALSE) {
  .Call('_pda_rcpp_aggregate', PACKAGE = 'pda', x, indices, simplify, cumulative, reversely)
}



split_data_by_spec <- function(dat, spec = c(train = 0.6, validate = 0, test = 0.4), cv = T, fold = 5) {
  if (!cv) {
    g <- sample(cut(seq(nrow(dat)), nrow(dat) * cumsum(c(0, spec)), labels = names(spec)))
    dat_split <- split(dat, g)
  } else {
    total_samples <- nrow(dat)
    test_prop <- spec[["test"]]
    train_prop <- 1 - test_prop
    test_size <- round(total_samples * test_prop)
    train_size <- total_samples - test_size
    dat_split <- list()
    
    # first get testing set:
    test_indices <- sample(1:total_samples, test_size)
    dat_split$test <- dat[test_indices, ]
    
    # Training set
    train_indices <- setdiff(1:total_samples, test_indices)
    train_indices <- sample(train_indices) # random shuffle
    p <- sample(1:fold, size = length(train_indices), replace = T, prob = rep(1 / fold, fold))
    train_folds <- split(train_indices, p)
    dat_split$train <- lapply(1:fold, function(f) dat[train_folds[[f]], ])
    # all_train set
    dat_split$all_train <- dat[train_indices, ]
  }
  return(dat_split)
}

# get design matrix for a given data_matrix
get_designX <- function(data_matrix) {
  model <- coxph(Surv(label_time, label_status) ~ ., data_matrix) 
  return(model.matrix(model))
}

# get cut list that match with the order of variable ranking:
get_ordered_cut <- function(input_vector, uni_cut) {
  matching_names <- intersect(input_vector, names(uni_cut))
  return(uni_cut[match(matching_names, names(uni_cut))])
}

get_fold <- function(mylist, ifold) {
  dat <- mylist[-ifold]
  dat <- do.call(rbind, dat)
  return(dat)
}

# keep records of the cutoff frequency table
get_cutoff_freq <- function(data_train, file_path, fold, uni_cut) {
  if (file.exists(file_path)) {
    file.remove(file_path)
  }
  file.create(file_path)
  for (i in 1:fold) {
    data_train_temp <- get_fold(data_train, i)
    data_train_temp <- transform_df_fixed(data_train_temp, uni_cut)
    capture.output(paste0("starting fold:", as.character(i)), file=file_path, append = TRUE)
    capture.output(summary(data_train_temp[, names(order_uni_cut)]), file = file_path, append = TRUE)
  }
  print(paste0("freq of cutoff is kept in ", file_path))
}

# For parsimony plot (all outcomes) -----

#' Internal function: Make parsimony plot
#' @inheritParams AutoScore_parsimony
#' @param AUC A vector of AUC values (or mAUC for ordinal outcomes).
#' @param variables A vector of variable names
#' @param num A vector of indices for AUC values to plot. Default is to plot all.
#' @param ylab Title of y-axis
#' @param title Plot title
#' @import ggplot2
plot_auc <- function(AUC, variables, num = seq_along(variables), auc_lim_min, auc_lim_max,
                     ylab = "Mean Area Under the Curve", title = "Parsimony plot on the validation set") {
  if (auc_lim_max == "adaptive") {
    auc_lim_max <- max(AUC)
  }
  dt <- data.frame(AUC = AUC, variables = factor(variables, levels = variables), num = num)
  p <- ggplot(data = dt, mapping = aes_string(x = "variables", y = "AUC")) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_cartesian(ylim = c(auc_lim_min, auc_lim_max)) +
    theme_bw() +
    labs(x = "", y = ylab, title = title) +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # Add number of variables to bar (resolves previous warning regarding dt$num):
  if (nrow(dt) >= 100) {
    p + geom_text(data = dt, aes_string(label = "num"), vjust = 1.5, colour = "white", angle = 90)
  } else {
    p + geom_text(data = dt, aes_string(label = "num"), vjust = 1.5, colour = "white")
  }
}

## get fed auc_file
extract <- function(input_file) {
  content <- readLines(input_file)
  
  lines1 <- grep("Integrated AUC ", content, value = TRUE)
  lines2 <- grep("generating auc", content, value = TRUE)
  
  res <- "AUC,CI,Site,Model,fold"
  
  for (i in seq_along(lines1)) {
    line1 <- lines1[i]
    splitted1 <- strsplit(line1, " ")[[1]]
    auc <- splitted1[length(splitted1) - 1]
    ci <- substr(splitted1[length(splitted1)], 2, nchar(splitted1[length(splitted1)]) - 1)
    
    line2 <- lines2[i]
    splitted2 <- strsplit(line2, " ")[[1]]
    site <- splitted2[6]
    model <- substr(splitted2[10], 1, nchar(splitted2[10]) - 1)
    fold <- substr(splitted2[12], 1, nchar(splitted2[12]) - 1)
    
    res <- paste(res, paste(auc, ci, site, model, fold, sep = ","), sep = "\n")
  }
  
  return(res)
}

extract_multi <- function(input_file) {
  content <- readLines(input_file)
  lines1 <- grep("local model ", content, value = TRUE)
  lines2 <- grep("Integrated AUC ", content, value = TRUE)
  lines3 <- grep("C_index:", content, value = TRUE)
  res <- data.frame(matrix(ncol = 7, nrow = 0))
  lines2 <- lines2[3:length(lines2)]
  lines3 <- lines3[3:length(lines3)]
  for (i in seq_along(lines1)) {
    line1 <- lines1[i]
    splitted1 <- strsplit(line1, " ")[[1]]
    model_site <- splitted1[5]
    test_site <- splitted1[10]
    var <- splitted1[14]
    
    line2 <- lines2[i]
    splitted2 <- strsplit(line2, " ")[[1]]
    auc <- splitted2[length(splitted2) - 1]
    auc_ci <- substr(splitted2[length(splitted2)], 2, nchar(splitted2[length(splitted2)]) - 1)
    
    line3 <- lines3[i]
    splitted3 <- strsplit(line3, " ")[[1]]
    c_index <- splitted3[length(splitted3) - 1]
    c_index_ci <- substr(splitted3[length(splitted3)], 2, nchar(splitted3[length(splitted3)]) - 1)
    print(c(model_site, test_site, var, auc, auc_ci, c_index, c_index_ci))
    res <- rbind(res, c(model_site, test_site, var, auc, auc_ci, c_index, c_index_ci))
  }
  colnames(res) <-  c("model_site", "test_site", "var", "auc", "auc_ci", "c_index", "c_index_ci")
  return(res)
}

extract_final <- function(input_file, output_file) {
  content <- readLines(input_file)
  
  lines1 <- grep("Integrated AUC ", content, value = TRUE)
  lines2 <- grep("generating auc", content, value = TRUE)
  lines3 <- grep("C_index:", content, value = TRUE)
  
  res <- data.frame(matrix(ncol = 3, nrow = 0))
  for (i in seq_along(lines1)) {
    line1 <- lines1[i]
    splitted1 <- strsplit(line1, " ")[[1]]
    auc <- splitted1[length(splitted1) - 1]
    ci <- substr(splitted1[length(splitted1)], 2, nchar(splitted1[length(splitted1)]) - 1)
    
    line2 <- lines2[i]
    splitted2 <- strsplit(line2, " ")[[1]]
    site <- splitted2[6]
    
    line3 <- lines3[i]
    splitted3 <- strsplit(line3, " ")[[1]]
    c_index <- splitted3[length(splitted3) - 1]
    c_ci <- substr(splitted3[length(splitted3)], 2, nchar(splitted3[length(splitted3)]) - 1)
    
    res <-rbind(res, c(site, auc, ci, c_index, c_ci))
  }
  colnames(res) <- c("Site","AUC","CI","Cindex","Cindex_CI")
  write.csv(res, output_file, row.names = FALSE)
  return(res)
}

## parsimonious plot
parsimony_plot_by_site <- function(auc_df, site_name, nfold = "no_cv", output_path) {
  pdf(file = output_path)
  tmp <- auc_df %>% filter(., Site == site_name)
  if (nfold != "no_cv") {
    opar <- par(no.readonly = TRUE)
    par(mfrow = c(3, 2))
    for(j in 1:nfold) {
      AUC <- tmp %>% filter(., fold==j) %>% select(AUC)
      AUC <- AUC$AUC
      names(AUC) <- 1:length(AUC)
      plot(
        AUC,
        main = paste("Parsimony plot (cross validation) for fold", j),
        xlab = "Number of Variables",
        ylab = "Area Under the Curve",
        col = "#2b8cbe",
        lwd = 2,
        type = "o",
        ylim = c(0.5, 1)
      )
    }
    par(opar[c("mfrow", "mar", "oma", "mgp", "cex", "col", "font")])
    par()
  }
  tmp_df_avg <- auc_df %>% group_by(Site, Model) %>% summarize(mean_auc = mean(AUC), .groups = "drop") %>% arrange(Site, Model)
  tmp <- tmp_df_avg %>% filter(., Site == site_name)
  AUC <- tmp$mean_auc
  plot(plot_auc(AUC = AUC, variables = var.rank,
                num = 1:length(AUC),
                auc_lim_min = 0,
                auc_lim_max = "adaptive",
                ylab = "Area Under the Curve",
                title = "Parsimony plot on the validation set"))
  dev.off()
}


eva_performance_iauc <- function(score, validation_set, print = TRUE) {
  Surv.rsp.new <- Surv(validation_set$label_time, validation_set$label_status)

  km_fit_test <- survfit(Surv.rsp.new ~ 1, data = validation_set)
  km_time <- summary(km_fit_test)$time
  km_survival <- summary(km_fit_test)$surv
  km_time <- km_time[-length(km_time)]
  km_survival <- km_survival[-length(km_survival)]
  
  AUC_uno <- AUC.uno(Surv.rsp.new, Surv.rsp.new, lpnew = score, times = km_time)
  km_survival_w <- c(1, km_survival[-length(km_survival)])
  km_survival_sum <- sum(km_survival_w - km_survival)
  weight <- (km_survival_w - km_survival) / km_survival_sum
  iAUC <- sum(weight * AUC_uno$auc)
  if (print) {
    cat(iAUC)
    cat("\n")
  }
  return(iAUC)
}





#' @title Internal Function: Change Reference category after first-step logistic regression (part of AutoScore Module 3)
#' @param df A \code{data.frame} used for logistic regression
#' @param coef_vec Generated from logistic regression
#' @return Processed \code{data.frame} after changing reference category
change_reference <- function(df, coef_vec) {
  # delete label first
  df_tmp <- subset(df, select = names(df)[!names(df) %in% c("label", "label_time", "label_status")])
  names(coef_vec) <- gsub("[`]", "", names(coef_vec)) # remove the possible "`" in the names
  
  # one loops to go through all variable
  for (i in (1:length(df_tmp))) {
    var_name <- names(df_tmp)[i]
    var_levels <- levels(df_tmp[, i])
    var_coef_names <- paste0(var_name, var_levels)
    coef_i <- coef_vec[which(names(coef_vec) %in% var_coef_names)]
    if (min(coef_i) < 0) {
      ref <- var_levels[which(var_coef_names == names(coef_i)[which.min(coef_i)])]
      print(ref)
      df_tmp[, i] <- relevel(df_tmp[, i], ref = ref)
    }
  }
  
  # add label again
  if(!is.null(df$label))  df_tmp$label <- df$label
  else{
    df_tmp$label_time <- df$label_time
    df_tmp$label_status <- df$label_status}
  return(df_tmp)
}

find_reference <- function(sitenameL, site_list, data_train, var, cutoff, modelIndex, heterogeneity = T) {
  N <- length(site_list) - 1 
  site_other <- site_list[site_list != sitenameL]
  data_transformed <- transform_df_fixed(data_train, cutoff)
  mat <- data_transformed[, c('label_status','label_time', var)]
  for (i in 1:N) {
    sitenameR <- site_other[i]
    Lbeta <- read_beta(sitename = sitenameL, modelIndex = modelIndex, i_fold = "NA") %>% unlist # get info from site_beta.json
    # Calculate D
    if (heterogeneity) {
      Dev <- get_dev_ODACH(mat, modelIndex, sitenameR, site_list, i_fold = "NA")
    } else {
      Dev <- get_dev_ODAC(mat, modelIndex, sitenameR, site_list, i_fold = "NA")
    }
    D1 <- as.vector(Dev$logL_D1)
    D2 <- as.matrix(Dev$logL_D2)
    dimnames(D1) <- NULL
    dimnames(D2) <- NULL
    return(Dev$bbar)
  }
}

change_to_newref <- function(train_set, variable_list, cut_vec, ref) {
  train_set_1 <- transform_df_fixed(train_set, cut_vec)
  train_set_1 <- train_set_1[, c("label_time", "label_status", variable_list)]
  train_set_2 <- change_reference(train_set_1, ref)
  return(train_set_2)
}

create_score_tbl <- function(train_set, max_score, variable_list, coef_vec) {
  # rounding for final scoring table "score_table"
  coef_vec_tmp <- round(coef_vec / min(coef_vec))
  score_table <- add_baseline(train_set, coef_vec_tmp)
  
  # normalization according to "max_score" and regenerate score_table
  total_max <- max_score
  total <- 0
  for (i in 1:length(variable_list))
    total <- total + max(score_table[grepl(variable_list[i], names(score_table))])
  score_table <- round(score_table / (total / total_max))
  return(score_table)
}

#' @title Internal Function: Automatically assign scores to each subjects given new data set and scoring table (Used for intermediate and final evaluation)
#' @param df A \code{data.frame} used for testing, where variables keep before categorization
#' @param score_table A \code{vector} containing the scoring table
#' @return Processed \code{data.frame} with assigned scores for each variables
assign_score <- function(df, score_table) {
  for (i in setdiff(names(df), c("label", "label_time", "label_status"))) {
    score_table_tmp <- score_table[grepl(i, names(score_table))]
    df[, i] <- as.character(df[, i])
    for (j in 1:length(names(score_table_tmp))) {
      df[, i][df[, i] %in% gsub(i, "", names(score_table_tmp)[j])] <- score_table_tmp[j]
    }
    df[, i] <- as.numeric(df[, i])
  }
  return(df)
}

show_score_table <- function(test_set,
                           coef_vec,
                           final_variables,
                           max_score = 100) {
  score_tbl <- create_score_tbl(test_set, max_score, final_variables, coef_vec)
  print_scoring_table(scoring_table = score_tbl, final_variable = var_final)
  return(score_tbl)
}
FedScore_survival_testing <- function(test_set, score_tbl, final_variables, time_point = c(1, 3, 7, 14, 30)) {
  set.seed(42)
  # Using "assign_score" to generate score based on new dataset and Scoring table "score_table"
  print_scoring_table(scoring_table = score_tbl, final_variable = var_final)
  test_score <- assign_score(test_set, score_tbl)
  test_score$total_score <- rowSums(subset(test_score, select = names(test_score)[names(test_score) != c("label_status") 
                                                                                  & names(test_score) != c("label_time")]))
  print_performance_ci_survival(test_score$total_score, test_score, time_point)
}

Get_auc_bySite <- function(dat1_train, nsite, var_rank, uni_cut, nmodel, nfold = "no_cv", site_index = 0, heterogeneity = T) {
  if (nfold == "no_cv") {
    nfold <- 1
    dat1_train[[1]] <- dat1_train
  }
  tmp_matrix <- matrix(0, nrow = nmodel*nfold, ncol = 4)
  colnames(tmp_matrix) <- c("AUC","Site","Model","fold")
  t <- 0
  for (i in 1:nmodel) {
    var <- var_rank[1:i]
    # loop through all folds:
    for (j in 1:nfold) {
      t <- t + 1
      print(sprintf("generating auc for Site %d validation for model %d, fold %d", nsite, i, j))
      validate_temp <- transform_df_fixed(dat1_train[[j]], uni_cut)
      if (heterogeneity) {
        csv_name <- sprintf("../../output/survival/fedsurv_hete/coef_global_hete_trainCV5_model%d_fold%d.csv", i, j)
      } else {
        csv_name <- sprintf("../../output/survival/fedsurv/coef_global_homo_trainCV5_model%d_fold%d.csv", i, j)
      }
      coef <- read.csv(csv_name)
      design_valid_temp <- validate_temp[, c('label_time', 'label_status', var)]
      designX_valid <- get_designX(design_valid_temp)
      iAUC <- eva_performance_iauc(designX_valid %*% coef$V1, validate_temp, print = FALSE)
      tmp_matrix[t,] <- c(round(iAUC, 3), paste0("Site", nsite), i, j)
    }
  }
  tmp_df <- as.data.frame(tmp_matrix)
  tmp_df <- tmp_df %>% mutate(across(c(AUC, Model, fold), as.numeric))
  return(tmp_df)
}

AutoScore_testing_Survival_prc <- function(test_set, final_variables, cut_vec, scoring_table, threshold = "best", 
                                           with_label = TRUE, time_point = c(1, 3, 7, 14, 30, 60, 90)) {
    if (with_label) {
      # prepare test set: categorization and "assign_score"
      test_set_1 <- test_set[, c(final_variables, c("label_time","label_status"))]
      test_set_2 <- transform_df_fixed(test_set_1, cut_vec = cut_vec)
      test_set_3 <- assign_score(test_set_2, scoring_table)
      test_set_3$total_score <- rowSums(subset(test_set_3, 
                                               select = names(test_set_3)[names(test_set_3) != c("label_status") & 
                                                                            names(test_set_3) != c("label_time")]))
      test_set_3$total_score[which(is.na(test_set_3$total_score))] <- 0
      
      # Final evaluation based on testing set
      cat("***Performance using AutoScore (based on unseen test Set):\n")
      print_performance_ci_survival(test_set_3$total_score, test_set_3, time_point)
      Modelprc <- pr.curve(test_set_3$total_score[which(y_test == 1)], test_set_3$total_score[which(y_test == 0)], curve = TRUE)
      values <- coords(model_roc, "best", ret = c("specificity", "sensitivity", "accuracy", "npv", "ppv", "precision"), transpose = TRUE)
      pred_score <- data.frame(pred_score = test_set_3$total_score, label_time = test_set_3$label_time, label_status = test_set_3$label_status)
      return(pred_score)
    } else {
      test_set_1 <- test_set[, c(final_variables)]
      test_set_2 <- transform_df_fixed(test_set_1, cut_vec = cut_vec)
      test_set_3 <- assign_score(test_set_2, scoring_table)
      test_set_3$total_score <- rowSums(subset(test_set_3, 
                                               select = names(test_set_3)[names(test_set_3) != c("label_status") &
                                                                            names(test_set_3) != c("label_time")]))
      test_set_3$total_score[which(is.na(test_set_3$total_score))] <- 0
      pred_score <- data.frame(pred_score = test_set_3$total_score, label_time = NA, label_status = NA)
      return(pred_score)
    }
  }

rank_csv_to_Rvec <- function(tmp_rank) {
  tmp_rk <- tmp_rank$imp
  names(tmp_rk) <- tmp_rank$var
  return(tmp_rk)
}
