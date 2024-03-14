# Internal functions for federated AutoScore

# reformat and write output of AutoScore_rank() to csv file
# rank: vector, output by AutoScore_rank()
# name: character, example: "imp1", "imp2",... assigned for each site
# rank: character, example: "rank1", "rank2".... assigned for each site accordingly
write_rank <- function(data, csv_name) {
  var <- names(data) %>% as.vector()
  tbl <- as.data.frame(cbind(var,data))
  tbl$rank = c(1:length(data))
  colnames(tbl) <- c("var", "imp", "rank")
  write.csv(tbl, file = csv_name, row.names = F)
  print(sprintf("reformatted rank has been written into %s", csv_name))
}


# get scores consist of integers from a coefficient 
transform_score <- function(vec, myrange = 100) {
  c1 <- range(vec)
  l <- c1[2] - c1[1]
  c2 <- vec * myrange / l
  return(round(c2))
}


# transform score tbl to all positive, by adding a positive integer
transform_scoretbl <- function(vec) {
  return(vec - min(vec))
}


# get unified cut off for all sites
# read from cutoff.json 
# return unified cutoff for variables
get_uni_cut <- function(site_list, path = "cutoff.json") {
  cut_all <- vector(mode = "list", length = length(site_list))
  for (i in 1:length(site_list)) {
    cutoff <- read_cutoff(site_list[i], path)
    cut_all[[i]] <- cutoff
    list_name <- names(cutoff)
  }
  cut_uni <- vector(mode = "list", length = length(list_name))
  for (i in 1:length(list_name)) {
    N <- length(cut_all[[1]][[i]])
    myvec <- rep(0, N)
     for (j in 1:length(cut_all)) {
       if (typeof(cut_all[[j]][[i]]) == "character") {
         cut_uni[[i]] <- "let_binary"
         break
       } else {
         myvec <- myvec + as.vector(cut_all[[j]][[i]])
       }
       cut_uni[[i]] <- floor(myvec / length(cut_all)) 
     }
  }
  names(cut_uni) <- list_name
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
  var_other <- variableList[!variableList %in% var_cts]
  df_temp <- dat[, var_cts]
  df <- c()
  for (i in 1:length(var_cts)) {
    if (typeof(cutoff[[i]]) == "character") {
      df  <- cbind(df, df_temp[, i])
    } else {
      vec <- cut(df_temp[,i], breaks = cutoff[[i]], include.lowest = T)
      df <- cbind(df,vec)
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
           max_cluster = 5,
           categorize = "quantile") {
    # Generate cut_vec for downstream usage
    cut_vec <- list()
    for (i in setdiff(names(df), c("label", "label_time", "label_status"))) {
      # for factor variable
      if (is.factor(df[, i])) {
        if (length(levels(df[, i])) < 10) {
          (next)()
        } else {
          warning("WARNING: The number of categories should be less than 10", i)
        }
      }
      
      ## mode 1 - quantile
      if (categorize == "quantile") {
        cut_off_tmp <- quantile(df[, i], quantiles, na.rm = T)
        cut_off_tmp <- unique(cut_off_tmp)
        cut_off <- signif(cut_off_tmp, 3)
        
      ## mode 2 k-means clustering
      } else if (categorize == "k_means") {
        clusters <- kmeans(na.omit(df[, i]), max_cluster)
        cut_off_tmp <- c()
        for (j in unique(clusters$cluster)) {
          cut_off_tmp <- append(cut_off_tmp, min(df[, i][clusters$cluster == j], na.rm = T))
        }
        cut_off_tmp <- append(cut_off_tmp, max(df[, i], na.rm = T))
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
    if (length(cut_vec) != 0) { ## in case all the variables are categorical
      for (i in 1:length(cut_vec)) {
        if (length(cut_vec[[i]]) <= 2)
          cut_vec[[i]] <- c("let_binary")
        else
          cut_vec[[i]] <- cut_vec[[i]][2:(length(cut_vec[[i]]) - 1)]
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
  # for loop going through all variables
  for (i in setdiff(names(df), c("label", "label_time", "label_status"))) {
    if (is.factor(df[, i])) {
      df[, i] <- factor(car::recode(var = df[, i], recodes = "NA = 'Unknown'"))
      if (length(levels(df[, i])) < 10)
        (next)()
      else
        stop("ERROR: The number of categories should be less than 9")
    }
    
    ## make conresponding cutvec for validation_set: cut_vec_new
    vec <- df[, i]
    cut_vec_new <- cut_vec[[j]]
    
    if (cut_vec_new[1] == "let_binary") {
      vec[vec != getmode(vec)] <- paste0("not_", getmode(vec))
      vec <- as.factor(vec)
      df[, i] <- vec
    } else{
      if (min(vec, na.rm=T) < cut_vec[[j]][1])
        cut_vec_new <- c(floor(min(df[, i], na.rm = T)) - 100, cut_vec_new)
      if (max(vec, na.rm=T) >= cut_vec[[j]][length(cut_vec[[j]])])
        cut_vec_new <- c(cut_vec_new, ceiling(max(df[, i], na.rm = T) + 100))
      
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
  # Remove items in coef_vec that are not meant to be in coef_vec_all (i.e., the intercept)
  coef_vec_core <-
    coef_vec[which(names(coef_vec) %in% names(coef_vec_all))]
  i_coef <-
    match(x = names(coef_vec_core),
          table = names(coef_vec_all))
  coef_vec_all[i_coef] <- coef_vec_core
  coef_vec_all
}


##### read site_beta.json (local beta)
##### input: sitename
##### return the optim beta by the site that is selected as local site
read_beta <- function(sitename = "SG", modelIndex) {
  if (file.exists("../../output/binary/site_beta.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/binary/site_beta.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$site == sitename & all_site_beta[[i]]$model_index == modelIndex) {
        return(all_site_beta[[i]]$beta_local)
      }
    }
  }
  return()
}

##### read site_beta.json (local Sigma)
##### input: sitename
##### return the optim Sigma by the site that is selected as local site
read_Sigma <- function(sitename = "SG", modelIndex) {
  if (file.exists("../../output/binary/site_beta.json")) {
    all_site_beta <- rjson::fromJSON(file = "../../output/binary/site_beta.json")
  } else {
    all_site_beta <- list()
  }
  K <- length(all_site_beta)
  if (K > 0) {
    for (i in 1:K) {
      if (all_site_beta[[i]]$site == sitename & all_site_beta[[i]]$model_index == modelIndex) {
        return(all_site_beta[[i]]$Sigma_local)
      }
    }
  }
  return()
}

