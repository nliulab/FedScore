# This is a demo using two subsets of `sample_data` from the AutoScore package.
# The outcome used in this demo is inpatient mortality (binary).

library(AutoScore)
library(tidyverse)
library(ggplot2)
library(mle.tools)
library(rjson)
library(rstudioapi)

setwd(dirname(getSourceEditorContext()$path))
source("../../code/R/fedscore_binary.R")
source("../../code/R/helpers_binary.R")

##### Import data
# Read data for two sites (split into train, validation and test sets)
dat1_split <- readRDS("../../data/binary/site1.rds")
dat2_split <- readRDS("../../data/binary/site2.rds")


##### Module 1: Obtain federated ranking
# each site first obtain and write ranks independently:
rank.s1 <- AutoScore::AutoScore_rank(dat1_split$train)
write_rank(rank.s1, "../../output/binary/ranks/rank1.csv")
rank.s2 <- AutoScore::AutoScore_rank(dat2_split$train)
write_rank(rank.s2, "../../output/binary/ranks/rank2.csv")
# then get global rank:
rank.fed <- AutoScore_Frank(K = 2, path = "../../output/binary/ranks/") # by using results generated by RF at each site

# take a look at the rank:
rank.fed[,c("var", "rank1", "rank2", "rank_fed")]

# get the variable list in the order of federated ranking:
var.rank <- rank.fed$var[1:9]


##### Module 2: Obtain unified cutoff vectors
# each site independently write cut off into cutoff.json
write_cutoff(df = dat1_split$train, sitename = "homo_S1", path = "../../output/binary/cutoff.json")
write_cutoff(df = dat2_split$train, sitename = "homo_S2", path = "../../output/binary/cutoff.json")
# decide uniform cutoff:
uni_cut <- get_uni_cut(c("homo_S1", "homo_S2"), path = "../../output/binary/cutoff.json")
uni_cut


##### Module 3: Conduct federated logistic regression
# Use site 1 as local site
names_all_site <- c("homo_S1", "homo_S2")
cutoff <- uni_cut
# 1. write site_beta.json independently at each site:
write_Lbeta("homo_S1", names_all_site, dat1_split$train, var.rank, cutoff) # site1
write_Lbeta("homo_S2", names_all_site, dat2_split$train, var.rank, cutoff) # site2
# 2. write site_intermediate.json:
write_Ds("homo_S1", "homo_S2", names_all_site, dat2_split$train, var.rank, cutoff)
write_Ds("homo_S2", "homo_S1", names_all_site, dat1_split$train, var.rank, cutoff)
# 3. do federated LR
dat1_train_transformed <- transform_df_fixed(dat1_split$train, uni_cut)
# loop though all variables:
for (i in 1:9) {
  print(i)
  var_m <- var.rank[1:i]
  model <- get_federatedLR(method = 1, data_matrix = dat1_train_transformed[,c("label", var_m)], 
                           sitenameL = "homo_S1", sitenameR_list = "homo_S2", modelIndex = i)
  coef <- model$beta_global %>% as.data.frame()
  csv_name <- sprintf("../../output/binary/coef_global_homoS1_train_model%d.csv", i)
  write.csv(coef, csv_name)
}


##### Module 4 & 5: Read the coefficients estimated by federated logistic regression to get score tables. Then test the model on test sets.

# Reorder the cut vector. The order of cut vector must be consistent with the order of variable ranking
uni_cut_fed <- list(uni_cut$Lab_H, uni_cut$Age, uni_cut$Lab_B, uni_cut$Lab_K, uni_cut$Vital_E, 
                    uni_cut$Lab_I, uni_cut$Vital_B, uni_cut$Lab_A, uni_cut$Lab_J)
names(uni_cut_fed) <- c("Lab_H", "Age", "Lab_B", "Lab_K", "Vital_E", "Lab_I", "Vital_B", "Lab_A", "Lab_J")

# get validation results
sink("../../output/binary/homoS1_auc_validation_fed.txt")
# loop through all 9 variables
for (i in 1:9) {
  var <- var.rank[1:i]
  csv_name <- sprintf("../../output/binary/coef_global_homoS1_train_model%d.csv", i)
  coef <- read.csv(csv_name)
  coef_vec <- transform_score(coef$.)
  names(coef_vec) <- coef$X
  df <- dat1_train_transformed[, var] %>% as.data.frame()
  colnames(df) <- var
  score_tbl <- add_baseline(df, transform_scoretbl(coef_vec))
  print(sprintf("generating auc for Site1 validation for model %d", i))
  pred <- AutoScore::AutoScore_testing(dat1_split$validate, var, uni_cut_fed, score_tbl, threshold = "best", with_label = TRUE)
  print(sprintf("generating auc for Site2 validation for model %d", i))
  pred <- AutoScore::AutoScore_testing(dat2_split$validate, var, uni_cut_fed, score_tbl, threshold = "best", with_label = TRUE)
}
sink()

# Use extract_fed.py to extract results from text to csv
# Sample command: python3 ../../code/Python/extract_fed.py ../../output/binary/homoS1_auc_validation_fed.txt ../../output/binary/fed.csv

### Plot parsimony plot from fed.csv
### site1
df <- read.csv("../../output/binary/fed.csv")
# We offer following two options for visualizations. 
# Users may set up other ways; see more discussion in the manuscript

# plot 1
df %>% group_by(Model) %>% 
  summarize(AUC_mean = mean(AUC)) %>% 
  ggplot(aes(x=Model, y = AUC_mean)) +  
  geom_line() + 
  geom_point()

# plot 2
df_new <- df %>% group_by(Model) %>% summarize(AUC_mean = mean(AUC))
df_new <- cbind(df_new, var.rank[1:9])
colnames(df_new) <- c("Model", "AUC_mean", "Var")

# use bar plots to stay consistent with AutoScore
var_names <- factor((var.rank)[1:9], levels = (var.rank)[1:9])
dt <- data.frame(AUC_mean = df_new$AUC_mean, variables = var_names, num = 1:9)
auc_lim_max <- ceiling(max(df$AUC)*10)/10
p <- ggplot(data = dt, mapping = aes_string(x = "variables", y = "AUC_mean")) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_cartesian(ylim = c(0, auc_lim_max)) +
  theme_bw() +
  labs(x = "", y = "Area Under the Curve", title = "Parsimony plot on the validation set") +
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(size=25, hjust = 0.5))
print(p + geom_text(aes(label = dt$num), vjust = 1.5, colour = "white"))

# Suppose we decide to select the first 6 variables for final model
# This is a subjective decision made by users, as discussed in the manuscript.

# Now fit new (final) logistic regression based on results of parsimony plot
sink("../../output/binary/fed_AUC_S1_with_scoretbl_var6.txt")

var_fed <- var.rank[1:6]
coef <- read.csv("../../output/binary/coef_global_homoS1_train_model6.csv")
coef_vec <- transform_score(coef$.)
names(coef_vec) <- coef$X
df <- dat1_train_transformed[, var_fed] %>% as.data.frame()
colnames(df) <- var_fed
score_tbl_fed <- add_baseline(df, transform_scoretbl(coef_vec))
print("printing score tbl")
score_tbl_fed %>% print

# test on testing sets independently:
print("generating auc, testing set = site1")
pred_site1 <- AutoScore::AutoScore_testing(dat1_split$test, var_fed, uni_cut_fed, score_tbl_fed, threshold = "best", with_label = TRUE)
print("generating auc, testing set = site2")
pred_site2 <- AutoScore::AutoScore_testing(dat2_split$test, var_fed, uni_cut_fed, score_tbl_fed, threshold = "best", with_label = TRUE)
sink()
