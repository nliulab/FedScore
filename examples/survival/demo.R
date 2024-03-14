set.seed(42)

library(tidyverse)
library(doParallel)
library(foreach)
library(dplyr)
library(survival)
library(data.table)
library(pda)
library(AutoScore)
library(survAUC)
library(ggplot2)
library(rstudioapi)

setwd(dirname(getSourceEditorContext()$path))
source("../../code/R/helpers_survival.R")
source("../../code/R/fedscore_survival.R")

split_all_site <- readRDS("../../data/survival/demo.rds")
n_site <- length(split_all_site)
names_all_site <- paste0("homo_S", 1:n_site)

for (i in 1:n_site) {
 tmp_dat <- split_all_site[[i]]
 assign(paste0("dat", i, "_split"), tmp_dat)
}

colnames(split_all_site[[1]]$test)

subfolder_name <- "../../output/survival/rank"
if (!file.exists(subfolder_name)) {
 dir.create(subfolder_name)
 cat(paste("Subfolder '", subfolder_name, "' created successfully.\n", sep = ""))
} else {
 cat(paste("Subfolder '", subfolder_name, "' already exists.\n", sep = ""))
}

rank_all_site <- paste0("rank.s", 1:n_site)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
foreach (i = 1:n_site) %dopar% {
 library(tidyverse)
 rank_all_site <- paste0("rank.s", 1:n_site)
 tmp_dat_all <- readRDS("../../data/survival/demo.rds")
 tmp_dat <- tmp_dat_all[[i]]
 tmp_rank <- AutoScore::AutoScore_rank_Survival(tmp_dat$all_train)
 assign(rank_all_site[i], tmp_rank)
 tmp_file <- paste0("../../output/survival/rank/rank", i, ".csv")
 write_rank(tmp_rank, tmp_file)
}
stopCluster(cl)

rank_all_site <- paste0("rank.s", 1:n_site)
for (i in 1:n_site) {
 tmp_rank <- read.csv(paste0("../../output/survival/rank/rank", i, ".csv"))
 assign(rank_all_site[i], tmp_rank$var)
}

weight <- c()
for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "_split"))
 weight <- c(weight, nrow(tmp_dat$all_train))
}
weight <- weight/sum(weight)
rank.fed <- AutoScore_Frank(weight = weight, K = n_site, path = "../../output/survival/rank/")
var.rank <- rank.fed$var

saveRDS(var.rank, "../../output/survival/rank/fed_rank.rds")
var.rank <- readRDS("../../output/survival/rank/fed_rank.rds")

names_all_site <- paste0("homo_S", 1:n_site)
for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "_split"))
 write_cutoff(df = tmp_dat$all_train, sitename = names_all_site[i], path = "../../output/survival/cutoff.json",
              quantiles = c(0, 0.2, 0.4, 0.6, 0.8, 1))
}
all_cutoff <- read_all_cutoff(names_all_site, path = "../../output/survival/cutoff.json")

uni_cut <- get_uni_cut(names_all_site, path = "../../output/survival/cutoff.json", weight = weight)
order_uni_cut <- get_ordered_cut(var.rank, uni_cut)
saveRDS(uni_cut, "../../output/survival/uni_cut.rds")
uni_cut <- readRDS("../../output/survival/uni_cut.rds")

for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "_split"))
 write_Lbeta(names_all_site[i], tmp_dat$train, var.rank, uni_cut, fold = 5)
}

for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "_split"))
 write_Ds(names_all_site[i], names_all_site, tmp_dat$train, var.rank, uni_cut, heterogeneity = F, fold = 5)
}

for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "_split"))
 write_surr(tmp_dat$train, names_all_site[i], names_all_site, var.rank, heterogeneity = "homo", uni_cut, fold = 5)
}

subfolder_name <- "../../output/survival/fedsurv"
if (!file.exists(subfolder_name)) {
 dir.create(subfolder_name)
 cat(paste("Subfolder '", subfolder_name, "' created successfully.\n", sep = ""))
} else {
 cat(paste("Subfolder '", subfolder_name, "' already exists.\n", sep = ""))
}

get_federatedCox(data_train_all = dat1_split$train,
                site_list = names_all_site,
                var_list = var.rank,
                heterogeneity = "homo",
                uni_cut = uni_cut,
                fold = 5,
                csv_name = "../../output/survival/fedsurv/coef_global_homo_trainCV5")

auc_df <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(auc_df) <- c("AUC", "Site", "Model", "fold")
nmodel <- ncol(dat1_split$test) - 2
for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "_split"))
 tmp_df <- Get_auc_bySite(tmp_dat$train, nsite = i, var.rank, uni_cut, nmodel, nfold = 5, heterogeneity = F)
 auc_df <- rbind(auc_df, tmp_df)
}
tmp_df_avg <- auc_df %>% group_by(Model, fold) %>%
 summarize(AUC = mean(AUC), .groups = "drop") %>%
 arrange(Model, fold) %>%
 mutate(Site = "Site0")

parsimony_plot_by_site(tmp_df_avg, site_name = "Site0", output_path = "../../output/survival/fedsurv/rawPara_fed_plot_all.pdf", nfold = 5)

var_max <- 6
var_final <- var.rank[1:var_max]

subfolder_name <- sprintf("../../output/survival/v%d", var_max)
if (!file.exists(subfolder_name)) {
 dir.create(subfolder_name)
 cat(paste("Subfolder '", subfolder_name, "' created successfully.\n", sep = ""))
} else {
 cat(paste("Subfolder '", subfolder_name, "' already exists.\n", sep = ""))
}


for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "_split"))
 write_Lbeta_final(names_all_site[i], tmp_dat$all_train, var_final, uni_cut, modelIndex = "change_ref")
}

ref <- find_reference("homo", names_all_site, dat1_split$all_train, var_final, uni_cut, modelIndex = "change_ref", heterogeneity = F)
write.csv(data.frame(name = names(ref), v1 = ref), "../../output/survival/ref.csv", row.names = FALSE)

ref0 <- read.csv("../../output/survival/ref.csv")
ref <- ref0$v1
names(ref) <- ref0$name

for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "_split"))
 tmp_dat1 <- change_to_newref(tmp_dat$all_train, var_final, uni_cut, ref)
   assign(paste0("dat", i, "transformed_train_final"), tmp_dat1)
}

for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "transformed_train_final"))
 write_Lbeta_final(names_all_site[i], tmp_dat, var_final, uni_cut, modelIndex = "final")
 }

for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "transformed_train_final"))
 write_Ds_final(names_all_site[i], names_all_site, tmp_dat, var_final, uni_cut, modelIndex = "final", heterogeneity = F)
 }

for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "transformed_train_final"))
 write_surr_final(tmp_dat, names_all_site[i], names_all_site, modelIndex = "final", heterogeneity = "homo")
 }

get_federatedCox_final(dat1transformed_train_final, names_all_site, var = length(var_final), modelIndex = "final",
                        heterogeneity = "homo", i_fold="NA", csv_name = "../../output/survival/fedsurv/coef_global_homo_trainCV")

time_point <- c(1:29)

csv_name <- sprintf("../../output/survival/fedsurv/coef_global_homo_trainCVvar%d_model%s_fold%s.csv", var = length(var_final), "final", "NA")
coef <- read.csv(csv_name)
coef_vec <- coef$V1
names(coef_vec) <- coef$X

score_tbl <- show_score_table(dat1transformed_train_final, coef_vec, final_variables = var_final, max_score = 100)

for (i in 1:n_site) {
   tmp_dat <- get(paste0("dat", i, "_split"))
   tmp_dat1 <- change_to_newref(tmp_dat$test, var_final, uni_cut, ref)
   assign(paste0("dat", i, "transformed_test_final"), tmp_dat1)
}

sink(sprintf("../../output/survival/v%d/Final_fed_test_max_v%d.txt", var_max, var_max))
for (i in 1:n_site) {
 tmp_dat <- get(paste0("dat", i, "transformed_test_final"))
 print(sprintf("(homo) generating auc for Site%d validation for model fed", i))
 FedScore_survival_testing(tmp_dat, score_tbl, var_final, time_point = time_point)
}
sink()

input <- sprintf("../../output/survival/v%d/Final_fed_test_max_v%d.txt", var_max, var_max)
output <- sprintf("../../output/survival/v%d/Final_fed_testAUC_max_v%d.csv", var_max, var_max)
res <- extract_final(input, output)

