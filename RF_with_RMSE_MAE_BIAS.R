### =========================================================
### Random Forest prediction using SKM
### Parallel computing
### Fixed RF parameters: 500 trees, node size = 5
### SNP, INDEL, SV, Combined, Combined + GxE
### With cor, RMSE, MAE, and Bias
### =========================================================

rm(list = ls())
gc()
graphics.off()

library(plyr)
library(dplyr)
library(tidyr)
library(SKM)
library(foreach)
library(doParallel)

setwd("C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/RMSE")

### =========================================================
### 1. Read Genomic Relationship Matrices
### =========================================================

G1 <- read.csv("GRM.GAPIT.csv", row.names = 1) %>% as.matrix()
G2 <- read.csv("GRM_indel.csv", row.names = 1) %>% as.matrix()
G3 <- read.csv("GRM_sv.csv",    row.names = 1) %>% as.matrix()

### =========================================================
### 2. Load and process phenotype
### =========================================================

pheno <- read.csv("C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/RMSE/pedigree_env_blups_fixedG.csv")

pheno <- pheno %>%
  dplyr::select(
    Env,
    Pedigree,
    FW = Predicted_BLUP,
    Species
  )

pheno$Species <- recode(
  pheno$Species,
  "cerasiforme" = "SLC",
  "pimpinellipolium" = "SP",
  "lycopersicum" = "SLL"
)

pheno <- pheno %>%
  filter(Pedigree %in% rownames(G1)) %>%
  arrange(Env, Pedigree) %>%
  as.data.frame()

pheno$Env <- as.factor(pheno$Env)
pheno$Pedigree <- as.factor(pheno$Pedigree)
pheno$Species <- as.factor(pheno$Species)

### =========================================================
### 3. Function to convert GRM into RF usable features
### =========================================================

make_features_from_grm <- function(G, pheno, prefix = "G") {
  
  common_ped <- intersect(as.character(pheno$Pedigree), rownames(G))
  common_ped <- unique(common_ped)
  
  G <- G[common_ped, common_ped]
  
  eig <- eigen(G)
  
  values <- eig$values
  values[values < 0] <- 0
  
  scores <- eig$vectors %*% diag(sqrt(values))
  rownames(scores) <- rownames(G)
  
  Z <- model.matrix(~ 0 + Pedigree, data = pheno)
  colnames(Z) <- gsub("Pedigree", "", colnames(Z))
  
  scores <- scores[colnames(Z), , drop = FALSE]
  
  Xg <- Z %*% scores
  colnames(Xg) <- paste0(prefix, "_PC", seq_len(ncol(Xg)))
  
  return(Xg)
}

### =========================================================
### 4. Build genomic features
### =========================================================

X_snp   <- make_features_from_grm(G1, pheno, "SNP")
X_indel <- make_features_from_grm(G2, pheno, "INDEL")
X_sv    <- make_features_from_grm(G3, pheno, "SV")

X_env  <- model.matrix(~ 0 + Env, data = pheno)
X_spec <- model.matrix(~ 0 + Species, data = pheno)

### =========================================================
### 5. Function to create G x E features
### =========================================================

make_gxe <- function(Xg, Xenv, prefix = "GxE") {
  
  out <- list()
  
  for (e in colnames(Xenv)) {
    temp <- Xg * as.numeric(Xenv[, e])
    colnames(temp) <- paste0(prefix, "_", e, "_", colnames(Xg))
    out[[e]] <- temp
  }
  
  do.call(cbind, out)
}

X_snp_env   <- make_gxe(X_snp,   X_env, "SNPxENV")
X_indel_env <- make_gxe(X_indel, X_env, "INDELxENV")
X_sv_env    <- make_gxe(X_sv,    X_env, "SVxENV")

### =========================================================
### 6. Define RF models
### =========================================================

models <- list(
  
  M1 = cbind(
    X_env,
    X_snp
  ),
  
  M2 = cbind(
    X_env,
    X_indel
  ),
  
  M3 = cbind(
    X_env,
    X_sv
  ),
  
  M4 = cbind(
    X_env,
    X_snp,
    X_indel,
    X_sv
  ),
  
  M5 = cbind(
    X_env,
    X_snp,
    X_snp_env,
    X_indel,
    X_indel_env,
    X_sv,
    X_sv_env
  )
)

### Optional species fixed effect
# models <- lapply(models, function(x) cbind(x, X_spec))

### =========================================================
### 7. Response variable
### =========================================================

y <- as.numeric(pheno$FW)

### =========================================================
### 8. Load 100 train-test splits
### =========================================================

load("C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/RMSE/Split_results_100x.RData")

envs <- unique(as.character(pheno$Env))

### =========================================================
### 9. Safe metric functions
### =========================================================

safe_cor <- function(obs, pred) {
  
  good <- complete.cases(obs, pred)
  obs <- obs[good]
  pred <- pred[good]
  
  if (length(obs) < 3) return(NA_real_)
  if (sd(obs) == 0 || sd(pred) == 0) return(NA_real_)
  
  cor(obs, pred)
}

safe_metrics <- function(obs, pred) {
  
  good <- complete.cases(obs, pred)
  obs <- obs[good]
  pred <- pred[good]
  
  if (length(obs) < 3) {
    return(data.frame(
      cor  = NA_real_,
      RMSE = NA_real_,
      MAE  = NA_real_,
      Bias = NA_real_
    ))
  }
  
  data.frame(
    cor  = safe_cor(obs, pred),
    RMSE = sqrt(mean((pred - obs)^2)),
    MAE  = mean(abs(pred - obs)),
    Bias = mean(pred - obs)
  )
}

### =========================================================
### 10. Parallel Random Forest CV loop
### =========================================================

n_cores <- parallel::detectCores() - 1

cl <- makeCluster(n_cores)
registerDoParallel(cl)

task_grid <- expand.grid(
  iteration = 1:100,
  Env_out = envs,
  model = names(models),
  stringsAsFactors = FALSE
)

parallel_results <- foreach(
  task = 1:nrow(task_grid),
  .packages = c("dplyr", "SKM"),
  .combine = rbind
) %dopar% {
  
  j <- task_grid$iteration[task]
  e <- task_grid$Env_out[task]
  m <- task_grid$model[task]
  
  train_ped <- Split_results[[j]]$train
  test_ped  <- Split_results[[j]]$test
  
  tested_env <- setdiff(envs, e)
  
  X <- models[[m]]
  df <- pheno
  df$predicted.FW <- NA
  
  train_index <- which(
    df$Pedigree %in% train_ped &
      df$Env != e &
      !is.na(df$FW)
  )
  
  test_index <- seq_len(nrow(df))
  
  X_training <- X[train_index, , drop = FALSE]
  y_training <- y[train_index]
  X_testing  <- X[test_index, , drop = FALSE]
  
  rf_model <- SKM::random_forest(
    X_training,
    y_training,
    trees_number = 500,
    node_size = 5
  )
  
  pred <- predict(rf_model, X_testing)
  df$predicted.FW <- pred$predicted
  
  CV1_temp <- df %>%
    filter(Pedigree %in% train_ped & Env %in% tested_env) %>%
    group_by(Env, Species) %>%
    dplyr::summarise(
      safe_metrics(FW, predicted.FW),
      .groups = "drop"
    ) %>%
    mutate(
      CV = "CV1",
      model = m,
      iteration = j,
      Env_out = e
    )
  
  CV2_temp <- df %>%
    filter(Pedigree %in% test_ped & Env %in% tested_env) %>%
    group_by(Env, Species) %>%
    dplyr::summarise(
      safe_metrics(FW, predicted.FW),
      .groups = "drop"
    ) %>%
    mutate(
      CV = "CV2",
      model = m,
      iteration = j,
      Env_out = e
    )
  
  CV0_temp <- df %>%
    filter(Pedigree %in% train_ped & Env == e) %>%
    group_by(Env, Species) %>%
    dplyr::summarise(
      safe_metrics(FW, predicted.FW),
      .groups = "drop"
    ) %>%
    mutate(
      CV = "CV0",
      model = m,
      iteration = j,
      Env_out = e
    )
  
  CV00_temp <- df %>%
    filter(Pedigree %in% test_ped & Env == e) %>%
    group_by(Env, Species) %>%
    dplyr::summarise(
      safe_metrics(FW, predicted.FW),
      .groups = "drop"
    ) %>%
    mutate(
      CV = "CV00",
      model = m,
      iteration = j,
      Env_out = e
    )
  
  bind_rows(
    CV1_temp,
    CV2_temp,
    CV0_temp,
    CV00_temp
  )
}

stopCluster(cl)
registerDoSEQ()

### =========================================================
### 11. Save CV results for each iteration
### =========================================================

summary_all <- parallel_results

CV1.out  <- summary_all %>% filter(CV == "CV1")
CV2.out  <- summary_all %>% filter(CV == "CV2")
CV0.out  <- summary_all %>% filter(CV == "CV0")
CV00.out <- summary_all %>% filter(CV == "CV00")

write.csv(CV1.out,  "RF_CV1_metrics_each_iteration.csv",  row.names = FALSE)
write.csv(CV2.out,  "RF_CV2_metrics_each_iteration.csv",  row.names = FALSE)
write.csv(CV0.out,  "RF_CV0_metrics_each_iteration.csv",  row.names = FALSE)
write.csv(CV00.out, "RF_CV00_metrics_each_iteration.csv", row.names = FALSE)

write.csv(summary_all, "RF_all_CV_metrics_each_iteration.csv", row.names = FALSE)

### =========================================================
### 12. Overall summary
### =========================================================

summary_mean <- summary_all %>%
  group_by(CV, model, Species) %>%
  dplyr::summarise(
    mean_cor  = mean(cor, na.rm = TRUE),
    sd_cor    = sd(cor, na.rm = TRUE),
    
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_RMSE   = sd(RMSE, na.rm = TRUE),
    
    mean_MAE  = mean(MAE, na.rm = TRUE),
    sd_MAE    = sd(MAE, na.rm = TRUE),
    
    mean_Bias = mean(Bias, na.rm = TRUE),
    sd_Bias   = sd(Bias, na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(summary_mean, "RF_summary_mean_cor_RMSE_MAE_Bias.csv", row.names = FALSE)

print(summary_mean)
