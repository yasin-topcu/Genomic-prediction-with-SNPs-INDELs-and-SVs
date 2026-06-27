## =========================================================
## BGLR prediction
## Parallel computing
## SNP, INDEL, SV, Combined, Combined + GxE
## With cor, RMSE, MAE, and Bias
## =========================================================

rm(list = ls())
gc()
graphics.off()

library(plyr)
library(dplyr)
library(tidyr)
library(BGLR)
library(fastmatrix)
library(ggplot2)
library(foreach)
library(doParallel)

setwd("C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/RMSE")

## =========================================================
## 1) Read genomic relationship matrices
## =========================================================

G1 <- read.csv("GRM.GAPIT.csv", row.names = 1) %>% as.matrix()
G2 <- read.csv("GRM_indel.csv", row.names = 1) %>% as.matrix()
G3 <- read.csv("GRM_sv.csv",    row.names = 1) %>% as.matrix()

## =========================================================
## 2) Load and process phenotype
## =========================================================

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
  pivot_wider(
    names_from = Env,
    values_from = FW
  ) %>%
  arrange(Pedigree) %>%
  as.data.frame()

pheno <- pheno[pheno$Pedigree %in% rownames(G1), ]

pheno <- pheno %>%
  pivot_longer(
    cols = -c(Pedigree, Species),
    names_to = "Env",
    values_to = "FW"
  ) %>%
  arrange(Env, Pedigree) %>%
  as.data.frame()

pheno$Env <- as.factor(pheno$Env)
pheno$Pedigree <- as.factor(pheno$Pedigree)
pheno$Species <- as.factor(pheno$Species)

## =========================================================
## 3) Build design matrices
## =========================================================

ZE    <- model.matrix(~ pheno$Env - 1)
ZPed  <- model.matrix(~ pheno$Pedigree - 1)
ZSpec <- model.matrix(~ pheno$Species - 1)

## =========================================================
## 4) Construct kernel matrices
## =========================================================

K1 <- ZPed %*% G1 %*% t(ZPed)
K2 <- ZPed %*% G2 %*% t(ZPed)
K3 <- ZPed %*% G3 %*% t(ZPed)

ZEZE <- tcrossprod(ZE)

K11 <- K1 * ZEZE
K22 <- K2 * ZEZE
K33 <- K3 * ZEZE

## =========================================================
## 5) Eigendecomposition
## =========================================================

Z_K1.eig  <- eigen(K1)
Z_K2.eig  <- eigen(K2)
Z_K3.eig  <- eigen(K3)

Z_K11.eig <- eigen(K11)
Z_K22.eig <- eigen(K22)
Z_K33.eig <- eigen(K33)

## =========================================================
## 6) Define multi-kernel BGLR models
## =========================================================

M1 <- list(
  ENV = list(X = ZE, model = "BRR"),
  SNP = list(V = Z_K1.eig$vectors, d = Z_K1.eig$values, model = "RKHS")
)

M2 <- list(
  ENV = list(X = ZE, model = "BRR"),
  INDEL = list(V = Z_K2.eig$vectors, d = Z_K2.eig$values, model = "RKHS")
)

M3 <- list(
  ENV = list(X = ZE, model = "BRR"),
  SV = list(V = Z_K3.eig$vectors, d = Z_K3.eig$values, model = "RKHS")
)

M4 <- list(
  ENV = list(X = ZE, model = "BRR"),
  SNP = list(V = Z_K1.eig$vectors, d = Z_K1.eig$values, model = "RKHS"),
  INDEL = list(V = Z_K2.eig$vectors, d = Z_K2.eig$values, model = "RKHS"),
  SV = list(V = Z_K3.eig$vectors, d = Z_K3.eig$values, model = "RKHS")
)

M5 <- list(
  ENV = list(X = ZE, model = "BRR"),
  SNP = list(V = Z_K1.eig$vectors, d = Z_K1.eig$values, model = "RKHS"),
  SNPxENV = list(V = Z_K11.eig$vectors, d = Z_K11.eig$values, model = "RKHS"),
  INDEL = list(V = Z_K2.eig$vectors, d = Z_K2.eig$values, model = "RKHS"),
  INDELxENV = list(V = Z_K22.eig$vectors, d = Z_K22.eig$values, model = "RKHS"),
  SV = list(V = Z_K3.eig$vectors, d = Z_K3.eig$values, model = "RKHS"),
  SVxENV = list(V = Z_K33.eig$vectors, d = Z_K33.eig$values, model = "RKHS")
)

models <- list(
  M1 = M1,
  M2 = M2,
  M3 = M3,
  M4 = M4,
  M5 = M5
)

## =========================================================
## 7) Load 100 training-testing splits
## =========================================================

load("C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/RMSE/Split_results_100x.RData")

envs <- unique(pheno$Env) %>% as.character()

## =========================================================
## 8) Safe metric functions
## =========================================================

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

## =========================================================
## 9) Parallel BGLR cross-validation loop
## =========================================================

n_cores <- parallel::detectCores() - 1

## If RAM is limited, use this instead:
# n_cores <- 4

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
  .packages = c("dplyr", "BGLR"),
  .combine = rbind
) %dopar% {
  
  j <- task_grid$iteration[task]
  e <- task_grid$Env_out[task]
  m <- task_grid$model[task]
  
  train_ped <- Split_results[[j]]$train
  test_ped  <- Split_results[[j]]$test
  
  tested_env <- setdiff(envs, e)
  
  df <- pheno
  df$Y2 <- df$FW
  
  ## Remove FW values for test genotypes and environment-out
  df$Y2[df$Env == e | df$Pedigree %in% test_ped] <- NA
  
  ## Important for parallel BGLR:
  ## unique saveAt prefix prevents workers overwriting each other's files
  save_prefix <- paste0(
    tempdir(),
    "/BGLR_iter_", j,
    "_env_", e,
    "_model_", m,
    "_task_", task,
    "_"
  )
  
  set.seed(123 + task)
  
  fit <- BGLR(
    y = df$Y2,
    ETA = models[[m]],
    nIter = 5000,
    burnIn = 1000,
    verbose = FALSE,
    saveAt = save_prefix
  )
  
  df$predicted.FW <- fit$yHat
  
  CV1_temp <- df %>%
    dplyr::filter(Pedigree %in% train_ped & Env %in% tested_env) %>%
    dplyr::group_by(Env, Species) %>%
    dplyr::summarise(
      safe_metrics(FW, predicted.FW),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      CV = "CV1",
      model = m,
      iteration = j,
      Env_out = e
    )
  
  CV2_temp <- df %>%
    dplyr::filter(Pedigree %in% test_ped & Env %in% tested_env) %>%
    dplyr::group_by(Env, Species) %>%
    dplyr::summarise(
      safe_metrics(FW, predicted.FW),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      CV = "CV2",
      model = m,
      iteration = j,
      Env_out = e
    )
  
  CV0_temp <- df %>%
    dplyr::filter(Pedigree %in% train_ped & Env == e) %>%
    dplyr::group_by(Env, Species) %>%
    dplyr::summarise(
      safe_metrics(FW, predicted.FW),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      CV = "CV0",
      model = m,
      iteration = j,
      Env_out = e
    )
  
  CV00_temp <- df %>%
    dplyr::filter(Pedigree %in% test_ped & Env == e) %>%
    dplyr::group_by(Env, Species) %>%
    dplyr::summarise(
      safe_metrics(FW, predicted.FW),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      CV = "CV00",
      model = m,
      iteration = j,
      Env_out = e
    )
  
  dplyr::bind_rows(
    CV1_temp,
    CV2_temp,
    CV0_temp,
    CV00_temp
  )
}

stopCluster(cl)
registerDoSEQ()

## =========================================================
## 10) Save CV results
## =========================================================

summary_all <- parallel_results

CV1.out  <- summary_all %>% dplyr::filter(CV == "CV1")
CV2.out  <- summary_all %>% dplyr::filter(CV == "CV2")
CV0.out  <- summary_all %>% dplyr::filter(CV == "CV0")
CV00.out <- summary_all %>% dplyr::filter(CV == "CV00")

write.csv(CV1.out,  "BGLR_CV1_metrics_each_iteration.csv",  row.names = FALSE)
write.csv(CV2.out,  "BGLR_CV2_metrics_each_iteration.csv",  row.names = FALSE)
write.csv(CV0.out,  "BGLR_CV0_metrics_each_iteration.csv",  row.names = FALSE)
write.csv(CV00.out, "BGLR_CV00_metrics_each_iteration.csv", row.names = FALSE)

write.csv(summary_all, "BGLR_all_CV_metrics_each_iteration.csv", row.names = FALSE)

## =========================================================
## 11) Summarize CV results
## =========================================================

cv_summary <- summary_all %>%
  dplyr::group_by(CV, model, Species) %>%
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

write.csv(cv_summary, "BGLR_summary_mean_cor_RMSE_MAE_Bias.csv", row.names = FALSE)

print(cv_summary)