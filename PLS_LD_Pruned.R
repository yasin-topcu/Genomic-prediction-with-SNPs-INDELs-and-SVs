### =========================================================
### Partial Least Squares prediction using SKM
### LD-pruned GRMs 0.1 to 0.9
### Models M1 to M9
### =========================================================

rm(list = ls())
gc()
graphics.off()
library(plyr)
library(dplyr)
library(tidyr)
library(SKM)

setwd("C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/GRM_LD_pruning_0.1_to_0.9")

### =========================================================
### 1. Read LD-pruned Genomic Relationship Matrices
### =========================================================

grm_dir <- "C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/GRM_LD_pruning_0.1_to_0.9"

grm_files <- file.path(
  grm_dir,
  paste0("GRM_LD_", sprintf("%02d", 1:9), ".csv")
)

names(grm_files) <- paste0("M", 1:9)

G_list <- lapply(grm_files, function(f) {
  G <- read.csv(f, row.names = 1, check.names = FALSE) %>%
    as.matrix()
  
  rownames(G) <- as.character(rownames(G))
  colnames(G) <- as.character(colnames(G))
  
  return(G)
})

### =========================================================
### 2. Load and process phenotype
### =========================================================

pheno <- read.csv("C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/GRM_LD_pruning_0.1_to_0.9/pedigree_env_blups_fixedG.csv")

pheno <- pheno %>%
  dplyr::select(
    Env,
    Pedigree,
    FW = Predicted_BLUP,
    Species
  )

pheno$Species <- dplyr::recode(
  pheno$Species,
  "cerasiforme" = "SLC",
  "pimpinellipolium" = "SP",
  "lycopersicum" = "SLL"
)

common_ped_all <- Reduce(
  intersect,
  lapply(G_list, rownames)
)

pheno <- pheno %>%
  filter(Pedigree %in% common_ped_all) %>%
  arrange(Env, Pedigree) %>%
  as.data.frame()

pheno$Env <- as.factor(pheno$Env)
pheno$Pedigree <- as.factor(pheno$Pedigree)
pheno$Species <- as.factor(pheno$Species)

### =========================================================
### 3. Function to convert GRM into PLS usable features
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
### 4. Build LD-pruned genomic features
### =========================================================

X_list <- lapply(names(G_list), function(m) {
  make_features_from_grm(
    G = G_list[[m]],
    pheno = pheno,
    prefix = m
  )
})

names(X_list) <- names(G_list)

### Environment fixed effect
X_env <- model.matrix(~ 0 + Env, data = pheno)

### Optional species effect
X_spec <- model.matrix(~ 0 + Species, data = pheno)

### =========================================================
### 5. Define PLS models M1 to M9
### =========================================================

models <- lapply(names(X_list), function(m) {
  cbind(
    X_env,
    X_list[[m]]
  )
})

names(models) <- names(X_list)

### If you want to include Species as fixed effect, use this:
# models <- lapply(models, function(x) cbind(x, X_spec))

### =========================================================
### 6. Response variable
### =========================================================

y <- as.numeric(pheno$FW)

### =========================================================
### 7. Load 100 train-test splits
### =========================================================

load("C:/Users/yasin/OneDrive/Desktop/fw/Revision_1/GRM_LD_pruning_0.1_to_0.9/Split_results_100x.RData")

envs <- unique(as.character(pheno$Env))

CV1 <- list()
CV2 <- list()
CV0 <- list()
CV00 <- list()

index <- 1

safe_cor <- function(obs, pred) {
  
  good <- complete.cases(obs, pred)
  obs <- obs[good]
  pred <- pred[good]
  
  if (length(obs) < 3) return(NA_real_)
  if (sd(obs) == 0 || sd(pred) == 0) return(NA_real_)
  
  cor(obs, pred)
}

### =========================================================
### 8. PLS CV loop
### =========================================================

for (j in 1:100) {
  
  train_ped <- Split_results[[j]]$train
  test_ped  <- Split_results[[j]]$test
  
  for (e in envs) {
    
    tested_env <- setdiff(envs, e)
    
    for (m in names(models)) {
      
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
      
      pls_model <- SKM::partial_least_squares(
        X_training,
        y_training,
        method = "kernel",
        scale = FALSE,
        seed = 123,
        verbose = FALSE
      )
      
      pred <- predict(pls_model, X_testing)
      df$predicted.FW <- pred$predicted
      
      CV1[[index]] <- df %>%
        filter(Pedigree %in% train_ped & Env %in% tested_env) %>%
        group_by(Env, Species) %>%
        dplyr::summarise(
          cor = safe_cor(FW, predicted.FW),
          .groups = "drop"
        ) %>%
        mutate(model = m, iteration = j, Env_out = e)
      
      CV2[[index]] <- df %>%
        filter(Pedigree %in% test_ped & Env %in% tested_env) %>%
        group_by(Env, Species) %>%
        dplyr::summarise(
          cor = safe_cor(FW, predicted.FW),
          .groups = "drop"
        ) %>%
        mutate(model = m, iteration = j, Env_out = e)
      
      CV0[[index]] <- df %>%
        filter(Pedigree %in% train_ped & Env == e) %>%
        group_by(Env, Species) %>%
        dplyr::summarise(
          cor = safe_cor(FW, predicted.FW),
          .groups = "drop"
        ) %>%
        mutate(model = m, iteration = j, Env_out = e)
      
      CV00[[index]] <- df %>%
        filter(Pedigree %in% test_ped & Env == e) %>%
        group_by(Env, Species) %>%
        dplyr::summarise(
          cor = safe_cor(FW, predicted.FW),
          .groups = "drop"
        ) %>%
        mutate(model = m, iteration = j, Env_out = e)
      
      cat(
        ">>> Iteration:", j,
        "| Model:", m,
        "| Env out:", e,
        "| CV1 mean:", round(mean(CV1[[index]]$cor, na.rm = TRUE), 2),
        "| CV2 mean:", round(mean(CV2[[index]]$cor, na.rm = TRUE), 2),
        "| CV0 mean:", round(mean(CV0[[index]]$cor, na.rm = TRUE), 2),
        "| CV00 mean:", round(mean(CV00[[index]]$cor, na.rm = TRUE), 2),
        "\n"
      )
      
      index <- index + 1
    }
  }
}

### =========================================================
### 9. Save results
### =========================================================

CV1.out  <- plyr::ldply(CV1, data.frame)
CV2.out  <- plyr::ldply(CV2, data.frame)
CV0.out  <- plyr::ldply(CV0, data.frame)
CV00.out <- plyr::ldply(CV00, data.frame)

write.csv(CV1.out,  "PLS_LD_CV1.out.csv",  row.names = FALSE)
write.csv(CV2.out,  "PLS_LD_CV2.out.csv",  row.names = FALSE)
write.csv(CV0.out,  "PLS_LD_CV0.out.csv",  row.names = FALSE)
write.csv(CV00.out, "PLS_LD_CV00.out.csv", row.names = FALSE)

### =========================================================
### 10. Overall summary
### =========================================================

summary_all <- bind_rows(
  CV1.out  %>% mutate(CV = "CV1"),
  CV2.out  %>% mutate(CV = "CV2"),
  CV0.out  %>% mutate(CV = "CV0"),
  CV00.out %>% mutate(CV = "CV00")
)

summary_mean <- summary_all %>%
  group_by(CV, model, Species) %>%
  dplyr::summarise(
    mean_cor = mean(cor, na.rm = TRUE),
    sd_cor = sd(cor, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_mean, "PLS_LD_summary_mean.csv", row.names = FALSE)

print(summary_mean)
