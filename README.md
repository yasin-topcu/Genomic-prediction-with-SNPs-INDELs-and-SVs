Multi-Environment Genomic Prediction of Tomato Fruit Weight Using SNPs, INDELs, and Structural Variants
Overview

This repository contains the data and R scripts used for the genomic prediction analyses presented in the manuscript:

"Multi-Variant GWAS and Genomic Prediction Dissect the Genetic Architecture Underlying Tomato Fruit Weight"

The repository includes phenotype data, genomic relationship matrices (GRMs), structural variant data, cross-validation splits, and prediction scripts implementing multiple genomic prediction algorithms. Model performance was evaluated using prediction accuracy, root mean square error (RMSE), mean absolute error (MAE), and prediction bias.

####Repository Contents

Phenotypic Data
1- Blues.Pheno.corrected.csv	Corrected phenotypic BLUPs for fruit weight used for genomic prediction.
2- pedigree_env_blups_fixedG.csv	Multi-environment BLUPs used in genomic prediction analyses.

Genomic Relationship Matrices
1- GRM.GAPIT.csv	SNP-based genomic relationship matrix.
2- GRM_indel.csv	INDEL-based genomic relationship matrix.
3- GRM_sv.csv	Structural variant (SV)-based genomic relationship matrix.

Structural Variant Genotypes
1- Hapmap.SVs.hmp.zip	HapMap-formatted structural variant genotype dataset used for SV-based analyses.

Cross-Validation Files
1- geno_70_30_splits_100x.RData	One hundred independent 70% training / 30% testing population splits used for model evaluation.
2- Split_results_100x.RData	Saved cross-validation partitions and intermediate results used across prediction methods.

Prediction Scripts
1- BGLR_with_RMSE_MAE_BIAS.R	Bayesian Genomic Linear Regression (BGLR).
2- PLS_with_RMSE_MAE_BIAS.R	Partial Least Squares (PLS) genomic prediction.
3- PLS_LD_Pruned.R	PLS genomic prediction using LD-pruned marker sets.
4- RF_with_RMSE_MAE_BIAS.R	Random Forest genomic prediction.
5- DL_with_RMSE_MAE_BIAS.R	Deep Learning genomic prediction.

LD-pruned PLS script
PLS_LD_Pruned.R
GRM_LD_01.csv
GRM_LD_02.csv
GRM_LD_03.csv
GRM_LD_04.csv
GRM_LD_05.csv
GRM_LD_06.csv
GRM_LD_07.csv
GRM_LD_08.csv
GRM_LD_09.csv
GRM_LD_01.csv to GRM_LD_09.csv LD-pruned genomic relationship matrices generated using LD thresholds from 0.1 to 0.9.

Five genomic prediction models were evaluated:Model	Description

M1	Environment + SNPs
M2	Environment + INDELs
M3	Environment + Structural Variants
M4	Environment + SNPs + INDELs + Structural Variants
M5	M4 + Marker × Environment interaction kernels

Cross-Validation Strategy
Genomic prediction was evaluated using repeated random cross-validation.
Training population: 70%
Validation population: 30%
Number of replicates: 100 independent random splits

Prediction performance was assessed under four scenarios:Scenario	Description

CV1	Known genotypes evaluated in known environments
CV2	New genotypes evaluated in known environments
CV0	Known genotypes evaluated in new environments
CV00	New genotypes evaluated in new environments

The identical train-test partitions were used across all prediction algorithms to ensure fair model comparison.

Model Evaluation:Prediction performance was assessed using:

Pearson correlation between observed and predicted phenotypes (prediction accuracy)
Root Mean Square Error (RMSE)
Mean Absolute Error (MAE)
Prediction bias (regression slope of observed on predicted values)

Software Requirements: The analyses were performed in R (version 4.5.3, 4.6.0, 4.5.2, 4.3.3, 4.1.2).

Major R packages include:

BGLR
pls
randomForest
keras / tensorflow
data.table
tidyverse
caret
Matrix
fastmatrix
ggplot2
plyr
dplyr
tidyr
SKM
foreach
doParallel


If you use these data or scripts, please cite:

Topcu, Y., Adak, A., Sapkota, M., et al. Multi-Variant GWAS and Genomic Prediction Dissect the Genetic Architecture Underlying Tomato Fruit Weight. (Manuscript submitted).
