#--------------------------------------------------
# Functions for random forest and comparative analyses
#--------------------------------------------------
# Prepare input data (Filter data, combine X and Y variables into one data.frame)
  # df: table containing readthrough efficiency and RPKM information. Final output from rt_efficiency.R
  # cds_rpkm_cutoff: transcripts with RPKM of the CDS less than the specified number will be discarded
  # utr3_rpkm_cutoff: transcripts with RPKM of the utr3_region less than the specified number will be discarded
  # group_param: for classification; a vector taking in 
  #               1) a mode to group transcripts to High and Low group by 
  #                   1. standard deviation ("std")
  #                   2. percentile ("percentile")
  #                   3. no grouping ("none")
  #               2) number associated with mode choice
  # feature_file: data table containing features (X variables) associated with each transcript
  # utr3_region: choices are full 3'-UTR length ("utr3"), extension ("ext"), or fixed part of 3'-UTR ("utr3p")
  # logbase: base of log transform readthrough efficiency. Default is log2.
  # genes_to_analyze: list of genes to include in analysis. Default is NULL, meaning all genes in df are included
prep <- function(df, cds_rpkm_cutoff = 5, utr3_rpkm_cutoff = 0.5, group_param = c("std", 1), feature_file, 
                 utr3_region = "ext", logbase = 2, genes_to_analyze = NULL) {
  require(data.table)
  require(dplyr)
  # 1-2. Keep mRNAs with 1. UTR annotations and 2. non-overlapping sequence
  if (length(genes_to_analyze) > 0) {
    message("Subset transcripts")
    df <- df[which(df$transcript %in% genes_to_analyze), ]
  }
  # 3. Filter out un-analyzable readthrough
  df <- df[df$rpkm_cds > cds_rpkm_cutoff & df[[paste0("rpkm_", utr3_region)]] > utr3_rpkm_cutoff, ]
  # 4. log10-transformed readthrough efficiency
  df$log_rte <- log(x = df[[paste0("rte_", utr3_region)]], base = logbase)
  # 5. Break data into high and low readthrough efficiency groups
  if (group_param[1] == "std") { # 5.1 Using standard deviation
    message(paste0("Group genes using standard deviation of ", group_param[2]))
    df$z <- data.table(scale(df$log_rte, center = TRUE, scale = TRUE))
    std <- as.integer(group_param[2])
    df <- data.table(df[z < -std | z > std, ])
    df$Group <- df$z > std
    df$Group <- gsub("FALSE", "Low", df$Group)
    df$Group <- gsub("TRUE", "High", df$Group)
  } else if (group_param[1] == "percentile") {  # 5.2 Using percentile
    message(paste0("Group genes using top and bottom ", group_param[2]), "th percentile")
    df$percentile <- data.table(dplyr::percent_rank(df$log_rte))
    pt <- as.integer(group_param[2])/100
    df <- data.table(df[percentile < pt | percentile > 1-pt])
    df$Group <- df$percentile > pt
    df$Group <- gsub("FALSE", "Low", df$Group)
    df$Group <- gsub("TRUE", "High", df$Group)
  } else if (group_param[1] == "none") {  # 5.3 No grouping
    message("Data will not be grouped")
  } else {
    message("Invalid grouping parameter")
  }
  # 6. Combine data with feature_file
  col_df <- colnames(df)
  col_ff <- colnames(feature_file)
  common_cols <- intersect(col_df, col_ff)
  keep <- feature_file[, c("transcript", col_ff[!(col_ff %in% common_cols)]), with = FALSE]
  df <- left_join(df, keep, by = "transcript")
  # 7. Convert any character variables to factor
  df <- df %>% mutate_if(is.character, as.factor)
  df <- data.table(df)
  return(df)
}
#--------------------------------------------------
# Normalized RMSE for regression (https://www.marinedatascience.co/blog/2019/01/07/normalizing-the-rmse/)
nrmse_calc <-  function(mse, obs) {
  rmse <- sqrt(mse)
  nrmse <- data.frame(row.names = 1, rmse = rmse,                                  # not normalized
                      nrmse_sd = rmse/sd(obs),                                     # normalized by standard deviation
                      nrmse_mean = rmse/mean(obs),                                 # normalized by mean
                      nrmse_range = rmse/(max(obs) - min(obs)),                    # normalized by range
                      nrmse_iq = rmse/(quantile(obs, 0.75) - quantile(obs, 0.25))) # normalized by interquartile range
  return(nrmse)
}
#--------------------------------------------------
# Random forest regression
rp_reg <- function(df, y_var = "log_rte", ntree = 100, cv_fold = 5, nrep = 1000) {
  require(caret)
  require(randomForest)
  require(rfPermute)
  message("\nStart")
  names(df)[names(df) == y_var] <- "y_var"
  folds <- createFolds(y = df$y_var, k = cv_fold)
  rp_imp <- list()
  rp_nmrse <- list()
  for (i in 1:cv_fold) {
    message("Fold ", i)
    train_df <- df[-folds[[i]], ]
    test_df <- df[folds[[i]], ]
    rp <- rfPermute(y_var ~ ., train_df, ntree = ntree, num.rep = nrep)
    imp <- data.frame(importance(rp))
    imp$feature <- rownames(imp)
    pred <- predict(rp, test_df)
    nmrse_df <- bind_rows(list(training = nrmse_calc(mse = rp$rf$mse[length(rp$rf$mse)], obs = train_df$log_rte),
                               testing = nrmse_calc(mse = mean((test_df$log_rte - pred)^2), obs = test_df$log_rte)), .id = "dataset")
    rp_imp[[i]] <- imp
    rp_nmrse[[i]] <- nmrse_df
  }
  rp_res <- list()
  rp_res[["imp"]] <- bind_rows(rp_imp, .id = "fold")
  rp_res[["metrics"]] <- bind_rows(rp_nmrse, .id = "fold")
  return(rp_res)
}
#--------------------------------------------------
# Random forest classification
rp_class <- function(df, y_var = "Group", ntree = 100, cv_fold = 5, nrep = 1000) {
  require(caret)
  require(randomForest)
  require(rfPermute)
  message("\nStart")
  names(df)[names(df) == y_var] <- "Group"
  df$Group <- as.factor(df$Group)
  df <- df %>% mutate_if(is.character, as.factor)
  folds <- createFolds(y = df$Group, k = cv_fold)
  rp_imp <- list()
  rp_auc <- list()
  for (i in 1:cv_fold) {
    message("Fold ", i)
    train_df <- df[-folds[[i]], ]
    test_df <- df[folds[[i]], ]
    rp <- rfPermute(Group ~ ., train_df, ntree = ntree, num.rep = nrep)
    imp <- data.frame(importance(rp))
    imp$feature <- rownames(imp)
    pred <- predict(rp, test_df, type = "prob")
    auc_df <- data.frame(dataset = c("traning", "testing"), 
                         auc = c(auc(roc(train_df$Group, rp$rf$votes[, 2])), auc(roc(test_df$Group, pred[, 2]))))
    rp_imp[[i]] <- imp
    rp_auc[[i]] <- auc_df
  }
  rp_res <- list()
  rp_res[["imp"]] <- bind_rows(rp_imp, .id = "fold")
  rp_res[["metrics"]] <- bind_rows(rp_auc, .id = "fold")
  return(rp_res)
}
#--------------------------------------------------
# Identify significant features and average importance scores across repeated CVs
ave_importance <- function(data_list, rf_type = "class", sig_mode = "across", pval_cutoff = 0.05, num_cutoff = 0.8) {
  # data_list: list of data.frames of "imp" from different repeats of rp_class or rp_reg
  # rf_type: type of random forest. "reg" for regression or "class" for classification.
  # sig_mode: mode of identifying significant features: "across" or "individual"
  # pval_cutoff: p-value cutoff for significant features
  # num_cutoff: frequency of the feature being significant across models to be considered finally significant
  # Example: data_list = 5 repeats of 5-fold cross-validation, sig_mode = "across", pval_cutoff = 0.05, num_cutoff = 0.8 
  #          The parameters would make sure that features that have p-value < 0.05 in at least 20 out of 25 models (0.8) are considered significant
  # Example: data_list = 5 repeats of 5-fold cross-validation, sig_mode = "individual", pval_cutoff = 0.05, num_cutoff = 0.8 
  #          The parameters would make sure that features that have p-value < 0.05 in at least 4 out of 5 models (0.8) in all 5 repeats are considered significant
  if (rf_type == "class") {
    data_list <- lapply(data_list, function(x) {
      x <- x[, c("fold", "feature", "MeanDecreaseAccuracy", "MeanDecreaseAccuracy.pval")]
      colnames(x) <- c("fold", "feature", "imp", "imp_pval")
      return(x)})
  } else if (rf_type == "reg") {
    data_list <- lapply(data_list, function(x) {
      x <- x[, c("fold", "feature", "X.IncMSE", "X.IncMSE.pval")]
      colnames(x) <- c("fold", "feature", "imp", "imp_pval")
      return(x)})
  }
  find_sig <- function(df, pval_cutoff, num_cutoff) {
    freqdf <- data.frame(table(df[df$imp_pval < pval_cutoff, ]$feature))  # Frequency of a feature being significant
    xx <- nrow(df) / length(unique(df$feature)) # Number of folds
    cn <- round(num_cutoff * xx, digits = 0)
    message(cn, " out of ", xx)
    sigf <- freqdf[which(freqdf$Freq >= cn), ]$Var1
    return(sigf)
  }
  if (sig_mode == "across") {
    df <- bind_rows(data_list, .id = "repeats")
    sigf <- find_sig(df, pval_cutoff = pval_cutoff, num_cutoff = num_cutoff)
    print(sigf)
  } else if (sig_mode == "individual") {
    siglist <- lapply(data_list, find_sig, pval_cutoff = pval_cutoff, num_cutoff = num_cutoff)
    print(siglist)
    sigf <- Reduce(intersect, siglist)
    print(sigf)
    df <- bind_rows(data_list, .id = "repeats")
  } else {
    message("Invalid mode")
  }
  avedf <- data.frame(df %>% group_by(feature) %>% summarize(imp_ave = mean(imp), imp_sd = sd(imp)))
  avedf$sig <- ifelse(test = avedf$feature %in% sigf, yes = paste0("p < ", pval_cutoff), no = "ns")
  return(avedf)
}
#--------------------------------------------------
# Compare readthrough efficiency between mRNA groups, divided into groups by categorical X variable
compare_rt_categorical <- function(dat, vars, y_var = "log_rte") {
  compare_rt <- function(dat, var, y_var) {
    df <- dat[, c(y_var, var), with = FALSE]
    names(df) <- c("y_var", "Var")
    cm <- ggpubr::compare_means(formula = y_var~Var, data = df, method = "wilcox.test", ref.group = ".all.", p.adjust.method = "BH")
    tb <- boxplot(y_var~Var, data = df, plot = FALSE)
    tbs <- data.frame(tb$stats)
    colnames(tbs) <- tb$names
    tbs <- data.frame(t(tbs))
    tbs$group2 <- row.names(tbs)
    res <- left_join(cm[, c("group2", "p", "p.adj")], tbs[, c("group2", "X3")], by = "group2")
    res$samp_median <- median(df$y_var)
    names(res) <- c("Var", "p", "p.adj", "median", "samp_median")
    res$median_diff <- res$median - res$samp_median
    res <- within(res, {
      p_sig = ifelse(p.adj < 0.05, "p < 0.05", "ns")
      p_hw = ifelse(p.adj < 0.05, 0.9, 0.6)
    })
    res$feature <- var
    res$p_sig <- factor(res$p_sig, levels = c("p < 0.05", "ns"))
    return(res)
  }
  res_list <- list()
  for (x in vars) {
    res_list[[x]] <- compare_rt(dat, x, y_var)
  }
  res_df <- data.table(bind_rows(res_list))
  return(res_df)
}
#--------------------------------------------------
