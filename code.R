options(repos = c(CRAN = "https://cloud.r-project.org"))
options(mc.cores = parallel::detectCores())
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages <- c(
  "readxl", "dplyr", "tidyr", "survival", "randomForestSRC", "glmnet",
  "RColorBrewer", "plsRcox", "superpc", "gbm", "survivalsvm", "CoxBoost",
  "ComplexHeatmap", "doParallel", "foreach", "Boruta", "timeROC", "rmda",
  "ggplot2", "grid", "gridExtra", "pROC", "data.table", "writexl"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("ComplexHeatmap", "survcomp", "mixOmics")) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}
library(parallel)

total_cores <- detectCores()
n_core <- max(1, total_cores - 6)
registerDoParallel(cores = n_core)

data_raw <- fread("your_data.csv", header = TRUE, stringsAsFactors = FALSE)
message("Data reading completed: ", nrow(data_raw), " rows × ", ncol(data_raw), " columns")

data_raw[, WHtR := WC / Height]
data_raw[, TG_HDL := TG / HDL]
data_raw[, CMI := TG_HDL * WHtR]

if (!"CKM" %in% colnames(data_raw)) stop("CKM staging column is missing in the data")

time_col <- "Follow_up"
event_main <- "Mortality_Status"
event_cvd <- "Death_CVD"
event_renal <- "Death_Renal"
event_ckm_high <- "CKM_high"
data_raw[, CKM_high := ifelse(CKM >= 3, 1, 0)]

weight_col <- "cycle_weight"
strata_col <- "SDMVSTRA"
cluster_col <- "SDMVPSU"

exclude_cols <- c("SEQN", "Year", "SDMVSTRA", "SDMVPSU", "cycle_weight",
                  "Mortality_Status", "Follow_up", "Death_CVD", "Death_Renal",
                  "CKD", "Cardiovascular_disease", "Cancer", "CKM", "CKM_high",
                  "Depression", "PHENOAGE", "KDM", "CMI")

feature_cols <- setdiff(colnames(data_raw), exclude_cols)
feature_cols <- unique(c(feature_cols, "WHtR", "TG_HDL"))
feature_cols <- setdiff(feature_cols, "CMI")

list(feature_cols)

feature_data_raw <- as.data.frame(data_raw[, feature_cols, with = FALSE])
surv_data <- as.data.frame(data_raw[, c(time_col, event_main, event_cvd, event_renal, 
                                        event_ckm_high, weight_col, strata_col, 
                                        cluster_col, "CKM"), with = FALSE])

cat_vars <- c("Sex", "Smoking", "Education", "Drinking", "Physical_Activity",
              "Ethnicity", "Marital_Status", "PIR", "Hypertension", "Diabetes",
              "Depression", "Cardiovascular_disease", "Cancer")
cat_vars <- intersect(cat_vars, feature_cols)
for (v in cat_vars) {
  feature_data_raw[[v]] <- as.factor(feature_data_raw[[v]])
}

original_feature <- feature_data_raw
original_surv <- surv_data

# ------------------------------ 2. Auxiliary functions------------------------------
quiet <- function(...) {
  sink(tempfile())
  on.exit(sink())
  eval(...)
}

prepare_event <- function(x, event_value = 1) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x <- ifelse(x > 0, 1, 0)
  x
}

preprocess_features <- function(feature_df, surv_df, cat_vars = NULL) {
  if (!is.null(cat_vars) && length(cat_vars) > 0) {
    num_vars <- setdiff(colnames(feature_df), cat_vars)
    if (length(num_vars) > 0) {
      num_mat <- as.matrix(feature_df[, num_vars, drop = FALSE])
      mode(num_mat) <- "numeric"
    } else {
      num_mat <- matrix(nrow = nrow(feature_df), ncol = 0)
    }
    cat_mat <- model.matrix(~ . - 1, data = feature_df[, cat_vars, drop = FALSE])
    expr_mat <- cbind(num_mat, cat_mat)
  } else {
    expr_mat <- as.matrix(feature_df)
    mode(expr_mat) <- "numeric"
  }
  list(expr = expr_mat, surv = surv_df)
}

scaleData_fast <- function(data, train_mean, train_sd) {
  sweep(sweep(data, 2, train_mean, "-"), 2, train_sd, "/")
}

# ------------------------------ 3. Model training functions ------------------------------
RunML_generic <- function(method, Train_expr, Train_surv, mode = "Model",
                          timeVar, statusVar, family = "cox") {
  method <- gsub(" ", "", method)
  method_name <- gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param <- gsub("(\\w+)\\[(.+)\\]", "\\2", method)
  method_param <- switch(method_name,
                         "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
                         "StepCox" = list("direction" = method_param), NULL)
  args <- list("Train_expr" = Train_expr, "Train_surv" = Train_surv, "mode" = mode,
               "timeVar" = timeVar, "statusVar" = statusVar, "family" = family)
  args <- c(args, method_param)
  do.call(what = paste0("Run", method_name), args = args)
}

RunEnet <- function(Train_expr, Train_surv, mode, timeVar, statusVar, alpha, family = "cox") {
  if (family == "cox") {
    y <- Surv(Train_surv[[timeVar]], Train_surv[[statusVar]])
    family_glmnet <- "cox"
  } else {
    y <- as.factor(Train_surv[[statusVar]])
    family_glmnet <- "binomial"
  }
  cv.fit <- tryCatch(
    cv.glmnet(x = Train_expr, y = y, family = family_glmnet, alpha = alpha, nfolds = 3,
              trace.it = 0, maxit = 100000, thresh = 1e-7),
    error = function(e) NULL)
  if (is.null(cv.fit)) stop("cv.glmnet 失败")
  lambda <- cv.fit$lambda.1se
  fit <- glmnet(x = Train_expr, y = y, family = family_glmnet, alpha = alpha,
                lambda = lambda, maxit = 100000, thresh = 1e-7)
  if (all(coef(fit) == 0) && cv.fit$lambda.min < lambda) {
    fit <- glmnet(x = Train_expr, y = y, family = family_glmnet, alpha = alpha,
                  lambda = cv.fit$lambda.min, maxit = 100000, thresh = 1e-7)
  }
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}
RunLasso <- function(...) RunEnet(..., alpha = 1)
RunRidge <- function(...) RunEnet(..., alpha = 0)

RunCoxPH <- function(Train_expr, Train_surv, mode, timeVar, statusVar, family = "cox") {
  df <- as.data.frame(Train_expr)
  df$time <- Train_surv[[timeVar]]
  df$event <- Train_surv[[statusVar]]
  fit <- coxph(Surv(time, event) ~ ., data = df)
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunStepCox <- function(Train_expr, Train_surv, mode, timeVar, statusVar, direction, family = "cox") {
  df <- as.data.frame(Train_expr)
  df$time <- Train_surv[[timeVar]]
  df$event <- Train_surv[[statusVar]]
  fit <- step(coxph(Surv(time, event) ~ ., data = df), direction = direction, trace = 0)
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunsurvivalSVM <- function(Train_expr, Train_surv, mode, timeVar, statusVar, family = "cox") {
  df <- as.data.frame(Train_expr)
  df$time <- Train_surv[[timeVar]]
  df$event <- Train_surv[[statusVar]]
  fit <- survivalsvm(Surv(time, event) ~ ., data = df, gamma.mu = 1, opt.meth = "ipop")
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunCoxBoost <- function(Train_expr, Train_surv, mode, timeVar, statusVar, family = "cox") {
  pen <- optimCoxBoostPenalty(time = Train_surv[[timeVar]], status = Train_surv[[statusVar]],
                              x = Train_expr, trace = FALSE, start.penalty = 500)
  cv.res <- cv.CoxBoost(time = Train_surv[[timeVar]], status = Train_surv[[statusVar]],
                        x = Train_expr, maxstepno = 50, K = 3, type = "verweij", penalty = pen$penalty)
  fit <- CoxBoost(time = Train_surv[[timeVar]], status = Train_surv[[statusVar]],
                  x = Train_expr, stepno = cv.res$optimal.step, penalty = pen$penalty)
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSuperPC <- function(Train_expr, Train_surv, mode, timeVar, statusVar, family = "cox") {
  data <- list(x = t(scale(Train_expr)), y = Train_surv[[timeVar]],
               censoring.status = Train_surv[[statusVar]], featurenames = colnames(Train_expr))
  fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5)
  cv.fit <- suppressWarnings(superpc.cv(fit, data, n.threshold = 5, n.fold = 3,
                                        n.components = 2, min.features = 5,
                                        max.features = nrow(data$x)))
  fit$threshold <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])]
  fit$data <- data
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunplsRcox <- function(Train_expr, Train_surv, mode, timeVar, statusVar, family = "cox") {
  data <- list(x = Train_expr, time = Train_surv[[timeVar]], status = Train_surv[[statusVar]])
  cv.plsRcox.res <- cv.plsRcox(data = data, nt = 3, verbose = FALSE)
  fit <- plsRcox(Xplan = data$x, time = data$time, event = data$status,
                 nt = as.numeric(cv.plsRcox.res[5]), verbose = FALSE, sparse = TRUE)
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunRSF <- function(Train_expr, Train_surv, mode, timeVar, statusVar, family = "cox") {
  df <- as.data.frame(Train_expr)
  df$time <- Train_surv[[timeVar]]
  df$event <- Train_surv[[statusVar]]
  fit <- rfsrc(Surv(time, event) ~ ., data = df, ntree = 300, nodesize = 5,
               splitrule = 'logrank', importance = TRUE, forest = TRUE)
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunGBM <- function(Train_expr, Train_surv, mode, timeVar, statusVar, family = "cox") {
  df <- as.data.frame(Train_expr)
  df$time <- Train_surv[[timeVar]]
  df$event <- Train_surv[[statusVar]]
  fit_cv <- gbm(Surv(time, event) ~ ., data = df, distribution = 'coxph',
                n.trees = 2000, interaction.depth = 3, n.minobsinnode = 10,
                shrinkage = 0.01, cv.folds = 3)
  best <- which.min(fit_cv$cv.error)
  fit <- gbm(Surv(time, event) ~ ., data = df, distribution = 'coxph',
             n.trees = best, interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.01)
  fit$subFeature <- colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

ExtractVar <- function(fit) {
  quiet(switch(class(fit)[1],
               "coxnet" = rownames(coef(fit))[which(coef(fit)[, 1] != 0)],
               "coxph" = names(coef(fit)),
               "survivalsvm" = fit$var.names,
               "CoxBoost" = names(coef(fit))[abs(coef(fit)) > 0],
               "superpc" = names(fit$feature.scores)[abs(fit$feature.scores) > 0.5],
               "plsRcoxmodel" = rownames(fit$Coeffs)[fit$Coeffs != 0],
               "rfsrc" = var.select(fit, verbose = FALSE)$topvars,
               "gbm" = rownames(summary.gbm(fit, plotit = FALSE))[summary.gbm(fit, plotit = FALSE)$rel.inf > 0]))
}

CalRiskScore <- function(fit, new_data) {
  common_features <- intersect(fit$subFeature, colnames(new_data))
  if (length(common_features) < length(fit$subFeature)) {
    new_data_full <- matrix(0, nrow = nrow(new_data), ncol = length(fit$subFeature))
    colnames(new_data_full) <- fit$subFeature
    idx <- match(common_features, fit$subFeature)
    new_data_full[, idx] <- as.matrix(new_data[, common_features])
    new_data <- new_data_full
  } else {
    new_data <- new_data[, fit$subFeature, drop = FALSE]
  }
  RS <- quiet(switch(class(fit)[1],
                     "coxnet" = predict(fit, type = 'link', as.matrix(new_data)),
                     "coxph" = predict(fit, type = 'lp', as.data.frame(new_data)),
                     "survivalsvm" = predict(fit, as.data.frame(new_data))$predicted,
                     "CoxBoost" = predict(fit, type = "lp", as.data.frame(new_data)),
                     "superpc" = superpc.predict(fit, fit$data, list(x = t(scale(as.matrix(new_data)))),
                                                 fit$threshold, 1)$v.pred,
                     "plsRcoxmodel" = predict(fit, type = "lp", as.data.frame(new_data)),
                     "rfsrc" = predict(fit, as.data.frame(new_data))$predicted,
                     "gbm" = predict(fit, type = 'link', as.data.frame(new_data))))
  as.vector(RS)
}

RunEval <- function(fit, Test_expr, Test_surv, timeVar, statusVar) {
  new_data <- Test_expr[, fit$subFeature, drop = FALSE]
  RS <- CalRiskScore(fit, new_data)
  eval_data <- data.frame(RS = RS,
                          time = Test_surv[[timeVar]],
                          event = Test_surv[[statusVar]])
  res <- unname(summary(coxph(Surv(time, event) ~ RS, data = eval_data))$concordance["C"])
  return(res)
}

CalAUC <- function(fit, new_data, new_surv, statusVar) {
  RS <- CalRiskScore(fit, new_data)
  if (all(is.na(RS))) return(NA)
  roc_obj <- tryCatch(pROC::roc(new_surv[[statusVar]], RS, quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(NA)
  as.numeric(roc_obj$auc)
}

# ------------------------------ 4. Configuration parameters ------------------------------
# Whether to skip stability selection (set to TRUE if results have been saved)
skip_stability_selection <- FALSE   # If TRUE, load previous stability selection results from file

# Define screening methods (for stability selection)
filter_methods <- c(
  "simple", "Lasso", "Ridge", "Enet[alpha=0.4]", "Enet[alpha=0.6]", "Enet[alpha=0.8]",
  "StepCox[both]", "StepCox[backward]", "StepCox[forward]",
  "CoxBoost", "SuperPC", "plsRcox", "RSF", "GBM"
)

# Define modeling methods (for the second stage)
model_methods <- c(
  "CoxPH", "Lasso", "Ridge", "Enet[alpha=0.4]", "Enet[alpha=0.6]", "Enet[alpha=0.8]",
  "StepCox[both]", "StepCox[backward]", "StepCox[forward]",
  "CoxBoost", "SuperPC", "plsRcox", "survivalSVM", "RSF", "GBM"
)

# ------------------------------ 5. Data preprocessing and splitting ------------------------------
prep <- preprocess_features(feature_data_raw, surv_data, cat_vars = cat_vars)
expr_raw <- prep$expr
surv_raw <- prep$surv

set.seed(123)
n <- nrow(expr_raw)
train_idx <- sample(1:n, size = round(0.7 * n))
test_idx <- setdiff(1:n, train_idx)

message("Training set size: ", length(train_idx), ", test set size: ", length(test_idx))

train_expr_raw <- expr_raw[train_idx, ]
train_surv <- surv_raw[train_idx, ]
test_expr_raw <- expr_raw[test_idx, ]
test_surv <- surv_raw[test_idx, ]

# Standardization (using training set mean and standard deviation)
train_mean <- apply(train_expr_raw, 2, mean, na.rm = TRUE)
train_sd <- apply(train_expr_raw, 2, sd, na.rm = TRUE)
train_expr <- scaleData_fast(train_expr_raw, train_mean, train_sd)
test_expr <- scaleData_fast(test_expr_raw, train_mean, train_sd)

# Remove features with zero variance
zero_var <- apply(train_expr, 2, function(x) var(x, na.rm = TRUE) == 0)
if (any(zero_var)) {
  message("Removing features with zero variance: ", paste(names(zero_var)[zero_var], collapse = ", "))
  train_expr <- train_expr[, !zero_var, drop = FALSE]
  test_expr <- test_expr[, !zero_var, drop = FALSE]
}

# Save standardization parameters for later use
std_params <- list(mean = train_mean, sd = train_sd)
saveRDS(std_params, file = "standard_params.rds")

# ------------------------------ 6. Stability screening and saving ------------------------------
stability_results_file <- "stability_results.rds"

if (!skip_stability_selection) {
  # Stability selection function (parallel)
  stability_selection_parallel <- function(expr, surv, timeVar, statusVar, family,
                                           filter_methods, n_bootstrap = 50, threshold = 0.6,
                                           n_cores = n_core) {
    n <- nrow(expr)
    set.seed(123)
    boot_indices <- replicate(n_bootstrap, sample(1:n, replace = TRUE), simplify = FALSE)
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    exports <- c("RunML_generic", "RunEnet", "RunLasso", "RunRidge", "RunCoxPH", "RunStepCox",
                 "RunsurvivalSVM", "RunCoxBoost", "RunSuperPC", "RunplsRcox", "RunRSF", "RunGBM",
                 "ExtractVar", "quiet")
    clusterExport(cl, exports, envir = environment())
    
    boot_res <- foreach(i = 1:n_bootstrap, .combine = function(x, y) {
      if (is.null(x)) return(y)
      if (is.null(y)) return(x)
      merged <- list()
      all_fms <- unique(c(names(x), names(y)))
      for (fm in all_fms) {
        x_vec <- if (is.null(x[[fm]])) character(0) else x[[fm]]
        y_vec <- if (is.null(y[[fm]])) character(0) else y[[fm]]
        merged[[fm]] <- c(x_vec, y_vec)
      }
      merged
    }, .packages = c("glmnet", "survival", "randomForestSRC", "CoxBoost", "plsRcox", "superpc", "gbm", "survivalsvm"),
    .export = exports) %dopar% {
      boot_data <- expr[boot_indices[[i]], , drop = FALSE]
      boot_surv <- surv[boot_indices[[i]], , drop = FALSE]
      selected_list <- list()
      for (fm in filter_methods) {
        tryCatch({
          if (fm == "simple") {
            vars <- colnames(boot_data)
          } else {
            vars <- RunML_generic(method = fm,
                                  Train_expr = boot_data,
                                  Train_surv = boot_surv,
                                  mode = "Variable",
                                  timeVar = timeVar,
                                  statusVar = statusVar,
                                  family = family)
            if (is.null(vars)) vars <- character(0)
          }
          selected_list[[fm]] <- vars
        }, error = function(e) {
          selected_list[[fm]] <<- character(0)
        })
      }
      selected_list
    }
    stopCluster(cl)
    
    all_features <- colnames(expr)
    freq <- setNames(rep(0, length(all_features)), all_features)
    for (fm in names(boot_res)) {
      selected_combined <- boot_res[[fm]]
      if (length(selected_combined) > 0) {
        tab <- table(selected_combined)
        freq[names(tab)] <- freq[names(tab)] + as.vector(tab)
      }
    }
    total_trials <- n_bootstrap * length(filter_methods)
    freq_norm <- freq / total_trials
    stable_vars <- names(freq_norm[freq_norm >= threshold])
    return(list(vars = stable_vars, freq = freq_norm))
  }
  
  message("\n=== Stability screening: Target A (all-cause mortality) ===")
  stab_A <- stability_selection_parallel(
    expr = train_expr, surv = train_surv,
    timeVar = time_col, statusVar = event_main,
    family = "cox",
    filter_methods = filter_methods,
    n_bootstrap = 50, threshold = 0.6,
    n_cores = n_core
  )
  vars_A <- stab_A$vars
  
  message("\n=== Stability screening: Target B (CKM high-risk) ===")
  stab_B <- stability_selection_parallel(
    expr = train_expr, surv = train_surv,
    timeVar = time_col, statusVar = event_ckm_high,
    family = "binomial",
    filter_methods = filter_methods,
    n_bootstrap = 50, threshold = 0.6,
    n_cores = n_core
  )
  vars_B <- stab_B$vars

  
  final_features <- unique(c(vars_A, vars_B))
  message("\nFinal candidate feature pool size: ", length(final_features))
  
  if (length(final_features) < 3) stop("Candidate features too few for subsequent modeling")
  
  train_expr_final <- train_expr[, final_features, drop = FALSE]
  test_expr_final <- test_expr[, final_features, drop = FALSE]
  
  # Save stability selection results
  stability_result <- list(
    stab_A = stab_A,
    stab_B = stab_B,
    final_features = final_features,
    train_expr_final = train_expr_final,
    test_expr_final = test_expr_final,
    train_surv = train_surv,
    test_surv = test_surv,
    train_expr_raw = train_expr_raw,
    test_expr_raw = test_expr_raw,
    cat_vars = cat_vars,
    train_idx = train_idx,
    test_idx = test_idx,
    time_col = time_col,
    event_main = event_main,
    event_cvd = event_cvd,
    event_renal = event_renal,
    event_ckm_high = event_ckm_high
  )
  saveRDS(stability_result, file = stability_results_file)
  message("\nStability screening results saved to: ", stability_results_file)
  
} else {
  # Load previously saved stability screening results
  if (!file.exists(stability_results_file)) {
    stop("File ", stability_results_file, " does not exist. Cannot skip stability screening.")
  }
  stability_result <- readRDS(stability_results_file)
  final_features <- stability_result$final_features
  train_expr_final <- stability_result$train_expr_final
  test_expr_final <- stability_result$test_expr_final
  train_surv <- stability_result$train_surv
  test_surv <- stability_result$test_surv
  cat_vars <- stability_result$cat_vars
  train_idx <- stability_result$train_idx
  test_idx <- stability_result$test_idx
  time_col <- stability_result$time_col
  event_main <- stability_result$event_main
  event_cvd <- stability_result$event_cvd
  event_renal <- stability_result$event_renal
  event_ckm_high <- stability_result$event_ckm_high
  message("\nLoaded stability screening results from: ", stability_results_file)
  message("Final candidate feature pool size: ", length(final_features))
}

# ------------------------------ 7. Multi-model training and evaluation ------------------------------
model_training_file <- "model_training_results.rds"
# If the file exists and the user wishes to skip training (optional), directly load and plot the heatmap, otherwise retrain
skip_model_training <- FALSE   # If TRUE, load training results from file without retraining

if (!skip_model_training) {
  message("\nStarting multi-model training (total ", length(filter_methods) * length(model_methods), " combinations)")
  
  combos <- expand.grid(filter = filter_methods, model = model_methods, stringsAsFactors = FALSE)
  combos <- unique(combos)
  model_names <- paste0(combos$filter, "+", combos$model)
  n_models <- length(model_names)
  
  pb <- txtProgressBar(min = 0, max = n_models, style = 3)
  
  results_list <- list()
  cindex_mat <- matrix(NA, nrow = n_models, ncol = 6)
  colnames(cindex_mat) <- c("Train_C", "Test_C", "Test_CVD", "Test_Renal", "Test_CKM_AUC", "Avg_C")
  rownames(cindex_mat) <- model_names
  
  for (i in 1:n_models) {
    filter_m <- combos$filter[i]
    model_m <- combos$model[i]
    combo_name <- model_names[i]
    
    tryCatch({
      if (filter_m == "simple") {
        selected_vars <- final_features
      } else {
        selected_vars <- RunML_generic(method = filter_m,
                                       Train_expr = train_expr_final,
                                       Train_surv = train_surv,
                                       mode = "Variable",
                                       timeVar = time_col,
                                       statusVar = event_main,
                                       family = "cox")
        selected_vars <- intersect(selected_vars, final_features)
      }
      
      if (length(selected_vars) < 2) {
        next
      }
      
      train_expr_sub <- train_expr_final[, selected_vars, drop = FALSE]
      test_expr_sub <- test_expr_final[, selected_vars, drop = FALSE]
      
      fit <- RunML_generic(method = model_m,
                           Train_expr = train_expr_sub,
                           Train_surv = train_surv,
                           mode = "Model",
                           timeVar = time_col,
                           statusVar = event_main,
                           family = "cox")
      
      c_train <- RunEval(fit, Test_expr = train_expr_sub, Test_surv = train_surv, timeVar = time_col, statusVar = event_main)
      c_test <- RunEval(fit, Test_expr = test_expr_sub, Test_surv = test_surv, timeVar = time_col, statusVar = event_main)
      c_cvd <- RunEval(fit, Test_expr = test_expr_sub, Test_surv = test_surv, timeVar = time_col, statusVar = event_cvd)
      c_renal <- RunEval(fit, Test_expr = test_expr_sub, Test_surv = test_surv, timeVar = time_col, statusVar = event_renal)
      auc_ckm <- CalAUC(fit, new_data = test_expr_sub, new_surv = test_surv, statusVar = event_ckm_high)
      avg_c <- mean(c(c_test, c_cvd, c_renal, auc_ckm), na.rm = TRUE)
      
      cindex_mat[i, ] <- c(c_train, c_test, c_cvd, c_renal, auc_ckm, avg_c)
      results_list[[combo_name]] <- list(fit = fit, selected_vars = selected_vars)
      
    }, error = function(e) {
      # Ignore errors
    })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Remove invalid rows
  valid_idx <- rowSums(is.na(cindex_mat[, c("Test_C", "Test_CKM_AUC")])) < 2
  cindex_mat <- cindex_mat[valid_idx, , drop = FALSE]
  results_list <- results_list[valid_idx]
  
  # Save model training results
  model_training_result <- list(
    cindex_mat = cindex_mat,
    results_list = results_list,
    combos = combos,
    filter_methods = filter_methods,
    model_methods = model_methods
  )
  saveRDS(model_training_result, file = model_training_file)
  message("\nModel training results saved to: ", model_training_file)
  
} else {
  if (!file.exists(model_training_file)) {
    stop("File ", model_training_file, " does not exist. Cannot skip model training.")
  }
  model_training_result <- readRDS(model_training_file)
  cindex_mat <- model_training_result$cindex_mat
  results_list <- model_training_result$results_list
  combos <- model_training_result$combos
  message("\nLoaded model training results from: ", model_training_file)
}

# ------------------------------ 8. Pareto optimal selection ------------------------------
plot_data <- data.frame(
  model = rownames(cindex_mat),
  c_index = cindex_mat[, "Test_C"],
  auc_ckm = cindex_mat[, "Test_CKM_AUC"]
)
plot_data <- plot_data[complete.cases(plot_data), ]

pareto_front <- function(x, y) {
  n <- length(x)
  is_pareto <- rep(TRUE, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j && x[j] >= x[i] && y[j] >= y[i] && (x[j] > x[i] || y[j] > y[i])) {
        is_pareto[i] <- FALSE
        break
      }
    }
  }
  is_pareto
}

is_pareto <- pareto_front(plot_data$c_index, plot_data$auc_ckm)
pareto_models <- plot_data[is_pareto, ]
pareto_models$sum <- pareto_models$c_index + pareto_models$auc_ckm
best_model_name <- pareto_models$model[which.max(pareto_models$sum)]
message("\nPareto optimal model: ", best_model_name)

best_fit <- results_list[[best_model_name]]$fit
best_vars <- results_list[[best_model_name]]$selected_vars
best_idx <- which(rownames(cindex_mat) == best_model_name)
best_cindex <- cindex_mat[best_idx, "Test_C"]
best_auc <- cindex_mat[best_idx, "Test_CKM_AUC"]
cat("\nBest model: ", best_model_name, "\n")
cat("Test set C-index: ", round(best_cindex, 3), "\n")
cat("Test set CKM AUC: ", round(best_auc, 3), "\n")

# ------------------------------ 9. Extract weights, compute CKMAI, and save results ------------------------------
# Extract feature weights
if (inherits(best_fit, c("coxnet", "coxph", "glmnet", "lognet"))) {
  if (inherits(best_fit, "coxnet") || inherits(best_fit, "lognet")) {
    coef_vec <- as.matrix(coef(best_fit))[, 1]
  } else {
    coef_vec <- coef(best_fit)
  }
  coef_vec <- coef_vec[best_vars]
  coef_vec[is.na(coef_vec)] <- 0
  weights <- abs(coef_vec)
} else if (inherits(best_fit, "rfsrc")) {
  vimp <- best_fit$importance
  weights <- abs(vimp[best_vars])
  weights[is.na(weights)] <- 0
} else if (inherits(best_fit, "gbm")) {
  imp <- summary(best_fit, plotit = FALSE)
  weights <- rep(0, length(best_vars))
  names(weights) <- best_vars
  idx <- match(imp$var, best_vars)
  weights[na.omit(idx)] <- imp$rel.inf[!is.na(idx)]
} else {
  weights <- rep(1, length(best_vars))
  names(weights) <- best_vars
}
weights_norm <- weights / max(weights, na.rm = TRUE)

# Compute CKMAI
train_expr_final <- train_expr[, best_vars, drop = FALSE]
test_expr_final <- test_expr[, best_vars, drop = FALSE]

train_score_raw <- as.vector(train_expr_final %*% weights_norm)
train_min <- min(train_score_raw)
train_max <- max(train_score_raw)

calc_CKMAI <- function(expr_matrix, vars, weights_norm, train_min, train_max) {
  score <- as.vector(expr_matrix[, vars, drop = FALSE] %*% weights_norm)
  score_01 <- (score - train_min) / (train_max - train_min)
  return(score_01 * 100)
}

train_CKMAI <- calc_CKMAI(train_expr_final, best_vars, weights_norm, train_min, train_max)
test_CKMAI <- calc_CKMAI(test_expr_final, best_vars, weights_norm, train_min, train_max)

# Save Excel files
library(writexl)

train_original_feature <- original_feature[train_idx, , drop = FALSE]
train_original_surv <- original_surv[train_idx, , drop = FALSE]
train_data <- cbind(train_original_feature, train_original_surv)
train_data$RowID <- rownames(train_data)
write_xlsx(list(TrainingSet = train_data), path = "CKMAI_Training_Data.xlsx")

test_original_feature <- original_feature[test_idx, , drop = FALSE]
test_original_surv <- original_surv[test_idx, , drop = FALSE]
test_data <- cbind(test_original_feature, test_original_surv)
test_data$RowID <- rownames(test_data)
write_xlsx(list(TestSet = test_data), path = "CKMAI_Test_Data.xlsx")

model_info <- data.frame(
  Model = best_model_name,
  Feature = best_vars,
  Weight = weights_norm,
  AbsoluteWeight = weights,
  stringsAsFactors = FALSE
)
eval_metrics <- as.data.frame(cindex_mat)
eval_metrics$Model <- rownames(eval_metrics)
write_xlsx(list(ModelWeights = model_info, EvaluationMetrics = eval_metrics), 
           path = "CKMAI_Model_Summary.xlsx")

ckmai_data <- data.frame(
  RowID = rownames(train_data),
  CKMAI = train_CKMAI,
  Set = "Train"
)
ckmai_data_test <- data.frame(
  RowID = rownames(test_data),
  CKMAI = test_CKMAI,
  Set = "Test"
)
ckmai_all <- rbind(ckmai_data, ckmai_data_test)
write_xlsx(list(CKMAI_Scores = ckmai_all), path = "CKMAI_Scores.xlsx")

# ------------------------------ 10. Plot heatmap ------------------------------
SimpleHeatmap <- function(Cindex_mat = NULL,
                          avg_Cindex = NULL,
                          CohortCol = NULL,
                          barCol = NULL,
                          col = c("#4195C1", "#FFFFFF", "#CB5746"),
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns = FALSE, cluster_rows = FALSE) {
  if (is.null(dim(Cindex_mat))) {
    Cindex_mat <- t(as.matrix(Cindex_mat))
    rownames(Cindex_mat) <- "Model"
  }
  col_ha <- columnAnnotation("Cohort" = colnames(Cindex_mat),
                             col = list("Cohort" = CohortCol),
                             show_annotation_name = FALSE)
  if (!is.null(avg_Cindex)) {
    row_ha <- rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                               gp = gpar(fill = barCol, col = NA),
                                               add_numbers = TRUE, numbers_offset = unit(-10, "mm"),
                                               axis_param = list("labels_rot" = 0),
                                               numbers_gp = gpar(fontsize = 9, col = "white"),
                                               width = unit(3, "cm")),
                            show_annotation_name = FALSE)
  } else {
    row_ha <- NULL
  }
  hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
                right_annotation = row_ha,
                top_annotation = col_ha,
                col = col,
                rect_gp = gpar(col = "black", lwd = 1),
                cluster_columns = cluster_columns, cluster_rows = cluster_rows,
                show_column_names = FALSE,
                show_row_names = TRUE,
                row_names_side = "left",
                width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
                height = unit(cellheight * nrow(Cindex_mat), "cm"),
                column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
                column_title = NULL,
                cell_fun = function(j, i, x, y, w, h, col) {
                  grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                            x, y, gp = gpar(fontsize = 10))
                })
  return(hm)
}

heat_mat <- cindex_mat[, c("Test_C", "Test_CVD", "Test_Renal", "Test_CKM_AUC"), drop = FALSE]
colnames(heat_mat) <- c("All-cause", "CVD", "Renal", "CKM AUC")
complete_rows <- complete.cases(heat_mat[, c("All-cause", "CKM AUC")])
heat_mat <- heat_mat[complete_rows, , drop = FALSE]
avg_c <- rowMeans(heat_mat, na.rm = TRUE)
cohort_col <- c("All-cause" = "#E41A1C", "CVD" = "#377EB8", "Renal" = "#4DAF4A", "CKM AUC" = "#984EA3")
pdf("model_performance_train_test.pdf", width = 12, height = max(6, nrow(heat_mat)*0.3))
hm <- SimpleHeatmap(Cindex_mat = heat_mat,
                    avg_Cindex = avg_c,
                    CohortCol = cohort_col,
                    barCol = "#FF7F00",
                    col = c("#053061", "#F7F7F7", "#B2182B"),
                    cellwidth = 1.2, cellheight = 0.5,
                    cluster_columns = FALSE, cluster_rows = TRUE)
draw(hm)
dev.off()
message("Heatmap saved as model_performance_train_test.pdf")

# Final output
cat("\n========== FINAL RESULTS ==========\n")
cat("Best model combination: ", best_model_name, "\n")
cat("Number of selected features: ", length(best_vars), "\n")
cat("Training set CKMAI range: ", round(train_min, 3), " - ", round(train_max, 3), "\n")
cat("Test set C-index for all-cause mortality: ", round(best_cindex, 3), "\n")
cat("Test set AUC for CKM high-risk: ", round(best_auc, 3), "\n")

message("\nAll results have been saved as Excel files:")
message(" - CKMAI_Training_Data.xlsx")
message(" - CKMAI_Test_Data.xlsx")
message(" - CKMAI_Model_Summary.xlsx")
message(" - CKMAI_Scores.xlsx")
message("Heatmap saved as model_performance_train_test.pdf")

# Save final results for subsequent analyses
final_result <- list(
  best_fit = best_fit,
  best_name = best_model_name,
  best_vars = best_vars,
  weights = weights_norm,
  train_min = train_min,
  train_max = train_max,
  train_mean = train_mean,
  train_sd = train_sd,
  train_idx = train_idx,
  test_idx = test_idx,
  train_CKMAI = train_CKMAI,
  test_CKMAI = test_CKMAI,
  cindex_mat = cindex_mat
)
saveRDS(final_result, file = "CKMAI_results.rds")

write.csv(data.frame(
  RowID = rownames(train_data),
  CKMAI = train_CKMAI,
  Set = "Train"
), file = "CKMAI_train_scores.csv", row.names = FALSE)
write.csv(data.frame(
  RowID = rownames(test_data),
  CKMAI = test_CKMAI,
  Set = "Test"
), file = "CKMAI_test_scores.csv", row.names = FALSE)

message("Part 1 results saved to CKMAI_results.rds and standard_params.rds")