stopifnot(file.exists("O2MCA.Rproj"))
required <- c("glmnet", "digest", "FactoMineR")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install packages first: ", paste(missing, collapse = ", "))
repo <- normalizePath(".")

helpers <- new.env(parent = globalenv())
source_file <- if (file.exists("05_sensitivities_and_internal_test.R")) {
  "05_sensitivities_and_internal_test.R"
} else {
  "functions.R"
}
for (e in parse(source_file)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, helpers)
}
needed <- c("scale_fit", "scale_apply", "feature_matrix", "project_drug", "project_omics")
stopifnot(all(vapply(needed, exists, logical(1), envir = helpers, inherits = FALSE)))

protocol <- readRDS("prepared_inputs/protocol_config.rds")
omics <- readRDS("prepared_inputs/omics_verified_headers.rds")
fp <- as.matrix(read.csv("prepared_inputs/PubchemFP.csv", row.names = 1, check.names = FALSE))
rownames(fp) <- trimws(rownames(fp))

a <- protocol$cell_pair_split
pairs <- protocol$analysis_pairs
stopifnot(identical(a$pair_id, pairs$pair_id))
dev <- pairs[a$split == "development", ]
folds <- a$inner_fold[a$split == "development"]
stopifnot(length(unique(folds)) == 5L)

cv_spec <- readRDS("results/02_cv/run_specification.rds")
stopifnot(cv_spec$rank == 20L, !cv_spec$single_reference,
          identical(cv_spec$pair_ids, dev$pair_id), identical(cv_spec$folds, folds))
for (name in names(cv_spec$inputs)) {
  path <- file.path("prepared_inputs", basename(name))
  stopifnot(file.exists(path),
            digest::digest(file = path, algo = "sha256") == cv_spec$inputs[[name]])
}

out <- tempfile("15_linear_baselines_", tmpdir = "results")
dir.create(out, recursive = TRUE)
wc <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)

alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)
lambda_grid <- 10^seq(3, -5, length.out = 81)

writeLines(c(
  "Post hoc reviewer-requested comparison after original test outcomes were inspected.",
  "O2PLS-MCA only: the same 20 omics scores plus 20 MCA coordinates as the primary XGBoost route.",
  "Original fold-specific transformations reused. All scaling and tuning use development cells only.",
  "OLS has no tuning. Ridge alpha=0; elastic-net family alpha in {0.25,0.5,0.75,1}, including lasso.",
  "Fixed lambda grid: 10^seq(3,-5,length.out=81). Minimize pooled OOF RMSE.",
  "Ties within 1e-12 prefer larger lambda, then smaller alpha. No test-based retuning.",
  "This is an additional descriptive analysis, not new independent confirmatory evidence."
), file.path(out, "ANALYSIS_PLAN.txt"))
saveRDS(list(alpha = alpha_grid, lambda = lambda_grid, input_hashes = cv_spec$inputs,
             created = Sys.time()), file.path(out, "specification.rds"))

load_checked <- function(path) {
  stopifnot(file.exists(path), file.exists(paste0(path, ".sha256")))
  stopifnot(digest::digest(file = path, algo = "sha256") ==
              trimws(readLines(paste0(path, ".sha256"), n = 1)))
  readRDS(path)
}

features <- function(tr, va, shared_path, model_path) {
  s <- load_checked(shared_path)
  m <- load_checked(model_path)
  tc <- sort(unique(tr$cell_name))
  vc <- sort(unique(va$cell_name))
  td <- sort(unique(tr$drug_key))
  vd <- sort(unique(va$drug_key))
  stopifnot(identical(tc, s$train_cells), identical(td, s$train_drugs),
            !any(tc %in% vc), all(vd %in% td), m$method == "O2PLS")
  scores <- function(cells) helpers$project_omics(
    m,
    helpers$scale_apply(omics$expression[cells, , drop = FALSE], s$expression_scaler),
    helpers$scale_apply(omics$methylation[cells, , drop = FALSE], s$methylation_scaler)
  )
  x <- helpers$feature_matrix(tr, helpers$project_drug(s, fp, td), scores(tc), 1:20)
  v <- helpers$feature_matrix(va, helpers$project_drug(s, fp, vd), scores(vc), 1:20)
  sc <- helpers$scale_fit(x, FALSE)
  x <- helpers$scale_apply(x, sc)
  v <- helpers$scale_apply(v, sc)
  stopifnot(ncol(x) == 40L, identical(colnames(x), colnames(v)),
            all(is.finite(x)), all(is.finite(v)))
  list(x = x, v = v, scaler = sc)
}

ols_predict <- function(x, y, v) {
  fit <- lm.fit(cbind(Intercept = 1, x), y)
  if (fit$rank != ncol(x) + 1L) stop("OLS design is rank deficient; inspect before proceeding.")
  list(fit = fit, prediction = drop(cbind(Intercept = 1, v) %*% fit$coefficients))
}

sse <- matrix(0, length(alpha_grid), length(lambda_grid))
ols_sse <- 0
for (k in 1:5) {
  message("Development fold ", k, "/5")
  tr <- dev[folds != k, ]
  va <- dev[folds == k, ]
  z <- features(tr, va,
                sprintf("results/02_cv/fold%d_shared.rds", k),
                sprintf("results/02_cv/fold%d_O2PLS_transform.rds", k))
  ols <- ols_predict(z$x, tr$response_value, z$v)
  ols_sse <- ols_sse + sum((ols$prediction - va$response_value)^2)
  for (j in seq_along(alpha_grid)) {
    fit <- glmnet::glmnet(z$x, tr$response_value, alpha = alpha_grid[j], lambda = lambda_grid,
                          standardize = FALSE, intercept = TRUE, thresh = 1e-10, maxit = 1000000)
    stopifnot(fit$jerr == 0, length(fit$lambda) == length(lambda_grid))
    pred <- predict(fit, z$v, s = lambda_grid)
    sse[j, ] <- sse[j, ] + colSums((pred - va$response_value)^2)
  }
}

grid <- do.call(rbind, lapply(seq_along(alpha_grid), function(j) {
  data.frame(alpha = alpha_grid[j], lambda = lambda_grid, OOF_RMSE = sqrt(sse[j, ] / nrow(dev)))
}))
wc(grid, "development_grid.csv")

choose <- function(d) {
  candidates <- d[d$OOF_RMSE <= min(d$OOF_RMSE) + 1e-12, ]
  candidates[order(-candidates$lambda, candidates$alpha), ][1, ]
}
selected <- rbind(cbind(model = "Ridge", choose(grid[grid$alpha == 0, ])),
                  cbind(model = "ElasticNet", choose(grid[grid$alpha > 0, ])))
wc(selected, "selected_on_development.csv")
wc(data.frame(model = "OLS", OOF_RMSE = sqrt(ols_sse / nrow(dev))), "OLS_development.csv")
if (any(selected$lambda %in% range(lambda_grid))) {
  warning("Selected lambda at grid boundary. Report this; do not retune using test results.")
}

test <- pairs[a$split == "test", ]
z <- features(dev, test, "results/06_internal/shared.rds", "results/06_internal/O2PLS_transform.rds")
pred <- data.frame(pair_id = test$pair_id, cell = test$cell_name, drug = test$drug_key,
                   observed = test$response_value)

ols <- ols_predict(z$x, dev$response_value, z$v)
pred$OLS <- ols$prediction
saveRDS(ols$fit, file.path(out, "OLS_fit.rds"))
for (i in seq_len(nrow(selected))) {
  fit <- glmnet::glmnet(z$x, dev$response_value, alpha = selected$alpha[i], lambda = lambda_grid,
                        standardize = FALSE, intercept = TRUE, thresh = 1e-10, maxit = 1000000)
  stopifnot(fit$jerr == 0, length(fit$lambda) == length(lambda_grid))
  pred[[selected$model[i]]] <- drop(predict(fit, z$v, s = selected$lambda[i]))
  saveRDS(fit, file.path(out, paste0(selected$model[i], "_fit.rds")))
}

reference <- read.csv("results/06_internal/predictions.csv", check.names = FALSE)
stopifnot(all(c("pair_id", "observed", "s15_O2PLS", "drug_mean", "drug_tissue_shrinkage10") %in%
                names(reference)))
ix <- match(pred$pair_id, reference$pair_id)
stopifnot(!anyNA(ix), !anyDuplicated(reference$pair_id), all(pred$observed == reference$observed[ix]))
for (name in c("s15_O2PLS", "drug_mean", "drug_tissue_shrinkage10")) pred[[name]] <- reference[[name]][ix]

metric <- function(y, p) {
  c(RMSE = sqrt(mean((y - p)^2)),
    MAE = mean(abs(y - p)),
    R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
}
models <- names(pred)[5:ncol(pred)]
stopifnot(all(is.finite(as.matrix(pred[, models]))))

summary <- do.call(rbind, lapply(models, function(m) {
  data.frame(model = m, t(metric(pred$observed, pred[[m]])))
}))
per_drug <- do.call(rbind, lapply(models, function(m) {
  do.call(rbind, lapply(split(pred, pred$drug), function(d) {
    data.frame(model = m, drug = d$drug[1], n = nrow(d), t(metric(d$observed, d[[m]])),
               Spearman = if (sd(d[[m]]) > 0 && sd(d$observed) > 0) {
                 cor(d$observed, d[[m]], method = "spearman")
               } else {
                 NA_real_
               })
  }))
}))

wc(pred, "test_predictions.csv")
wc(summary, "test_metrics.csv")
wc(per_drug, "per_drug_metrics.csv")
saveRDS(z$scaler, file.path(out, "component_scaler.rds"))
capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
writeLines("COMPLETE: descriptive post hoc comparison. No statistical superiority test performed.",
           file.path(out, "STATUS.txt"))
print(selected)
print(summary)
cat("\nOutputs:", normalizePath(out), "\n")
