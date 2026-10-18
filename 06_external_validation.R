local({
    bundle <- readRDS("data/data.rds")
    for (name in names(bundle)) {
        path <- file.path("data", name)
        if (!file.exists(path)) {
            dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
            writeBin(bundle[[name]], path)
        }
    }
    inputs <- c("raw_data/Gene_expression.csv", "raw_data/Methylation.csv",
                "raw_data/Drug_response_IC50.csv", "raw_data/PubchemFP.csv",
                "prepared_inputs/external_expression_fit_inputs.rds",
                "prepared_inputs/external_reduced_fit_inputs.rds")
    if (!all(file.exists(inputs))) {
        parts <- sort(list.files("data", pattern = "^inputs\\.rds\\.part[0-9]+$", full.names = TRUE))
        joined <- tempfile(fileext = ".rds")
        con <- file(joined, "wb")
        for (p in parts) writeBin(readBin(p, "raw", file.info(p)$size), con)
        close(con)
        bundle <- readRDS(joined)
        unlink(joined)
        for (name in names(bundle)) {
            if (!file.exists(name)) writeBin(bundle[[name]], name)
        }
    }
})

stopifnot(file.exists("O2MCA.Rproj"))
repo <- normalizePath(".")

result_dir <- function(name) file.path(repo, "results", name)

compare_external_refits <- function(repo, out) {
    for (kind in c("expression", "reduced")) {
        d <- read.csv(file.path(repo, "data", paste0("external_", kind, "_predictions.csv")), check.names = FALSE)
        ms <- if (kind == "expression") 
            c("MCA", "RAW")
        else c("O2PLS", "PCA_concat", "Expression")
        bundle <- readRDS(file.path(repo, "prepared_inputs", paste0("external_", kind, "_fit_inputs.rds")))
        means <- tapply(bundle$train_pairs$response_value, bundle$train_pairs$drug_key, mean)
        d$DRUG_MEAN <- as.numeric(means[d$drug_key])
        stopifnot(all(is.finite(d$DRUG_MEAN)))
        comparisons <- list()
        for (m in ms) {
            route <- if (m == "Expression") 
                "MCA"
            else m
            a <- read.csv(file.path(out, route, paste0(route, "_external_predictions_without_targets.csv")))
            ix <- match(d$external_pair_id, a$external_pair_id)
            stopifnot(!anyNA(ix), !anyDuplicated(a$external_pair_id), all(a$cell_name[ix] == d$cell_name), 
                all(a$drug_key[ix] == d$drug_key))
            comparisons[[m]] <- data.frame(model = m, max_abs_prediction_difference = max(abs(a$prediction[ix] - 
                d[[m]])), refit_RMSE = sqrt(mean((a$prediction[ix] - d$response_value)^2)), recorded_RMSE = sqrt(mean((d[[m]] - 
                d$response_value)^2)))
            d[[m]] <- a$prediction[ix]
        }
        write.csv(d, file.path(out, paste0(kind, "_refit_predictions_with_targets.csv")), row.names = FALSE)
        z <- do.call(rbind, comparisons)
        z$RMSE_difference <- z$refit_RMSE - z$recorded_RMSE
        write.csv(z, file.path(out, paste0(kind, "_comparison_with_manuscript.csv")), row.names = FALSE)
        print(z)
        per <- do.call(rbind, lapply(c(ms, "DRUG_MEAN"), function(m) do.call(rbind, lapply(unique(d$drug_key), 
            function(dr) {
                k <- d$drug_key == dr
                data.frame(model = m, drug = dr, n = sum(k), RMSE = sqrt(mean((d[k, m] - d$response_value[k])^2)))
            }))))
        write.csv(per, file.path(out, paste0(kind, "_refit_per_drug_RMSE.csv")), row.names = FALSE)
    }
}

run_external_expression <- function(repo, root, route, mode = "final") {
    args <- c(mode, route, file.path(repo, "prepared_inputs/external_expression_fit_inputs.rds"), 
        root)
    input <- normalizePath(args[3])
    root <- path.expand(args[4])
    if (file.exists(root)) 
        stop("Use a new output directory")

    require_model_packages()
    expected <- jsonlite::fromJSON(file.path(repo, "config/prepared_input_hashes.json"))
    stopifnot(sha(input) == expected[["external_expression_fit_inputs.rds"]])
    if (mode == "check") {
        cat("PASS: archived input hash and model dependencies\n")
        return(invisible(TRUE))
    }
    dir.create(root, recursive = TRUE)
    dir.create(file.path(root, "models"))
    b <- readRDS(input)
    grid <- read.csv(file.path(repo, "data/tuning_grid.csv"))
    if (mode == "final") 
        file.copy(file.path(repo, "data", paste0("external_expression_", route, "_selected.csv")), 
            file.path(root, paste0(route, "_selected.csv")))
    seed <- b$seed
    params <- function(g) {
        z <- xgb_params(g)
        z$seed <- seed
        z$nthread <- 2L
        z
    }
    fit_features <- function(tc, vc, tr, va, external = FALSE) {
        x <- b$expression[tc, , drop = FALSE]
        sx <- scale_fit(x)
        xs <- scale_apply(x, sx)
        set.seed(seed)
        pc <- prcomp(xs, center = FALSE, scale. = FALSE, rank. = 20)
        os <- pc$x[, 1:20, drop = FALSE]
        v <- if (external) 
            b$external_expression[vc, , drop = FALSE]
        else b$expression[vc, , drop = FALSE]
        vs <- predict(pc, newdata = scale_apply(v, sx))[, 1:20, drop = FALSE]
        drugs <- sort(unique(tr$drug_key))
        fp <- b$fp[drugs, , drop = FALSE]
        k <- apply(fp, 2, function(z) length(unique(z)) > 1)
        fp <- fp[, k, drop = FALSE]
        mca <- NULL
        if (route == "MCA") {
            fac <- factor_fp(fp)
            set.seed(seed)
            mca <- FactoMineR::MCA(fac, ncp = 20, graph = FALSE)
            ds <- mca$ind$coord[, 1:20, drop = FALSE]
        }
        else ds <- fp
        train <- cbind(os[tr$cell_name, , drop = FALSE], ds[tr$drug_key, , drop = FALSE])
        valid <- cbind(vs[va$cell_name, , drop = FALSE], ds[va$drug_key, , drop = FALSE])
        colnames(train) <- colnames(valid) <- paste0("F", seq_len(ncol(train)))
        sc <- scale_fit(train)
        train <- scale_apply(train, sc)
        valid <- scale_apply(valid, sc)
        list(train = train, valid = valid, bundle = list(expression_scaler = sx, pca = pc, mca = mca, 
            drug_scores = ds, component_scaler = sc, train_cells = tc, train_drugs = drugs))
    }
    if (mode == "cv") {
        for (fold in 1:5) {
            tr <- b$train_pairs[b$folds[b$train_pairs$cell_name] != fold, ]
            va <- b$train_pairs[b$folds[b$train_pairs$cell_name] == fold, ]
            tc <- sort(unique(tr$cell_name))
            vc <- sort(unique(va$cell_name))
            stopifnot(!any(tc %in% vc))
            cat("PREP", route, fold, as.character(Sys.time()), "\n")
            flush.console()
            z <- fit_features(tc, vc, tr, va)
            dt <- xgboost::xgb.DMatrix(z$train, label = tr$response_value)
            dv <- xgboost::xgb.DMatrix(z$valid)
            rm(z)
            gc()
            for (base in sort(unique(grid$base_id))) {
                path <- file.path(root, "models", paste0(route, "_fold", fold, "_base", base, ".rds"))
                if (file.exists(path)) {
                  stopifnot(sha(path) == readLines(paste0(path, ".sha256")))
                  next
                }
                g <- grid[grid$base_id == base, ]
                set.seed(seed)
                t0 <- Sys.time()
                model <- xgboost::xgb.train(params = params(g[1, ]), data = dt, nrounds = max(g$nrounds), 
                  verbose = 0)
                pred <- sapply(g$nrounds, function(n) as.numeric(predict(model, dv, iterationrange = c(1L, 
                  as.integer(n)))))
                colnames(pred) <- g$config_id
                if (fold == 1 && base == 1) {
                  set.seed(seed)
                  short <- xgboost::xgb.train(params = params(g[1, ]), data = dt, nrounds = min(g$nrounds), 
                    verbose = 0)
                  stopifnot(max(abs(predict(short, dv) - pred[, 1])) < 1e-10)
                  rm(short)
                }
                saveRDS(list(pair_id = va$pair_id, predictions = pred, config_id = g$config_id, 
                  fold = fold, base = base, route = route), path)
                writeLines(sha(path), paste0(path, ".sha256"))
                cat("DONE", route, fold, base, "seconds", as.numeric(difftime(Sys.time(), t0, units = "secs")), 
                  "\n")
                flush.console()
                rm(model, pred)
                gc()
            }
            rm(dt, dv)
            gc()
        }
        sums <- setNames(numeric(nrow(grid)), grid$config_id)
        n <- setNames(integer(nrow(grid)), grid$config_id)
        for (fold in 1:5) for (base in 1:16) {
            q <- readRDS(file.path(root, "models", paste0(route, "_fold", fold, "_base", base, ".rds")))
            truth <- b$train_pairs$response_value[match(q$pair_id, b$train_pairs$pair_id)]
            for (id in q$config_id) {
                sums[id] <- sums[id] + sum((q$predictions[, id] - truth)^2)
                n[id] <- n[id] + length(truth)
            }
        }
        stopifnot(all(n == nrow(b$train_pairs)))
        grid$pooled_RMSE <- sqrt(sums[grid$config_id]/n[grid$config_id])
        best <- grid[grid$pooled_RMSE <= min(grid$pooled_RMSE) + 1e-08, ]
        best <- best[order(best$max_depth, best$nrounds, best$config_id), ][1, ]
        write.csv(grid, file.path(root, paste0(route, "_cv_metrics.csv")), row.names = FALSE)
        write.csv(best, file.path(root, paste0(route, "_selected.csv")), row.names = FALSE)
        cat("CV_COMPLETE", route, "\n")
    }
    else {
        g <- read.csv(file.path(root, paste0(route, "_selected.csv")))
        tr <- b$train_pairs
        va <- b$external_pairs
        z <- fit_features(sort(unique(tr$cell_name)), sort(unique(va$cell_name)), tr, va, TRUE)
        dt <- xgboost::xgb.DMatrix(z$train, label = tr$response_value)
        dv <- xgboost::xgb.DMatrix(z$valid)
        set.seed(seed)
        model <- xgboost::xgb.train(params = params(g), data = dt, nrounds = g$nrounds, verbose = 0)
        pred <- as.numeric(predict(model, dv))
        stopifnot(all(is.finite(pred)))
        saveRDS(z$bundle, file.path(root, "models", paste0(route, "_final_transform.rds")))
        xgboost::xgb.save(model, file.path(root, "models", paste0(route, "_final.json")))
        q <- va
        q$prediction <- pred
        write.csv(q, file.path(root, paste0(route, "_external_predictions_without_targets.csv")), 
            row.names = FALSE)
        cat("PREDICTIONS_SAVED", route, "\n")
    }
    capture.output(sessionInfo(), file = file.path(root, paste0(route, "_", mode, "_session.txt")))
}

run_external_inference <- function(repo, input_dir, out) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    rd <- function(n) {
        if (n == "external_expression_predictions.csv") 
            return(read.csv(file.path(input_dir, "expression_refit_predictions_with_targets.csv")))
        if (n == "external_reduced_predictions.csv") 
            return(read.csv(file.path(input_dir, "reduced_refit_predictions_with_targets.csv")))
        read.csv(file.path(repo, "data", n))
    }
    wc <- function(x, n) write.csv(x, file.path(out, n), row.names = FALSE)
    for (kind in c("expression", "reduced")) {
        message("Recalculating 20,000 external bootstrap samples: ", kind)
        e <- rd(paste0("external_", kind, "_predictions.csv"))
        ms <- if (kind == "expression") 
            c("MCA", "RAW")
        else c("O2PLS", "PCA_concat", "Expression")
        groups <- unique(e$registry_component)
        drugs <- unique(e$drug_key)
        cnt <- matrix(0, length(groups), length(drugs))
        ss <- array(0, c(length(groups), length(drugs), length(ms)))
        for (i in seq_along(groups)) for (j in seq_along(drugs)) {
            k <- e$registry_component == groups[i] & e$drug_key == drugs[j]
            cnt[i, j] <- sum(k)
            for (m in seq_along(ms)) ss[i, j, m] <- sum((e[k, ms[m]] - e$response_value[k])^2)
        }
        set.seed(20260923)
        draws <- matrix(NA_real_, 20000, length(ms) - 1)
        for (i in 1:20000) {
            w <- tabulate(sample.int(length(groups), length(groups), replace = TRUE), nbins = length(groups))
            nn <- colSums(cnt * w)
            if (all(nn > 0)) {
                v <- sapply(seq_along(ms), function(m) mean(sqrt(colSums(ss[, , m] * w)/nn)))
                draws[i, ] <- v[1] - v[-1]
            }
        }
        point <- sapply(ms, function(m) mean(vapply(drugs, function(dr) {
            k <- e$drug_key == dr
            sqrt(mean((e[k, m] - e$response_value[k])^2))
        }, 0)))
        tails <- if (kind == "expression") 
            c(0.05/6, 1 - 0.05/6)
        else c(0.0125, 0.9875)
        result <- do.call(rbind, lapply(seq_len(ncol(draws)), function(j) {
            ci <- quantile(draws[, j], tails, na.rm = TRUE, names = FALSE)
            data.frame(difference = point[1] - point[j + 1], lower = ci[1], upper = ci[2])
        }))
        ref <- rd(if (kind == "expression") 
            "external_expression_contrast.csv"
        else "external_reduced_contrasts.csv")
        result$contrast <- ref$contrast
        wc(result, paste0("recalculated_external_", kind, "_contrasts.csv"))
    }
    e <- rd("external_expression_predictions.csv")
    k <- e$below_max_tested_dose
    sensitivity <- do.call(rbind, lapply(c("MCA", "RAW", "DRUG_MEAN"), function(m) cbind(model = m, 
        metrics(e$response_value[k], e[k, m]))))
    wc(sensitivity, "below_max_dose_sensitivity.csv")
}

run_external_reduced <- function(repo, root, route, mode = "final") {
    args <- c(mode, route, file.path(repo, "prepared_inputs/external_reduced_fit_inputs.rds"), root)
    input <- normalizePath(args[3])
    root <- path.expand(args[4])
    if (file.exists(root)) 
        stop("Use a new output directory")

    require_model_packages()
    expected <- jsonlite::fromJSON(file.path(repo, "config/prepared_input_hashes.json"))
    stopifnot(sha(input) == expected[["external_reduced_fit_inputs.rds"]])
    if (mode == "check") {
        cat("PASS: archived input hash and model dependencies\n")
        return(invisible(TRUE))
    }
    dir.create(root, recursive = TRUE)
    dir.create(file.path(root, "models"))
    b <- readRDS(input)
    grid <- read.csv(file.path(repo, "data/tuning_grid.csv"))
    if (mode == "final") 
        file.copy(file.path(repo, "data", paste0("external_reduced_", route, "_selected.csv")), 
            file.path(root, paste0(route, "_selected.csv")))
    seed <- 20260923L
    params <- function(g) {
        z <- xgb_params(g)
        z$seed <- seed
        z$nthread <- 2L
        z
    }
    fill <- function(x, med) {
        for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- med[j]
        stopifnot(all(is.finite(x)))
        x
    }
    features <- function(tc, vc, tr, va, external = FALSE, label) {
        x <- b$expression[tc, , drop = FALSE]
        y <- b$methylation[tc, , drop = FALSE]
        med <- apply(y, 2, median, na.rm = TRUE)
        stopifnot(all(is.finite(med)))
        y <- fill(y, med)
        sx <- scale_fit(x)
        sy <- scale_fit(y)
        xx <- scale_apply(x, sx)
        yy <- scale_apply(y, sy)
        om <- fit_omics(route, xx, yy, file.path(root, "models", paste0(route, "_", label, "_convergence.txt")))
        os <- project_omics(om, xx, yy)
        vx <- if (external) 
            b$external_expression[vc, , drop = FALSE]
        else b$expression[vc, , drop = FALSE]
        vy <- if (external) 
            b$external_methylation[vc, , drop = FALSE]
        else b$methylation[vc, , drop = FALSE]
        vs <- project_omics(om, scale_apply(vx, sx), scale_apply(fill(vy, med), sy))
        drugs <- sort(unique(tr$drug_key))
        fp <- b$fp[drugs, , drop = FALSE]
        fp <- fp[, apply(fp, 2, function(z) length(unique(z)) > 1), drop = FALSE]
        set.seed(seed)
        mca <- FactoMineR::MCA(factor_fp(fp), ncp = 20, graph = FALSE)
        ds <- mca$ind$coord[, 1:20, drop = FALSE]
        train <- cbind(os[tr$cell_name, , drop = FALSE], ds[tr$drug_key, , drop = FALSE])
        valid <- cbind(vs[va$cell_name, , drop = FALSE], ds[va$drug_key, , drop = FALSE])
        colnames(train) <- colnames(valid) <- paste0("F", seq_len(ncol(train)))
        sc <- scale_fit(train)
        list(train = scale_apply(train, sc), valid = scale_apply(valid, sc), bundle = list(expression_scaler = sx, 
            methylation_scaler = sy, methylation_medians = med, omics = om, mca = mca, drug_scores = ds, 
            component_scaler = sc, train_cells = tc))
    }
    if (mode == "cv") {
        for (fold in 1:5) {
            tr <- b$train_pairs[b$folds[b$train_pairs$cell_name] != fold, ]
            va <- b$train_pairs[b$folds[b$train_pairs$cell_name] == fold, ]
            cat("PREP", route, fold, "\n")
            flush.console()
            z <- features(sort(unique(tr$cell_name)), sort(unique(va$cell_name)), tr, va, label = paste0("fold", 
                fold))
            dt <- xgboost::xgb.DMatrix(z$train, label = tr$response_value)
            dv <- xgboost::xgb.DMatrix(z$valid)
            rm(z)
            gc()
            for (base in 1:16) {
                path <- file.path(root, "models", paste0(route, "_fold", fold, "_base", base, ".rds"))
                if (file.exists(path)) {
                  stopifnot(sha(path) == readLines(paste0(path, ".sha256")))
                  next
                }
                g <- grid[grid$base_id == base, ]
                set.seed(seed)
                m <- xgboost::xgb.train(params = params(g[1, ]), data = dt, nrounds = max(g$nrounds), 
                  verbose = 0)
                pred <- sapply(g$nrounds, function(n) as.numeric(predict(m, dv, iterationrange = c(1L, 
                  as.integer(n)))))
                colnames(pred) <- g$config_id
                saveRDS(list(pair_id = va$pair_id, predictions = pred, config_id = g$config_id), 
                  path)
                writeLines(sha(path), paste0(path, ".sha256"))
                cat("DONE", route, fold, base, as.character(Sys.time()), "\n")
                flush.console()
                rm(m, pred)
                gc()
            }
        }
        sums <- setNames(numeric(nrow(grid)), grid$config_id)
        counts <- sums
        for (fold in 1:5) for (base in 1:16) {
            q <- readRDS(file.path(root, "models", paste0(route, "_fold", fold, "_base", base, ".rds")))
            truth <- b$train_pairs$response_value[match(q$pair_id, b$train_pairs$pair_id)]
            for (id in q$config_id) {
                sums[id] <- sums[id] + sum((q$predictions[, id] - truth)^2)
                counts[id] <- counts[id] + length(truth)
            }
        }
        stopifnot(all(counts == nrow(b$train_pairs)))
        grid$pooled_RMSE <- sqrt(sums[grid$config_id]/counts[grid$config_id])
        best <- grid[grid$pooled_RMSE <= min(grid$pooled_RMSE) + 1e-08, ]
        best <- best[order(best$max_depth, best$nrounds, best$config_id), ][1, ]
        write.csv(grid, file.path(root, paste0(route, "_cv_metrics.csv")), row.names = FALSE)
        write.csv(best, file.path(root, paste0(route, "_selected.csv")), row.names = FALSE)
    }
    else {
        g <- read.csv(file.path(root, paste0(route, "_selected.csv")))
        tr <- b$train_pairs
        va <- b$external_pairs
        z <- features(sort(unique(tr$cell_name)), sort(unique(va$cell_name)), tr, va, TRUE, "final")
        set.seed(seed)
        m <- xgboost::xgb.train(params = params(g), data = xgboost::xgb.DMatrix(z$train, label = tr$response_value), 
            nrounds = g$nrounds, verbose = 0)
        pred <- as.numeric(predict(m, xgboost::xgb.DMatrix(z$valid)))
        stopifnot(all(is.finite(pred)))
        saveRDS(z$bundle, file.path(root, "models", paste0(route, "_final_transform.rds")))
        xgboost::xgb.save(m, file.path(root, "models", paste0(route, "_final.json")))
        va$prediction <- pred
        write.csv(va, file.path(root, paste0(route, "_external_predictions_without_targets.csv")), 
            row.names = FALSE)
    }
    capture.output(sessionInfo(), file = file.path(root, paste0(route, "_", mode, "_session.txt")))
}

metrics <- function(y, p, min_cor_n = 2L) {
    stopifnot(length(y) == length(p), length(y) > 0L, all(is.finite(y)), all(is.finite(p)))
    n <- length(y)
    sse <- sum((y - p)^2)
    sst <- sum((y - mean(y))^2)
    cor_ok <- n >= min_cor_n && sd(y) > 0 && sd(p) > 0
    data.frame(n = n, RMSE = sqrt(sse/n), MAE = mean(abs(y - p)), R2 = if (n >= 2 && sst > 0) 
        1 - sse/sst
    else NA_real_, Pearson = if (cor_ok) 
        cor(y, p)
    else NA_real_, Spearman = if (cor_ok) 
        cor(y, p, method = "spearman")
    else NA_real_, R2_reason = if (n < 2) 
        "n_lt_2"
    else if (sst == 0) 
        "constant_observed"
    else "defined", correlation_reason = if (n < min_cor_n) 
        "insufficient_n"
    else if (sd(y) == 0 || sd(p) == 0) 
        "constant_vector"
    else "defined")
}

require_model_packages <- function() {
    for (p in c("OmicsPLS", "FactoMineR", "xgboost", "digest", "r.jive", "jsonlite")) stopifnot(requireNamespace(p, 
        quietly = TRUE))
}

sha <- function(p) digest::digest(file = p, algo = "sha256")

scale_fit <- function(x, drop_constant = TRUE) {
    s <- apply(x, 2, sd)
    k <- is.finite(s) & s > 0
    if (!drop_constant) 
        stopifnot(all(k))
    list(names = colnames(x)[k], center = colMeans(x)[k], scale = s[k])
}

scale_apply <- function(x, s) {
    stopifnot(all(s$names %in% colnames(x)))
    z <- sweep(sweep(x[, s$names, drop = FALSE], 2, s$center, "-"), 2, s$scale, "/")
    stopifnot(all(is.finite(z)))
    z
}

factor_fp <- function(x, lev = NULL) {
    z <- as.data.frame(lapply(seq_len(ncol(x)), function(j) {
        v <- paste0(colnames(x)[j], "=", x[, j])
        l <- if (is.null(lev)) 
            sort(unique(v))
        else lev[[j]]
        stopifnot(all(v %in% l))
        factor(v, levels = l)
    }), check.names = FALSE)
    names(z) <- colnames(x)
    rownames(z) <- rownames(x)
    z
}

joint_scores <- function(x, w, highd) {
    if (!highd) 
        return(x %*% w)
    z <- matrix(0, nrow(x), ncol(w))
    for (j in seq_len(ncol(w))) {
        z[, j] <- x %*% w[, j]
        x <- x - tcrossprod(z[, j], w[, j])
    }
    z
}

project_omics <- function(m, x, y) {
    stopifnot(identical(rownames(x), rownames(y)))
    z <- switch(m$method, O2PLS = {
        f <- m$fit
        tx <- x %*% f$W_Yosc
        uy <- y %*% f$C_Xosc
        cbind(joint_scores(x - tx %*% t(f$P_Yosc.), f$W., isTRUE(f$flags$highd)), tx, uy)
    }, PLS = joint_scores(x, m$fit$W., isTRUE(m$fit$flags$highd)), PCA_concat = predict(m$fit, newdata = cbind(x, 
        y)), PCA_expression = predict(m$fit, newdata = x), PCA_methylation = predict(m$fit, newdata = y), 
        PCA_separate = cbind(predict(m$gene, newdata = x), predict(m$meth, newdata = y)), JIVE = t(qr.coef(m$loading_qr, 
            t(cbind(x, y)))))
    rownames(z) <- rownames(x)
    colnames(z) <- paste0("Omics_", seq_len(ncol(z)))
    stopifnot(ncol(z) == 20L, all(is.finite(z)))
    z
}

fit_omics <- function(method, x, y, logfile) {
    set.seed(20260921L)
    m <- list(method = method)
    messages <- character()
    if (method == "PCA_concat") 
        m$fit <- prcomp(cbind(x, y), center = FALSE, scale. = FALSE, rank. = 20)
    if (method == "PCA_expression") 
        m$fit <- prcomp(x, center = FALSE, scale. = FALSE, rank. = 20)
    if (method == "PCA_methylation") 
        m$fit <- prcomp(y, center = FALSE, scale. = FALSE, rank. = 20)
    if (method == "PCA_separate") {
        m$gene <- prcomp(x, center = FALSE, scale. = FALSE, rank. = 10)
        m$meth <- prcomp(y, center = FALSE, scale. = FALSE, rank. = 10)
    }
    if (method %in% c("O2PLS", "PLS")) {
        m$fit <- withCallingHandlers(OmicsPLS::o2m(x, y, n = if (method == "O2PLS") 
            10
        else 20, nx = if (method == "O2PLS") 
            5
        else 0, ny = if (method == "O2PLS") 
            5
        else 0, stripped = TRUE, tol = 1e-10, max_iterations = 10000), message = function(z) {
            messages <<- c(messages, conditionMessage(z))
            invokeRestart("muffleMessage")
        })
        writeLines(messages, logfile)
        if (isTRUE(m$fit$flags$highd)) {
            z <- messages[grepl("stopped after", messages, fixed = TRUE)]
            ni <- as.integer(sub(".*stopped after ([0-9]+) iterations.*", "\\1", z))
            stopifnot(length(ni) > 0, all(ni < 10000))
        }
    }
    if (method == "JIVE") {
        messages <- capture.output(j <- r.jive::jive(list(t(x), t(y)), rankJ = 10, rankA = c(5, 
            5), method = "given", center = FALSE, scale = FALSE, orthIndiv = TRUE, est = TRUE, maxiter = 1000, 
            showProgress = TRUE))
        writeLines(messages, logfile)
        stopifnot(any(grepl("JIVE algorithm converged after", messages, fixed = TRUE)))
        jl <- svd(do.call(rbind, j$joint), nu = 10, nv = 0)$u
        i1 <- svd(j$individual[[1]], nu = 5, nv = 0)$u
        i2 <- svd(j$individual[[2]], nu = 5, nv = 0)$u
        m$loadings <- cbind(jl, rbind(i1, matrix(0, ncol(y), 5)), rbind(matrix(0, ncol(x), 5), i2))
        m$loading_qr <- qr(m$loadings)
        stopifnot(m$loading_qr$rank == 20)
    }
    st <- project_omics(m, x, y)
    target <- switch(method, O2PLS = cbind(m$fit$Tt, m$fit$T_Yosc, m$fit$U_Xosc), PLS = m$fit$Tt, 
        PCA_concat = m$fit$x, PCA_expression = m$fit$x, PCA_methylation = m$fit$x, PCA_separate = cbind(m$gene$x, 
            m$meth$x), JIVE = st)
    stopifnot(max(abs(st - target)) <= 1e-07 * max(1, max(abs(target))), qr(st)$rank == 20)
    m
}

xgb_params <- function(g) {
    c(list(objective = "reg:squarederror", tree_method = "hist", subsample = 1, colsample_bytree = 1, 
        alpha = 0, gamma = 0, nthread = 1, seed = 20260921L), list(max_depth = as.integer(g$max_depth), 
        eta = as.numeric(g$eta), min_child_weight = as.numeric(g$min_child_weight), lambda = as.numeric(g$lambda)))
}

for (route in c("MCA", "RAW")) {
  run_external_expression(repo, result_dir(paste0("12_external_cv/", route)), route, mode = "cv")
}
for (route in c("O2PLS", "PCA_concat")) {
  run_external_reduced(repo, result_dir(paste0("12_external_cv/", route)), route, mode = "cv")
}

for (route in c("MCA", "RAW")) {
  run_external_expression(repo, result_dir(paste0("07_external/", route)), route)
}
for (route in c("O2PLS", "PCA_concat")) {
  run_external_reduced(repo, result_dir(paste0("07_external/", route)), route)
}
compare_external_refits(repo, result_dir("07_external"))

run_external_inference(repo, result_dir("07_external"), result_dir("13_external_inference"))
