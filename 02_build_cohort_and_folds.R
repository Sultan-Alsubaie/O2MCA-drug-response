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

cluster_rmse <- function(y, pred, cell, maps) {
    stopifnot(nrow(pred) == length(y), length(cell) == length(y), all(is.finite(y)), all(is.finite(pred)), 
        !anyNA(cell))
    ids <- sort(unique(cell))
    stopifnot(nrow(maps) == length(ids), all(maps >= 1 & maps <= length(ids)))
    e <- sweep(pred, 1, y, "-")^2
    ss <- rowsum(e, group = factor(cell, levels = ids), reorder = TRUE)
    stopifnot(identical(rownames(ss), ids))
    counts <- as.numeric(table(factor(cell, levels = ids)))
    vapply(seq_len(ncol(maps)), function(i) {
        w <- tabulate(maps[, i], nbins = length(ids))
        sqrt(colSums(ss * w)/sum(counts * w))
    }, numeric(ncol(pred)))
}

infer_contrasts <- function(y, pred, cell, contrasts, B = 20000L, seed = 20260923L) {
    ids <- sort(unique(cell))
    set.seed(seed)
    maps <- matrix(sample.int(length(ids), length(ids) * B, replace = TRUE), nrow = length(ids))
    draws <- t(cluster_rmse(y, pred, cell, maps))
    colnames(draws) <- colnames(pred)
    point <- sqrt(colMeans(sweep(pred, 1, y, "-")^2))
    family <- length(contrasts)
    stopifnot(family == 6L)
    tab <- do.call(rbind, lapply(contrasts, function(c) {
        d <- draws[, c$left] - draws[, c$right]
        q <- quantile(d, c(0.025, 0.975, 0.05/(2 * family), 1 - 0.05/(2 * family)), type = 7, names = FALSE)
        data.frame(contrast = c$id, left = c$left, right = c$right, delta_RMSE = point[c$left] - 
            point[c$right], lower95 = q[1], upper95 = q[2], lower_adjusted = q[3], upper_adjusted = q[4], 
            adjusted_level = 1 - 0.05/family, B = B, conclusion = if (q[4] < 0) 
                "lower_RMSE_left"
            else if (q[3] > 0) 
                "lower_RMSE_right"
            else "inconclusive")
    }))
    list(contrasts = tab, draws = draws, maps = maps, cell_order = ids, seed = seed)
}

verify_inference <- function() {
    stopifnot(metrics(c(1, 2, 3), c(1, 2, 3))$RMSE == 0, metrics(c(1, 2, 3), c(1, 2, 3))$R2 == 1, 
        is.na(metrics(1, 2)$R2), is.na(metrics(rep(1, 5), 1:5)$Pearson), is.na(metrics(1:4, 1:4, 
            min_cor_n = 5)$Pearson))
    y <- c(1, 2, 3, 4, 5, 6)
    p <- cbind(A = y + c(1, 0, -1, 2, 0, 1), B = y + c(0, 1, 2, 0, 1, 0))
    cell <- c("a", "a", "b", "c", "c", "c")
    maps <- cbind(c(1, 1, 3), c(2, 3, 2), c(1, 2, 3))
    fast <- cluster_rmse(y, p, cell, maps)
    slow <- vapply(1:ncol(maps), function(b) {
        ix <- unlist(lapply(maps[, b], function(j) which(cell == sort(unique(cell))[j])))
        sqrt(colMeans((p[ix, , drop = FALSE] - y[ix])^2))
    }, numeric(2))
    stopifnot(max(abs(fast - slow)) < 1e-12, max(abs(cluster_rmse(y, p[, 2:1], cell, maps) - fast[2:1, 
        ])) < 1e-12)
    p2 <- cbind(A = p[, 1], B = p[, 1])
    cs <- lapply(1:6, function(i) list(id = paste0("C", i), left = "A", right = "B"))
    z <- infer_contrasts(y, p2, cell, cs, B = 200)
    stopifnot(all(z$contrasts$delta_RMSE == 0), all(z$contrasts$lower_adjusted == 0), all(z$contrasts$upper_adjusted == 
        0), abs(z$contrasts$adjusted_level[1] - (1 - 0.05/6)) < 1e-14)
    stopifnot(inherits(try(metrics(c(1, 2), c(NA, 2)), silent = TRUE), "try-error"))
    TRUE
}

run_internal <- function(repo, output, check_only = TRUE) {
    args <- c(file.path(repo, "prepared_inputs", c("protocol_config.rds", "omics_verified_headers.rds", 
        "PubchemFP.csv", "development_cell_metadata.csv")), output)
    if (check_only) 
        args <- c(args, "--check-only")


    require_model_packages()
    protocol_path <- normalizePath(args[1])
    omics_path <- normalizePath(args[2])
    fp_path <- normalizePath(args[3])
    metadata_path <- normalizePath(args[4])
    spec_path <- file.path(repo, "config/final_internal_choices.json")
    spec <- jsonlite::fromJSON(spec_path, simplifyVector = FALSE)
    expected <- jsonlite::fromJSON(file.path(repo, "config/internal_input_hashes.json"))
    inputs <- c(protocol = protocol_path, omics = omics_path, fingerprints = fp_path, metadata = metadata_path, 
        spec = spec_path)
    stopifnot(sha(protocol_path) == expected[["protocol_config.rds"]], sha(omics_path) == expected[["omics_verified_headers.rds"]], 
        sha(fp_path) == expected[["PubchemFP.csv"]], sha(metadata_path) == expected[["development_cell_metadata.csv"]])
    for (p in c("OmicsPLS", "FactoMineR", "xgboost", "digest", "r.jive", "jsonlite")) {
        if (as.character(packageVersion(p)) != spec$package_versions[[p]]) 
            stop("Recorded package version required: ", p, " ", spec$package_versions[[p]])
    }
    protocol <- readRDS(protocol_path)
    a <- protocol$cell_pair_split
    ap <- protocol$analysis_pairs
    stopifnot(identical(a$pair_id, ap$pair_id), !anyDuplicated(ap$pair_id))
    it <- which(a$split == "development")
    iv <- which(a$split == "test")
    tr <- ap[it, ]
    va <- ap[iv, setdiff(names(ap), "response_value"), drop = FALSE]
    tc <- sort(unique(tr$cell_name))
    vc <- sort(unique(va$cell_name))
    td <- sort(unique(tr$drug_key))
    vd <- sort(unique(va$drug_key))
    stopifnot(!any(tc %in% vc), !any(tr$cosmic_id %in% va$cosmic_id), all(vd %in% td), nrow(tr) == 
        99838L, nrow(va) == 25066L, length(tc) == 587L, length(vc) == 147L)
    if (length(args) == 6L) {
        cat("PASS: recorded data hashes, core package versions and original split. No fitting performed.\n")
        return(invisible(TRUE))
    }
    out <- path.expand(args[5])
    if (file.exists(out)) 
        stop("Choose a NEW output directory; existing runs are not overwritten.")
    dir.create(out, recursive = TRUE)
    out <- normalizePath(out)
    mode <- "final"
    before <- vapply(inputs, sha, "")
    jsonlite::write_json(list(purpose = "Reproduction of already reported results; not new confirmatory evidence", 
        inputs = as.list(before)), file.path(out, "REPRODUCTION_MANIFEST.json"), pretty = TRUE, 
        auto_unbox = TRUE)
    wc <- function(d, f) write.csv(d, file.path(out, f), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    checks <- data.frame(check = character(), passed = logical(), value = double())
    check <- function(n, b, value = NA_real_) {
        checks <<- rbind(checks, data.frame(check = n, passed = isTRUE(b), value = value))
        wc(checks, "validation_checks.csv")
        if (!isTRUE(b)) 
            stop(n)
    }
    near <- function(n, x, y, tol = 1e-07) {
        d <- if (identical(dim(x), dim(y))) 
            max(abs(x - y))
        else Inf
        check(n, all(is.finite(x)) && all(is.finite(y)) && d <= tol * max(1, max(abs(y))), d)
    }
    cache_save <- function(x, f) {
        p <- file.path(out, f)
        saveRDS(x, p, compress = FALSE)
        writeLines(sha(p), paste0(p, ".sha256"))
    }
    cache_load <- function(f) {
        p <- file.path(out, f)
        if (!file.exists(p)) 
            return(NULL)
        stopifnot(file.exists(paste0(p, ".sha256")), sha(p) == readLines(paste0(p, ".sha256"), n = 1))
        readRDS(p)
    }
    finish <- function(status) {
        writeLines(status, file.path(out, "STATUS.txt"))
        capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
        wc(data.frame(input = names(inputs), path = unname(inputs), sha256_before = unname(before), 
            sha256_after = unname(vapply(inputs, sha, ""))), "input_integrity.csv")
        fs <- list.files(out, full.names = TRUE)
        fs <- fs[basename(fs) != "output_manifest.csv"]
        wc(data.frame(file = basename(fs), sha256 = vapply(fs, sha, "")), "output_manifest.csv")
    }
    tryCatch({
        writeLines("RUNNING", file.path(out, "STATUS.txt"))
        message(mode, ": ", out)
        check("inference_known_answer_tests", verify_inference())
        check("cold_cell_and_cosmic_disjoint", !any(tc %in% vc) && !any(tr$cosmic_id %in% va$cosmic_id))
        check("known_drugs", all(vd %in% td))
        wc(rbind(data.frame(pair_id = tr$pair_id, cell = tr$cell_name, drug = tr$drug_key, role = "train"), 
            data.frame(pair_id = va$pair_id, cell = va$cell_name, drug = va$drug_key, role = mode)), 
            "pair_ids.csv")
        meta <- unique(ap[, c("cell_name", "cosmic_id", "primary_site")])
        check("one_metadata_record_per_cell", !anyDuplicated(meta$cell_name))
        tissue <- setNames(normal_tissue(meta$primary_site), meta$cell_name)
        m21 <- read.csv(metadata_path, stringsAsFactors = FALSE)
        check("development_tissue_matches_verified_S5", identical(unname(tissue[m21$cell]), normal_tissue(m21$primary_site)))
        lev <- sort(setdiff(unique(tissue[tc]), "__UNKNOWN__"))
        et <- encode_tissue(tc, tissue, lev)
        ev <- encode_tissue(vc, tissue, lev)
        check("unknown_tissue_all_zero", all(rowSums(ev)[!tissue[vc] %in% lev] == 0))
        check("synthetic_unseen_tissue", all(encode_tissue("sentinel", c(sentinel = "unseen"), lev) == 
            0))
        cache_save(list(levels = lev, train_cells = tc, unknown = "all_zero"), "tissue_encoder.rds")
        wc(data.frame(cell = vc, tissue = unname(tissue[vc]), known = rowSums(ev) == 1), "evaluation_tissue_coverage.csv")
        omics <- readRDS(omics_path)
        fp <- as.matrix(read.csv(fp_path, row.names = 1, check.names = FALSE))
        rownames(fp) <- trimws(rownames(fp))
        check("unique_fingerprint_ids", !anyDuplicated(rownames(fp)))
        shared <- cache_load("shared.rds")
        if (is.null(shared)) {
            message("fit training-only omics scalers and MCA")
            shared <- fit_shared(omics, fp, tc, td)
            cache_save(shared, "shared.rds")
        }
        check("shared_training_entities", identical(tc, shared$train_cells) && identical(td, shared$train_drugs))
        for (bl in c("expression", "methylation")) {
            s <- scale_fit(omics[[bl]][tc, , drop = FALSE])
            old <- shared[[paste0(bl, "_scaler")]]
            check(paste0(bl, "_feature_names"), identical(s$names, old$names))
            near(paste0(bl, "_center"), s$center, old$center)
            near(paste0(bl, "_scale"), s$scale, old$scale)
        }
        x <- scale_apply(omics$expression[tc, , drop = FALSE], shared$expression_scaler)
        y <- scale_apply(omics$methylation[tc, , drop = FALSE], shared$methylation_scaler)
        xv <- scale_apply(omics$expression[vc, , drop = FALSE], shared$expression_scaler)
        yv <- scale_apply(omics$methylation[vc, , drop = FALSE], shared$methylation_scaler)
        drt <- project_drug(shared, fp, td)
        drv <- project_drug(shared, fp, vd)
        near("MCA_training_projection", drt, shared$mca$ind$coord[td, , drop = FALSE])
        learned <- Filter(function(m) !is.null(m$selected_configuration), spec$models)
        routes <- lapply(learned, function(m) model_route(m$id))
        methods <- unique(vapply(routes, function(r) r$method, ""))
        scores <- list()
        for (method in methods) {
            message("prepare ", method)
            m <- cache_load(paste0(method, "_transform.rds"))
            if (is.null(m)) {
                m <- fit_omics(method, x, y, file.path(out, paste0(method, "_fit_messages.txt")))
                cache_save(m, paste0(method, "_transform.rds"))
            }
            st <- project_omics(m, x, y)
            sv <- project_omics(m, xv, yv)
            near(paste0(method, "_single_row"), project_omics(m, xv[1, , drop = FALSE], yv[1, , 
                drop = FALSE]), sv[1, , drop = FALSE])
            perm <- rev(seq_len(nrow(xv)))
            near(paste0(method, "_row_order"), project_omics(m, xv[perm, , drop = FALSE], yv[perm, 
                , drop = FALSE])[order(perm), , drop = FALSE], sv)
            check(paste0(method, "_training_rank"), qr(st)$rank == 20L)
            scores[[method]] <- list(train = st, evaluation = sv)
            rm(m)
            gc(FALSE)
        }
        pred <- matrix(NA_real_, nrow(va), length(spec$models), dimnames = list(as.character(va$pair_id), 
            vapply(spec$models, function(m) m$id, "")))
        for (m in learned) {
            r <- model_route(m$id)
            g <- m$selected_configuration
            os <- scores[[r$method]]
            pt <- feature_matrix(tr, drt, os$train, r$indices)
            pv <- feature_matrix(va, drv, os$evaluation, r$indices)
            sc <- scale_fit(pt, FALSE)
            cache_save(sc, paste0(m$id, "_component_scaler.rds"))
            pt <- scale_apply(pt, sc)
            pv <- scale_apply(pv, sc)
            if (r$tissue) {
                pt <- cbind(pt, et[tr$cell_name, , drop = FALSE])
                pv <- cbind(pv, ev[va$cell_name, , drop = FALSE])
            }
            check(paste0(m$id, "_features"), identical(colnames(pt), colnames(pv)) && ncol(pt) == 
                20 + length(r$indices) + (if (r$tissue) 
                  length(lev)
                else 0))
            res <- cache_load(paste0(m$id, "_predictions.rds"))
            bp <- file.path(out, paste0(m$id, ".json"))
            if (is.null(res)) {
                message("fit ", m$id, " ", g$config_id, " rounds ", g$nrounds)
                set.seed(20260921L)
                start <- proc.time()[3]
                b <- xgboost::xgb.train(params = xgb_params(g), data = xgboost::xgb.DMatrix(pt, 
                  label = tr$response_value, nthread = 1), nrounds = as.integer(g$nrounds), verbose = 0)
                v <- as.numeric(predict(b, xgboost::xgb.DMatrix(pv, nthread = 1)))
                xgboost::xgb.save(b, bp)
                res <- list(pair_id = va$pair_id, prediction = v, model_sha256 = sha(bp), seconds = unname(proc.time()[3] - 
                  start), configuration = g)
                cache_save(res, paste0(m$id, "_predictions.rds"))
                rm(b)
            }
            check(paste0(m$id, "_cached_ids_and_hash"), identical(res$pair_id, va$pair_id) && res$model_sha256 == 
                sha(bp))
            b <- xgboost::xgb.load(bp)
            again <- as.numeric(predict(b, xgboost::xgb.DMatrix(pv, nthread = 1)))
            near(paste0(m$id, "_reload_predictions"), again, res$prediction, 1e-10)
            pred[, m$id] <- again
            rm(pt, pv, b, os)
            gc(FALSE)
        }
        means <- tapply(tr$response_value, tr$drug_key, mean)
        bp <- as.numeric(means[va$drug_key])
        check("drug_mean_known", all(is.finite(bp)))
        pred[, "drug_mean"] <- bp
        kt <- paste(tr$drug_key, tissue[tr$cell_name], sep = "||")
        kv <- paste(va$drug_key, tissue[va$cell_name], sep = "||")
        ss <- tapply(tr$response_value, kt, sum)
        nn <- table(kt)
        s <- as.numeric(ss[kv])
        n <- as.numeric(nn[kv])
        missing <- is.na(s) | is.na(n) | tissue[va$cell_name] == "__UNKNOWN__"
        s[missing] <- 0
        n[missing] <- 0
        pred[, "drug_tissue_shrinkage10"] <- (s + 10 * bp)/(n + 10)
        cache_save(list(drug_means = means, tissue_sum = ss, tissue_n = nn, shrinkage = 10), "baseline_models.rds")
        check("all_predictions_finite", all(is.finite(pred)))
        check("15_paths", ncol(pred) == 15L)
        wc(data.frame(pair_id = va$pair_id, cell = va$cell_name, drug = va$drug_key, pred, check.names = FALSE), 
            "predictions_without_targets.csv")
        if (mode == "final" && !file.exists(file.path(out, "REPRODUCTION_TARGETS_READ_AT.txt"))) 
            writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), file.path(out, "REPRODUCTION_TARGETS_READ_AT.txt"))
        observed <- ap$response_value[iv]
        check("observed_targets_only", all(is.finite(observed)))
        wc(data.frame(pair_id = va$pair_id, cell = va$cell_name, drug = va$drug_key, tissue = unname(tissue[va$cell_name]), 
            observed, pred, check.names = FALSE), "predictions.csv")
        perdrug <- do.call(rbind, lapply(colnames(pred), function(id) do.call(rbind, lapply(sort(unique(va$drug_key)), 
            function(d) {
                ix <- va$drug_key == d
                cbind(model = id, drug = d, metrics(observed[ix], pred[ix, id], min_cor_n = 5L))
            }))))
        wc(perdrug, "per_drug_metrics.csv")
        pooled <- do.call(rbind, lapply(colnames(pred), function(id) {
            d <- metrics(observed, pred[, id])
            d$macro_drug_RMSE <- mean(perdrug$RMSE[perdrug$model == id])
            cbind(model = id, d)
        }))
        wc(pooled, "pooled_metrics.csv")
        calibration <- do.call(rbind, lapply(colnames(pred), function(id) {
            cf <- coef(lm(observed ~ pred[, id]))
            data.frame(model = id, intercept = unname(cf[1]), slope = unname(cf[2]))
        }))
        wc(calibration, "calibration.csv")
        by_tissue <- do.call(rbind, lapply(colnames(pred), function(id) do.call(rbind, lapply(sort(unique(tissue[va$cell_name])), 
            function(ti) {
                ix <- tissue[va$cell_name] == ti
                cbind(model = id, tissue = ti, metrics(observed[ix], pred[ix, id], 5L))
            }))))
        wc(by_tissue, "per_tissue_metrics.csv")
        message("20,000 paired cell bootstrap draws")
        inf <- infer_contrasts(observed, pred, va$cell_name, spec$primary_contrasts)
        wc(inf$contrasts, "primary_contrasts.csv")
        cache_save(inf, "bootstrap.rds")
        win <- do.call(rbind, lapply(spec$primary_contrasts, function(c) {
            l <- perdrug[perdrug$model == c$left, ]
            r <- perdrug[perdrug$model == c$right, ]
            delta <- l$RMSE - r$RMSE[match(l$drug, r$drug)]
            data.frame(contrast = c$id, drug = l$drug, delta_RMSE = delta, outcome = ifelse(abs(delta) <= 
                1e-08, "tie", ifelse(delta < 0, "left", "right")))
        }))
        wc(win, "per_drug_contrast_wins.csv")
        check("inputs_unchanged", identical(before, vapply(inputs, sha, "")))
        check("all_drugs_in_macro", all(table(perdrug$model) == length(vd)))
        finish(if (mode == "final") 
            "COMPLETE_REPRODUCTION_OF_REPORTED_HELDOUT_EVALUATION"
        else "DEVELOPMENT_PIPELINE_VALIDATED_TEST_UNTOUCHED")
        message("COMPLETE ", out)
    }, error = function(e) {
        finish(paste("FAILED_RESUMABLE", conditionMessage(e)))
        stop(e)
    })
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

fit_shared <- function(omics, fp, cells, drugs) {
    sx <- scale_fit(omics$expression[cells, , drop = FALSE])
    sy <- scale_fit(omics$methylation[cells, , drop = FALSE])
    f <- fp[drugs, , drop = FALSE]
    k <- apply(f, 2, function(v) length(unique(v)) > 1)
    f <- factor_fp(f[, k, drop = FALSE])
    set.seed(20260921L)
    list(expression_scaler = sx, methylation_scaler = sy, mca = FactoMineR::MCA(f, ncp = 20, graph = FALSE), 
        fp_features = names(f), fp_levels = lapply(f, levels), train_cells = cells, train_drugs = drugs)
}

project_drug <- function(s, fp, drugs) {
    f <- factor_fp(fp[drugs, s$fp_features, drop = FALSE], s$fp_levels)
    z <- predict(s$mca, newdata = f)$coord
    rownames(z) <- drugs
    stopifnot(ncol(z) == 20, all(is.finite(z)))
    z
}

normal_tissue <- function(x) {
    x <- trimws(x)
    x[is.na(x) | x == "" | toupper(x) == "NS"] <- "__UNKNOWN__"
    x
}

encode_tissue <- function(cells, tissue, lev) {
    z <- outer(unname(tissue[cells]), lev, "==") * 1
    rownames(z) <- cells
    colnames(z) <- paste0("Tissue_", seq_along(lev))
    stopifnot(all(is.finite(z)), all(rowSums(z) %in% 0:1))
    z
}

feature_matrix <- function(pairs, ds, os, indices) {
    p <- cbind(ds[pairs$drug_key, , drop = FALSE], os[pairs$cell_name, indices, drop = FALSE])
    colnames(p) <- c(paste0("MCA_", 1:20), paste0("Omics_", seq_along(indices)))
    p
}

xgb_params <- function(g) {
    c(list(objective = "reg:squarederror", tree_method = "hist", subsample = 1, colsample_bytree = 1, 
        alpha = 0, gamma = 0, nthread = 1, seed = 20260921L), list(max_depth = as.integer(g$max_depth), 
        eta = as.numeric(g$eta), min_child_weight = as.numeric(g$min_child_weight), lambda = as.numeric(g$lambda)))
}

model_route <- function(id) {
    label <- sub("^s[0-9]+_", "", id)
    if (startsWith(id, "s16_")) 
        return(list(method = "O2PLS", indices = switch(label, joint_only = 1:10, joint_expression = 1:15, 
            joint_methylation = c(1:10, 16:20), specific_only = 11:20), tissue = FALSE))
    list(method = label, indices = 1:20, tissue = startsWith(id, "s24_"))
}

objects <- readRDS("results/01_data/omics_from_csv.rds")
read_table <- function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
response <- as.matrix(read.csv("raw_data/Drug_response_IC50.csv", row.names = 1, check.names = FALSE))
fp <- as.matrix(read.csv("raw_data/PubchemFP.csv", row.names = 1, check.names = FALSE))
stopifnot(identical(rownames(response), rownames(objects$expression)), identical(rownames(response), rownames(objects$methylation)), all(fp %in% c(0, 1)), !anyNA(fp))
trace <- read_table("data/internal_response_pair_trace.csv")
stopifnot(identical(trace$cell_name, rep(rownames(response), times = ncol(response))), identical(trace$drug_key, rep(trimws(colnames(response)), each = nrow(response))), identical(trace$response_value, as.vector(response)))
trace$response_value <- as.vector(response)
has_fp <- trace$drug_key %in% trimws(rownames(fp))
keep <- has_fp & !is.na(trace$response_value)
stopifnot(identical(has_fp, trace$fingerprint_matched), all(trace$equal_ln_6dp_candidates[keep] == 1L))
a <- trace[keep, c("pair_id", "key", "cell_name", "cosmic_id", "drug_key", "response_value", "reference_candidates", "reference_drug_id_numeric_trace", "reference_result_id_numeric_trace")]
a$multiple_reference <- a$reference_candidates > 1L
a$target_scale <- "ln(IC50 / 1 micromolar)"
annotation <- read_table("data/source_metadata/cell_annotations.csv")
a$primary_site <- annotation$primary_site[match(a$cell_name, annotation$cell_name)]
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20260922L)
cells <- sort(unique(a$cell_name))
test_cells <- sample(cells, ceiling(length(cells) * 0.20))
cell_assign <- data.frame(cell_name = cells, split = ifelse(cells %in% test_cells, "test", "development"), inner_fold = NA_integer_)
dev <- which(cell_assign$split == "development")
cell_assign$inner_fold[dev] <- sample(rep(seq_len(5L), length.out = length(dev)))
cell_assign$primary_site <- annotation$primary_site[match(cell_assign$cell_name, annotation$cell_name)]
ci <- match(a$cell_name, cell_assign$cell_name)
assignments <- data.frame(pair_id = a$pair_id, split = cell_assign$split[ci], inner_fold = cell_assign$inner_fold[ci])
stopifnot(nrow(a) == 124904L, sum(assignments$split == "development") == 99838L, sum(assignments$split == "test") == 25066L)
actual <- c(objects, list(analysis_pairs = a, cell_pair_split = assignments, cell_split = cell_assign))
expected <- jsonlite::fromJSON("config/reconstructed_object_hashes.json")
checks <- data.frame(object = names(actual), passed = vapply(names(actual), function(n) digest::digest(actual[[n]], algo = "sha256") == expected[[n]], logical(1)))
dir.create("results/01_data", recursive = TRUE, showWarnings = FALSE)
write.csv(checks, "results/01_data/reconstruction_checks.csv", row.names = FALSE)
print(checks)
if (!all(checks$passed)) stop("Reconstructed objects differ from the manuscript inputs. Inspect reconstruction_checks.csv.")
saveRDS(objects, "prepared_inputs/omics_verified_headers.rds")
saveRDS(list(analysis_pairs = a, cell_pair_split = assignments, cell_split = cell_assign), "prepared_inputs/protocol_config.rds")
file.copy("raw_data/PubchemFP.csv", "prepared_inputs/PubchemFP.csv", overwrite = TRUE)
for (config in c("config/internal_input_hashes.json", "config/prepared_input_hashes.json")) {
    hashes <- jsonlite::fromJSON(config)
    for (n in c("protocol_config.rds", "omics_verified_headers.rds", "PubchemFP.csv")) hashes[[n]] <- digest::digest(file = file.path("prepared_inputs", n), algo = "sha256")
    jsonlite::write_json(as.list(hashes), config, pretty = TRUE, auto_unbox = TRUE)
}
write.csv(assignments, "results/01_data/pair_assignments.csv", row.names = FALSE)
write.csv(cell_assign, "results/01_data/cell_assignments.csv", row.names = FALSE)
write.csv(a, "results/01_data/analysis_pairs.csv", row.names = FALSE)
message("PASS: original CSV values, corrected feature names and recreated split match manuscript inputs exactly.")
run_internal(repo, result_dir("input_check"), check_only = TRUE)
run_external_expression(repo, result_dir("external_check"), "MCA", mode = "check")
run_external_reduced(repo, result_dir("external_check"), "O2PLS", mode = "check")
packages <- c("org.Hs.eg.db", "GO.db", "AnnotationDbi",
              "IlluminaHumanMethylation450kanno.ilmn12.hg19")
print(data.frame(package = packages, installed = vapply(packages,
  requireNamespace, logical(1), quietly = TRUE)))
capture.output(sessionInfo(), file = "results/session_info.txt")
