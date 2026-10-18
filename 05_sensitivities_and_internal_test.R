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

compare_internal_refit <- function(repo, out) {
    a <- read.csv(file.path(out, "predictions.csv"), check.names = FALSE)
    b <- read.csv(file.path(repo, "data/internal_predictions.csv"), check.names = FALSE)
    ix <- match(b$pair_id, a$pair_id)
    stopifnot(!anyNA(ix), !anyDuplicated(a$pair_id))
    models <- names(b)[6:ncol(b)]
    z <- do.call(rbind, lapply(models, function(m) data.frame(model = m, max_abs_prediction_difference = max(abs(a[ix, 
        m] - b[[m]])), refit_RMSE = sqrt(mean((a[ix, m] - b$observed)^2)), recorded_RMSE = sqrt(mean((b[[m]] - 
        b$observed)^2)))))
    z$RMSE_difference <- z$refit_RMSE - z$recorded_RMSE
    write.csv(z, file.path(out, "comparison_with_manuscript.csv"), row.names = FALSE)
    print(z)
}

project_ranked_omics <- function(model, x, y) {
    if (is.null(model$rank) || model$rank == 20L) 
        return(project_omics(model, x, y))
    stopifnot(identical(rownames(x), rownames(y)))
    fit <- model$fit
    scores <- switch(model$method, O2PLS = {
        specific_x <- x %*% fit$W_Yosc
        specific_y <- y %*% fit$C_Xosc
        joint <- joint_scores(x - specific_x %*% t(fit$P_Yosc.), fit$W., isTRUE(fit$flags$highd))
        cbind(joint, specific_x, specific_y)
    }, PLS = joint_scores(x, fit$W., isTRUE(fit$flags$highd)), PCA_concat = predict(fit, newdata = cbind(x, 
        y)), PCA_separate = cbind(predict(model$gene, newdata = x), predict(model$meth, newdata = y)), 
        JIVE = t(qr.coef(model$loading_qr, t(cbind(x, y)))))
    rownames(scores) <- rownames(x)
    colnames(scores) <- paste0("Omics_", seq_len(ncol(scores)))
    stopifnot(ncol(scores) == model$rank, all(is.finite(scores)))
    scores
}

fit_ranked_omics <- function(method, x, y, logfile, rank) {
    if (rank == 20L) 
        return(fit_omics(method, x, y, logfile))
    stopifnot(rank %in% c(12L, 28L))
    set.seed(20260921L)
    model <- list(method = method, rank = rank)
    if (method == "PCA_concat") 
        model$fit <- prcomp(cbind(x, y), center = FALSE, scale. = FALSE, rank. = rank)
    if (method == "PCA_separate") {
        model$gene <- prcomp(x, center = FALSE, scale. = FALSE, rank. = rank/2)
        model$meth <- prcomp(y, center = FALSE, scale. = FALSE, rank. = rank/2)
    }
    if (method %in% c("O2PLS", "PLS")) {
        messages <- character()
        model$fit <- withCallingHandlers(OmicsPLS::o2m(x, y, n = if (method == "O2PLS") 
            rank/2
        else rank, nx = if (method == "O2PLS") 
            rank/4
        else 0, ny = if (method == "O2PLS") 
            rank/4
        else 0, stripped = TRUE, tol = 1e-10, max_iterations = 10000), message = function(m) {
            messages <<- c(messages, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
        writeLines(messages, logfile)
        if (isTRUE(model$fit$flags$highd)) {
            lines <- messages[grepl("stopped after", messages, fixed = TRUE)]
            iterations <- as.integer(sub(".*stopped after ([0-9]+) iterations.*", "\\1", lines))
            stopifnot(length(iterations) > 0L, all(iterations < 10000))
        }
    }
    if (method == "JIVE") {
        messages <- capture.output(jive <- r.jive::jive(list(t(x), t(y)), rankJ = rank/2, rankA = c(rank/4, 
            rank/4), method = "given", center = FALSE, scale = FALSE, orthIndiv = TRUE, est = TRUE, 
            maxiter = if (rank == 28L) 
                10000
            else 1000, showProgress = TRUE))
        writeLines(messages, logfile)
        stopifnot(any(grepl("JIVE algorithm converged after", messages, fixed = TRUE)))
        joint <- svd(do.call(rbind, jive$joint), nu = rank/2, nv = 0)$u
        first <- svd(jive$individual[[1]], nu = rank/4, nv = 0)$u
        second <- svd(jive$individual[[2]], nu = rank/4, nv = 0)$u
        model$loadings <- cbind(joint, rbind(first, matrix(0, ncol(y), rank/4)), rbind(matrix(0, 
            ncol(x), rank/4), second))
        model$loading_qr <- qr(model$loadings)
        stopifnot(model$loading_qr$rank == rank)
    }
    scores <- project_ranked_omics(model, x, y)
    stopifnot(qr(scores)$rank == rank)
    model
}

run_cv <- function(repo, out, routes, rank = 20L, single_reference = FALSE) {
    require_model_packages()
    recorded <- jsonlite::fromJSON(file.path(repo, "config/final_internal_choices.json"), simplifyVector = FALSE)
    allowed_routes <- vapply(Filter(function(m) !is.null(m$selected_configuration), recorded$models), function(m) m$id, "")
    if (!all(routes %in% allowed_routes)) stop("Unknown manuscript route: ", paste(setdiff(routes, allowed_routes), collapse = ", "))
    protocol <- readRDS(file.path(repo, "prepared_inputs/protocol_config.rds"))
    omics <- readRDS(file.path(repo, "prepared_inputs/omics_verified_headers.rds"))
    fp <- as.matrix(read.csv(file.path(repo, "prepared_inputs/PubchemFP.csv"), row.names = 1, check.names = FALSE))
    rownames(fp) <- trimws(rownames(fp))
    a <- protocol$cell_pair_split
    pairs <- protocol$analysis_pairs
    stopifnot(identical(a$pair_id, pairs$pair_id))
    keep <- a$split == "development"
    if (single_reference) 
        keep <- keep & !pairs$multiple_reference
    a <- a[keep, ]
    pairs <- pairs[keep, ]
    grid <- read.csv(file.path(repo, "data/tuning_grid.csv"))
    stopifnot(nrow(grid) == 48L, length(unique(a$inner_fold)) == 5L)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    signature <- list(rank = rank, single_reference = single_reference, pair_ids = pairs$pair_id, 
        folds = a$inner_fold, grid = grid, code = digest::digest(deparse(sys.function()), algo = "sha256"), inputs = vapply(file.path(repo, "prepared_inputs", c("protocol_config.rds", 
            "omics_verified_headers.rds", "PubchemFP.csv")), sha, ""), versions = vapply(c("OmicsPLS", 
            "FactoMineR", "xgboost", "r.jive"), function(p) as.character(packageVersion(p)), ""))
    manifest <- file.path(out, "run_specification.rds")
    if (file.exists(manifest)) 
        stopifnot(identical(readRDS(manifest), signature))
    else saveRDS(signature, manifest)
    cached <- function(name, calculate) {
        path <- file.path(out, name)
        if (file.exists(path)) {
            stopifnot(file.exists(paste0(path, ".sha256")), sha(path) == readLines(paste0(path, 
                ".sha256")))
            return(readRDS(path))
        }
        value <- calculate()
        saveRDS(value, path)
        writeLines(sha(path), paste0(path, ".sha256"))
        value
    }
    sums <- setNames(lapply(routes, function(x) numeric(48)), routes)
    counts <- setNames(integer(length(routes)), routes)
    for (fold in 1:5) {
        message("Development fold ", fold, " / 5")
        tr <- pairs[a$inner_fold != fold, ]
        va <- pairs[a$inner_fold == fold, ]
        tc <- sort(unique(tr$cell_name))
        vc <- sort(unique(va$cell_name))
        td <- sort(unique(tr$drug_key))
        vd <- sort(unique(va$drug_key))
        stopifnot(!any(tc %in% vc), all(vd %in% td))
        shared <- cached(paste0("fold", fold, "_shared.rds"), function() fit_shared(omics, fp, tc, 
            td))
        x <- scale_apply(omics$expression[tc, , drop = FALSE], shared$expression_scaler)
        y <- scale_apply(omics$methylation[tc, , drop = FALSE], shared$methylation_scaler)
        xv <- scale_apply(omics$expression[vc, , drop = FALSE], shared$expression_scaler)
        yv <- scale_apply(omics$methylation[vc, , drop = FALSE], shared$methylation_scaler)
        drt <- project_drug(shared, fp, td)
        drv <- project_drug(shared, fp, vd)
        for (route in routes) {
            info <- model_route(route)
            method <- info$method
            message("  ", route)
            model <- cached(paste0("fold", fold, "_", method, "_transform.rds"), function() fit_ranked_omics(method, 
                x, y, file.path(out, paste0("fold", fold, "_", method, ".log")), rank))
            st <- project_ranked_omics(model, x, y)
            sv <- project_ranked_omics(model, xv, yv)
            indices <- if (rank == 20L) 
                info$indices
            else seq_len(rank)
            pt <- feature_matrix(tr, drt, st, indices)
            pv <- feature_matrix(va, drv, sv, indices)
            sc <- scale_fit(pt, FALSE)
            pt <- scale_apply(pt, sc)
            pv <- scale_apply(pv, sc)
            if (info$tissue) {
                meta <- unique(pairs[, c("cell_name", "primary_site")])
                tissue <- setNames(normal_tissue(meta$primary_site), meta$cell_name)
                levels <- sort(setdiff(unique(tissue[tc]), "__UNKNOWN__"))
                pt <- cbind(pt, encode_tissue(tr$cell_name, tissue, levels))
                pv <- cbind(pv, encode_tissue(va$cell_name, tissue, levels))
            }
            saveRDS(sc, file.path(out, paste0("fold", fold, "_", route, "_component_scaler.rds")))
            dt <- xgboost::xgb.DMatrix(pt, label = tr$response_value, nthread = 1)
            dv <- xgboost::xgb.DMatrix(pv, nthread = 1)
            for (base in unique(grid$base_id)) {
                rows <- grid[grid$base_id == base, ]
                name <- paste0("fold", fold, "_", route, "_base", base, ".rds")
                value <- cached(name, function() {
                  set.seed(20260921L)
                  booster <- xgboost::xgb.train(params = xgb_params(rows[1, ]), data = dt, nrounds = max(rows$nrounds), 
                    verbose = 0)
                  model_path <- file.path(out, paste0("fold", fold, "_", route, "_base", base, ".json"))
                  xgboost::xgb.save(booster, model_path)
                  writeLines(sha(model_path), paste0(model_path, ".sha256"))
                  predictions <- vapply(rows$nrounds, function(n) as.numeric(predict(booster, dv, 
                    iterationrange = c(1L, as.integer(n)))), numeric(nrow(va)))
                  colnames(predictions) <- rows$config_id
                  stopifnot(all(is.finite(predictions)))
                  list(pair_id = va$pair_id, observed = va$response_value, predictions = predictions)
                })
                stopifnot(identical(value$pair_id, va$pair_id))
                ix <- match(rows$config_id, grid$config_id)
                sums[[route]][ix] <- sums[[route]][ix] + colSums((value$predictions - va$response_value)^2)
            }
            counts[route] <- counts[route] + nrow(va)
            rm(model, st, sv, pt, pv, dt, dv)
            gc(FALSE)
        }
    }
    selections <- list()
    for (route in routes) {
        stopifnot(counts[route] == nrow(pairs))
        summary <- cbind(model = route, grid, RMSE = sqrt(sums[[route]]/counts[route]))
        write.csv(summary, file.path(out, paste0(route, "_cv_metrics.csv")), row.names = FALSE)
        best <- summary[summary$RMSE <= min(summary$RMSE) + 1e-08, ]
        best <- best[order(best$max_depth, best$nrounds, best$config_id), ][1, ]
        selections[[route]] <- best
        parts <- lapply(1:5, function(fold) {
            path <- paste0("fold", fold, "_", route, "_base", best$base_id, ".rds")
            z <- cached(path, function() stop("Missing CV result"))
            p <- pairs[match(z$pair_id, pairs$pair_id), ]
            data.frame(pair_id = z$pair_id, cell = p$cell_name, drug = p$drug_key, fold = fold, 
                observed = z$observed, prediction = z$predictions[, best$config_id])
        })
        oof <- do.call(rbind, parts)
        stopifnot(!anyDuplicated(oof$pair_id), nrow(oof) == nrow(pairs), abs(sqrt(mean((oof$observed - 
            oof$prediction)^2)) - best$RMSE) < 1e-10)
        write.csv(oof, file.path(out, paste0(route, "_selected_oof_predictions.csv")), row.names = FALSE)
    }
    path <- file.path(out, "selected_configurations.csv")
    previous <- if (file.exists(path)) 
        read.csv(path)
    else NULL
    if (!is.null(previous)) 
        previous <- previous[!previous$model %in% routes, ]
    selected <- rbind(previous, do.call(rbind, selections))
    write.csv(selected, path, row.names = FALSE)
    spec <- jsonlite::fromJSON(file.path(repo, "config/final_internal_choices.json"), simplifyVector = FALSE)
    comparison <- do.call(rbind, lapply(selections, function(s) {
        original <- Filter(function(m) m$id == s$model, spec$models)[[1]]$selected_configuration
        data.frame(model = s$model, rerun_config = s$config_id, recorded_config = original$config_id, 
            same_configuration = s$config_id == original$config_id)
    }))
    write.csv(comparison, file.path(out, paste0("selection_comparison_", paste(routes, collapse = "_"), 
        ".csv")), row.names = FALSE)
    capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
    print(do.call(rbind, selections))
    invisible(selected)
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

run_permutation <- function(repo, cv_dir, out) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    selected <- read.csv(file.path(cv_dir, "selected_configurations.csv"))
    selected <- selected[selected$model == "s15_O2PLS", ]
    stopifnot(nrow(selected) == 1L)
    reference <- read.csv(file.path(cv_dir, "s15_O2PLS_selected_oof_predictions.csv"))
    omics <- readRDS(file.path(repo, "prepared_inputs/omics_verified_headers.rds"))
    groups <- list(joint = 1:10, expression_specific = 11:15, methylation_specific = 16:20)
    results <- list()
    maps <- list()
    for (fold in 1:5) {
        message("Permutation fold ", fold, " / 5")
        shared <- readRDS(file.path(cv_dir, paste0("fold", fold, "_shared.rds")))
        model <- readRDS(file.path(cv_dir, paste0("fold", fold, "_O2PLS_transform.rds")))
        scaler <- readRDS(file.path(cv_dir, paste0("fold", fold, "_s15_O2PLS_component_scaler.rds")))
        d <- reference[reference$fold == fold, ]
        cells <- sort(unique(d$cell))
        x <- scale_apply(omics$expression[cells, , drop = FALSE], shared$expression_scaler)
        y <- scale_apply(omics$methylation[cells, , drop = FALSE], shared$methylation_scaler)
        scores <- project_omics(model, x, y)
        model_path <- file.path(cv_dir, paste0("fold", fold, "_s15_O2PLS_base", selected$base_id, 
            ".json"))
        stopifnot(sha(model_path) == readLines(paste0(model_path, ".sha256")))
        booster <- xgboost::xgb.load(model_path)
        predict_scores <- function(s) {
            features <- cbind(shared$mca$ind$coord[d$drug, , drop = FALSE], s[d$cell, , drop = FALSE])
            colnames(features) <- c(paste0("MCA_", 1:20), paste0("Omics_", 1:20))
            as.numeric(predict(booster, xgboost::xgb.DMatrix(scale_apply(features, scaler), nthread = 1), 
                iterationrange = c(1L, as.integer(selected$nrounds))))
        }
        stopifnot(max(abs(predict_scores(scores) - d$prediction)) < 1e-06)
        for (repetition in 1:50) {
            set.seed(20260922L + fold * 1000L + repetition)
            order <- sample.int(length(cells))
            maps[[length(maps) + 1L]] <- data.frame(fold, repetition, recipient_cell = cells, donor_cell = cells[order])
            for (group in names(groups)) {
                permuted <- scores
                permuted[, groups[[group]]] <- scores[order, groups[[group]], drop = FALSE]
                p <- predict_scores(permuted)
                results[[length(results) + 1L]] <- data.frame(fold, repetition, group, n = nrow(d), 
                  SSE = sum((p - d$observed)^2))
            }
        }
    }
    results <- do.call(rbind, results)
    pooled <- aggregate(cbind(SSE, n) ~ group + repetition, results, sum)
    pooled$RMSE <- sqrt(pooled$SSE/pooled$n)
    pooled$delta_RMSE <- pooled$RMSE - sqrt(mean((reference$observed - reference$prediction)^2))
    summary <- do.call(rbind, lapply(split(pooled, pooled$group), function(d) data.frame(group = d$group[1], 
        mean_delta_RMSE = mean(d$delta_RMSE), sd_delta_RMSE = sd(d$delta_RMSE), q025 = unname(quantile(d$delta_RMSE, 
            0.025)), q975 = unname(quantile(d$delta_RMSE, 0.975)))))
    write.csv(results, file.path(out, "permutation_fold_metrics.csv"), row.names = FALSE)
    write.csv(do.call(rbind, maps), file.path(out, "cell_permutation_maps.csv"), row.names = FALSE)
    write.csv(pooled, file.path(out, "permutation_pooled_metrics.csv"), row.names = FALSE)
    write.csv(summary, file.path(out, "permutation_summary.csv"), row.names = FALSE)
    print(summary)
}

run_registry_sensitivity <- function(repo, old, root) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    a <- read.csv(file.path(old, "predictions.csv"), check.names = FALSE)
    id <- read.csv(file.path(repo, "data/internal_test_relatedness_audit.csv"))
    excluded <- id$cell[tolower(as.character(id$related_to_development)) == "true"]
    original_n <- nrow(a)
    a <- a[!a$cell %in% excluded, ]
    groups <- setNames(id$registry_component, id$cell)
    stopifnot(!anyNA(groups[a$cell]))
    models <- setdiff(names(a), c("pair_id", "cell", "drug", "tissue", "observed"))
    err <- sweep(as.matrix(a[, models]), 1, a$observed, "-")
    sq <- err^2
    met <- data.frame(model = models, n = nrow(a), cells = length(unique(a$cell)), RMSE = sqrt(colMeans(sq)), 
        MAE = colMeans(abs(err)), R2 = 1 - colSums(sq)/sum((a$observed - mean(a$observed))^2))
    write.csv(met, file.path(root, "internal_registry_sensitivity_metrics.csv"), row.names = FALSE)
    ss <- rowsum(sq, groups[a$cell])
    nn <- table(groups[a$cell])[rownames(ss)]
    spec <- read.csv(file.path(old, "primary_contrasts.csv"))
    B <- 20000L
    boot <- matrix(0, B, nrow(spec))
    set.seed(20260923)
    for (i in 1:B) {
        w <- tabulate(sample.int(nrow(ss), nrow(ss), replace = TRUE), nbins = nrow(ss))
        rmse <- sqrt(colSums(ss * w)/sum(nn * w))
        boot[i, ] <- rmse[spec$left] - rmse[spec$right]
    }
    out <- data.frame(contrast = spec$contrast, left = spec$left, right = spec$right, delta_RMSE = met$RMSE[match(spec$left, 
        met$model)] - met$RMSE[match(spec$right, met$model)], lower_adjusted = apply(boot, 2, quantile, 
        probs = 0.05/12), upper_adjusted = apply(boot, 2, quantile, probs = 1 - 0.05/12), level = 1 - 
        0.05/6, B = B, posthoc = TRUE)
    write.csv(out, file.path(root, "internal_registry_sensitivity_contrasts.csv"), row.names = FALSE)
    jsonlite::write_json(list(posthoc = TRUE, excluded_cells = excluded, removed_pairs = original_n - 
        nrow(a), remaining_cells = length(unique(a$cell)), remaining_pairs = nrow(a), registry_components = nrow(ss), 
        all_adjusted_intervals_include_zero = all(out$lower_adjusted <= 0 & out$upper_adjusted >= 
            0), original_frozen_results_unchanged = TRUE, interpretation = "Conservative transitive registry graph may connect hybrid parents; not proof of shared donor"), 
        file.path(root, "internal_registry_sensitivity_summary.json"), pretty = TRUE, auto_unbox = TRUE)
    print(met[1:5, ])
    print(out)
}

routes <- c("s15_O2PLS", "s15_PCA_concat", "s15_PCA_separate", "s15_PLS", "s15_JIVE")
run_cv(repo, result_dir("05_source"), routes, single_reference = TRUE)
run_cv(repo, result_dir("05_rank12"), routes, rank = 12L)
run_cv(repo, result_dir("05_rank28"), routes, rank = 28L)

run_internal(repo, result_dir("06_internal"), check_only = FALSE)
compare_internal_refit(repo, result_dir("06_internal"))

run_permutation(repo, result_dir("02_cv"), result_dir("14_permutation"))

run_registry_sensitivity(repo, result_dir("06_internal"), result_dir("13_registry"))
