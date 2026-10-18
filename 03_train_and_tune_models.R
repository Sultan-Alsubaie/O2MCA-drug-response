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

run_cv(repo, result_dir("02_cv"), routes = c("s15_O2PLS", "s15_PCA_concat",
  "s15_PCA_separate", "s15_PLS", "s15_JIVE"))
