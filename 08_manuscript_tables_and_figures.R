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

run_fresh_report <- function(repo, out) {
    input_files <- c(internal_predictions.csv = "06_internal/predictions.csv", internal_contrasts.csv = "06_internal/primary_contrasts.csv", 
        external_expression_predictions.csv = "07_external/expression_refit_predictions_with_targets.csv", 
        external_reduced_predictions.csv = "07_external/reduced_refit_predictions_with_targets.csv", 
        feature_stability.csv = "08_features/feature_stability.csv", case_study_associations.csv = "09_biology/case_study_associations.csv", 
        case_study_data.csv = "09_biology/case_study_data.csv", probe_gene_annotation.csv = "09_biology/probe_gene_annotation.csv", 
        expression_gene_mapping.csv = "09_expression_GO/gene_identifier_mapping.csv")
    paths <- file.path(repo, "results", input_files)
    if (!all(file.exists(paths))) 
        stop("Run the preceding stages. Missing: ", paste(input_files[!file.exists(paths)], collapse = ", "))
    if (dir.exists(out)) 
        stop("Choose a new report directory.")
    dir.create(out, recursive = TRUE)
    staging <- file.path(out, "report_inputs")
    dir.create(staging)
    source_files <- list.files(file.path(repo, "data"), full.names = TRUE)
    source_files <- source_files[!file.info(source_files)$isdir]
    stopifnot(all(file.copy(source_files, staging)))
    stopifnot(all(file.copy(paths, file.path(staging, names(input_files)), overwrite = TRUE)))
    run_report(out, data_dir = staging, strict_reference = FALSE)
    write.csv(data.frame(report_input = names(input_files), source = input_files, sha256 = vapply(paths, 
        sha, "")), file.path(out, "fresh_result_sources.csv"), row.names = FALSE)
    writeLines(c("Predictive tables and Figures 2, 3, 5 and 6 use the new fitted predictions.", 
        "Table 4 and Figure 4 use the rerun biological calculations.", "Table S5 uses newly fitted fold loadings and rerun annotation mapping.", 
        "Cohort identities, targets and split references remain the frozen study inputs.", "Agreement with the manuscript is reported separately; it is not guaranteed."), 
        file.path(out, "RESULT_SOURCE.txt"))
}

sha <- function(p) digest::digest(file = p, algo = "sha256")

run_report <- function(out, data_dir = "data", strict_reference = TRUE) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    rd <- function(n) read.csv(file.path(data_dir, n), check.names = FALSE, stringsAsFactors = FALSE)
    wc <- function(d, n) write.csv(d, file.path(out, n), row.names = FALSE, na = "")
    metrics <- function(y, p) c(n = length(y), RMSE = sqrt(mean((y - p)^2)), MAE = mean(abs(y - 
        p)), R2 = if (var(y) > 0) 1 - sum((y - p)^2)/sum((y - mean(y))^2) else NA, Pearson = if (var(p) > 
        0 && var(y) > 0) cor(y, p) else NA, Spearman = if (var(p) > 0 && var(y) > 0) cor(y, p, method = "spearman") else NA)
    perdrug <- function(d, models, y, drug) do.call(rbind, lapply(models, function(m) do.call(rbind, 
        lapply(sort(unique(d[[drug]])), function(k) {
            a <- d[d[[drug]] == k, ]
            data.frame(model = m, drug = k, as.list(metrics(a[[y]], a[[m]])))
        }))))
    pooled <- function(d, models, y, pd) do.call(rbind, lapply(models, function(m) data.frame(model = m, 
        as.list(metrics(d[[y]], d[[m]])), macro_drug_RMSE = mean(pd$RMSE[pd$model == m]))))
    fig <- function(name, w, h, fun) {
        pdf(file.path(out, paste0(name, ".pdf")), width = w, height = h)
        fun()
        dev.off()
        png(file.path(out, paste0(name, ".png")), width = w * 240, height = h * 240, res = 240)
        fun()
        dev.off()
    }
    d <- rd("internal_predictions.csv")
    models <- names(d)[6:ncol(d)]
    stopifnot(nrow(d) == 25066, length(unique(d$cell)) == 147, length(unique(d$drug)) == 176, !anyDuplicated(d$pair_id), 
        all(is.finite(as.matrix(d[, c("observed", models)]))))
    pd <- perdrug(d, models, "observed", "drug")
    pm <- pooled(d, models, "observed", pd)
    ref <- rd("internal_pooled_reference.csv")
    delta <- max(abs(as.matrix(pm[, c("RMSE", "MAE", "R2", "macro_drug_RMSE")]) - as.matrix(ref[match(pm$model, 
        ref$model), c("RMSE", "MAE", "R2", "macro_drug_RMSE")])))
    if (strict_reference) 
        stopifnot(delta < 1e-10)
    pr <- rd("internal_per_drug_reference.csv")
    ii <- match(paste(pd$model, pd$drug), paste(pr$model, pr$drug))
    delta_drug <- max(abs(pd$RMSE - pr$RMSE[ii]))
    if (strict_reference) 
        stopifnot(delta_drug < 1e-10)
    wc(pm, "Table2_internal_metrics.csv")
    wc(pd, "internal_per_drug_metrics.csv")
    sel <- c("s15_O2PLS", "s23_PCA_expression", "s15_PCA_concat")
    labs <- c("O2PLS", "Expression PCA", "Concatenated PCA")
    colors <- c("#0072B2", "#009E73", "#D55E00")
    z <- pd[pd$model %in% sel, ]
    b <- pd[pd$model == "drug_mean", ]
    z$delta_vs_drug_mean <- z$RMSE - b$RMSE[match(z$drug, b$drug)]
    b <- pd[pd$model == "drug_tissue_shrinkage10", ]
    z$delta_vs_drug_tissue <- z$RMSE - b$RMSE[match(z$drug, b$drug)]
    z$spearman_eligible <- z$n >= 30 & is.finite(z$Spearman)
    summary <- do.call(rbind, lapply(seq_along(sel), function(i) {
        a <- z[z$model == sel[i], ]
        v <- a$Spearman[a$spearman_eligible]
        data.frame(model = sel[i], label = labs[i], drugs = nrow(a), n_min = min(a$n), n_max = max(a$n), 
            median_delta_drug_mean = median(a$delta_vs_drug_mean), q25_delta_drug_mean = quantile(a$delta_vs_drug_mean, 
                0.25), q75_delta_drug_mean = quantile(a$delta_vs_drug_mean, 0.75), improved_drug_mean = sum(a$delta_vs_drug_mean < 
                -1e-08), improved_drug_tissue = sum(a$delta_vs_drug_tissue < -1e-08), tied_drug_tissue = sum(abs(a$delta_vs_drug_tissue) <= 
                1e-08), spearman_drugs = length(v), median_spearman = median(v), q25_spearman = quantile(v, 
                0.25), q75_spearman = quantile(v, 0.75))
    }))
    wc(z, "TableS6_within_drug_detail.csv")
    wc(summary, "TableS7_within_drug_summary.csv")
    fig("Figure2_within_drug", 10, 8, function() {
        par(mfrow = c(2, 1), mar = c(4, 12, 3, 2), las = 1)
        boxplot(lapply(sel, function(m) z$delta_vs_drug_mean[z$model == m]), horizontal = TRUE, 
            names = labs, col = adjustcolor(colors, 0.3), outline = FALSE, xlab = "Per-drug RMSE difference from drug-mean baseline", 
            main = "A  Error change across 176 drugs", ylim = range(z$delta_vs_drug_mean))
        abline(v = 0, lty = 2, col = "gray40")
        for (i in 1:3) {
            v <- z$delta_vs_drug_mean[z$model == sel[i]]
            points(v, i + seq(-0.16, 0.16, length.out = length(v)), pch = 16, cex = 0.38, col = adjustcolor(colors[i], 
                0.5))
        }
        boxplot(lapply(sel, function(m) z$Spearman[z$model == m & z$spearman_eligible]), horizontal = TRUE, 
            names = paste0(labs, "\n", summary$improved_drug_tissue, "/176 improved"), col = adjustcolor(colors, 
                0.3), outline = FALSE, xlab = "Within-drug Spearman correlation", main = "B  Cell ranking within each drug", 
            ylim = c(-1, 1))
        abline(v = 0, lty = 2, col = "gray40")
        for (i in 1:3) {
            v <- z$Spearman[z$model == sel[i] & z$spearman_eligible]
            points(v, i + seq(-0.16, 0.16, length.out = length(v)), pch = 16, cex = 0.38, col = adjustcolor(colors[i], 
                0.5))
        }
    })
    for (kind in c("expression", "reduced")) {
        e <- rd(paste0("external_", kind, "_predictions.csv"))
        ms <- if (kind == "expression") 
            c("MCA", "RAW", "DRUG_MEAN")
        else c("O2PLS", "PCA_concat", "Expression", "DRUG_MEAN")
        ep <- perdrug(e, ms, "response_value", "drug_key")
        em <- pooled(e, ms, "response_value", ep)
        wc(ep, paste0("external_", kind, "_per_drug.csv"))
        wc(em, if (kind == "expression") 
            "Table5_external_expression.csv"
        else "Table6_external_reduced.csv")
        fig(if (kind == "expression") 
            "Figure5_external_expression"
        else "Figure6_external_reduced", 10, 5.5, function() {
            par(mar = c(7, 4, 3, 1))
            mat <- sapply(ms, function(m) ep$RMSE[ep$model == m])
            rownames(mat) <- ep$drug[ep$model == ms[1]]
            barplot(t(mat), beside = TRUE, names.arg = rownames(mat), las = 2, col = c(colors, "gray65")[seq_along(ms)], 
                ylab = "Within-drug RMSE", main = if (kind == "expression") 
                  "Expression-only CCLE evaluation: 122 cells"
                else "Reduced-panel CCLE evaluation: 80 cells", ylim = c(0, max(mat) * 1.3))
            legend("topright", ms, fill = c(colors, "gray65")[seq_along(ms)], bty = "n", cex = 0.8)
        })
    }
    c <- rd("internal_contrasts.csv")
    wc(c, "Table3_internal_contrasts.csv")
    fig("Figure3_internal_contrasts", 10, 5, function() {
        par(mar = c(4, 5, 3, 2))
        plot(c$delta_RMSE, 6:1, xlim = range(c$lower_adjusted, c$upper_adjusted) * 1.2, ylim = c(0.5, 
            6.5), yaxt = "n", pch = 19, col = colors[1], ylab = "", xlab = "Pooled RMSE difference (left minus right)", 
            main = "Frozen contrasts with 99.1667% paired cell-bootstrap intervals")
        segments(c$lower_adjusted, 6:1, c$upper_adjusted, 6:1, lwd = 2, col = colors[1])
        axis(2, at = 6:1, labels = c$contrast, las = 1)
        abline(v = 0, lty = 2)
    })
    bio <- rd("case_study_data.csv")
    wc(rd("case_study_associations.csv"), "Table4_case_studies.csv")
    fig("Figure4_biology", 10, 4.5, function() {
        par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
        for (g in c("MGMT", "SLFN11")) {
            a <- bio[bio$gene == g & is.finite(bio$expression) & is.finite(bio$response), ]
            a <- a[!duplicated(a$cell), ]
            plot(a$expression, a$response, pch = 16, cex = 0.55, col = adjustcolor(colors[1], 0.3), 
                xlab = paste(g, "expression"), ylab = "ln(IC50 / 1 micromolar)", main = paste(g, 
                  unique(a$drug), paste0("n=", nrow(a))))
            abline(lm(response ~ expression, a), col = colors[3], lwd = 2)
        }
    })
    f <- rd("feature_stability.csv")
    s <- f[f$stable_top5pct, ]
    a <- rd("probe_gene_annotation.csv")
    g <- rd("expression_gene_mapping.csv")
    s$loading_matrix <- c(expression_joint = "W.", methylation_joint = "C.", expression_specific = "P_Yosc.", 
        methylation_specific = "P_Xosc.")[s$group]
    s$gene_mapping <- vapply(seq_len(nrow(s)), function(i) {
        if (grepl("expression", s$group[i])) {
            v <- g$current_symbol[match(s$feature[i], g$feature)]
            if (is.na(v)) 
                ""
            else v
        }
        else paste(sort(unique(a$gene[a$CpG == s$feature[i]])), collapse = ";")
    }, "")
    s$mapping_status <- ifelse(grepl("expression", s$group), g$status[match(s$feature, g$feature)], 
        ifelse(nchar(s$gene_mapping) > 0, "450k_hg19_annotation", "no_gene_annotation"))
    wc(s, "TableS5_stable_features.csv")
    stopifnot(all(s$count_top5pct >= 4))
    if (strict_reference) 
        stopifnot(sum(s$group == "methylation_joint") == 249, sum(s$group == "methylation_specific") == 
            270)
    wc(unique(f[, c("group", "feature")]), "retained_feature_lists.csv")
    wc(data.frame(drug = sort(unique(d$drug))), "retained_drugs.csv")
    splits <- rd("cell_splits.csv")
    cells <- sort(unique(splits$cell_name))
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(20260922)
    test <- sample(cells, ceiling(length(cells) * 0.2))
    dev <- which(!cells %in% test)
    fold <- rep(NA_integer_, length(cells))
    fold[dev] <- sample(rep(1:5, length.out = length(dev)))
    j <- match(cells, splits$cell_name)
    stopifnot(identical(sort(test), sort(splits$cell_name[splits$split == "test"])), identical(fold[dev], 
        as.integer(splits$inner_fold[j][dev])))
    wc(data.frame(check = c("pooled_metrics_max_abs_difference", "per_drug_RMSE_max_abs_difference", 
        "split_and_folds_recreated", "stable_features"), value = c(delta, delta_drug, 1, nrow(s))), 
        "verification.csv")
    capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
    print(summary)
    e1 <- rd("external_expression_predictions.csv")
    e2 <- rd("external_reduced_predictions.csv")
    sp <- rd("cell_splits.csv")
    pa <- rd("internal_pair_assignments.csv")
    wc(data.frame(cohort = c("Development", "Internal test", "CCLE expression", "CCLE matched-site"), 
        cells = c(587, 147, length(unique(e1$cell_name)), length(unique(e2$cell_name))), pairs = c(99838, 
            nrow(d), nrow(e1), nrow(e2)), expression = c(8046, 8046, 7276, 7276), CpGs = c(8473, 
            8473, NA, 941)), "Table1_cohorts.csv")
    x <- read.csv(file.path(out, "external_expression_per_drug.csv"))
    y <- read.csv(file.path(out, "external_reduced_per_drug.csv"))
    tab <- data.frame(drug = sort(unique(x$drug)))
    for (m in c("MCA", "RAW")) {
        a <- x[x$model == m, ]
        tab[[m]] <- a$RMSE[match(tab$drug, a$drug)]
        tab$n_E <- a$n[match(tab$drug, a$drug)]
    }
    for (m in c("O2PLS", "PCA_concat", "Expression")) {
        a <- y[y$model == m, ]
        tab[[m]] <- a$RMSE[match(tab$drug, a$drug)]
        tab$n_M <- a$n[match(tab$drug, a$drug)]
    }
    wc(tab, "Table7_external_per_drug.csv")
    ca <- rd("case_study_associations.csv")
    rr <- numeric(nrow(ca))
    pp <- rr
    for (i in seq_len(nrow(ca))) {
        a <- bio[bio$gene == ca$gene[i], ]
        if (ca$association[i] == "expression_response") 
            a <- a[!duplicated(a$cell), ]
        if (ca$association[i] == "methylation_expression") {
            a <- rd("SLFN11_development_omics.csv")
            a$expression <- a$SLFN11
            a$methylation <- a$cg14380270
        }
        v <- switch(ca$association[i], expression_response = c("expression", "response"), methylation_response = c("methylation", 
            "response"), methylation_expression = c("methylation", "expression"))
        a <- a[complete.cases(a[, c(v, "tissue")]), ]
        D <- model.matrix(~factor(tissue), a)
        u <- lm.fit(D, rank(a[[v[1]]], ties.method = "average"))$residuals
        v2 <- lm.fit(D, rank(a[[v[2]]], ties.method = "average"))$residuals
        rr[i] <- cor(u, v2)
        df <- nrow(a) - qr(D)$rank - 1
        pp[i] <- 2 * pt(-abs(rr[i] * sqrt(df/(1 - rr[i]^2))), df)
    }
    stopifnot(max(abs(rr - ca$tissue_partial_rank_r)) < 1e-10, max(abs(p.adjust(pp, "BH") - ca$q_BH_case_family)) < 
        1e-10)
    if (strict_reference) 
        for (kind in c("expression", "reduced")) {
            a <- read.csv(file.path(out, if (kind == "expression") 
                "Table5_external_expression.csv"
            else "Table6_external_reduced.csv"))
            b <- rd(paste0("external_", kind, "_metrics_reference.csv"))
            ix <- match(a$model, b$model)
            stopifnot(!anyNA(ix), max(abs(a$RMSE - b$RMSE[ix])) < 1e-10, max(abs(a$macro_drug_RMSE - 
                b$mean_drug_RMSE[ix])) < 1e-10)
        }
    wc(data.frame(check = c("split_and_folds", "case_study_correlations_and_BH", "stable_selection_threshold"), 
        passed = TRUE), "final_checks.csv")
    writeLines(if (strict_reference) 
        "Archived predictions: reference equality checks enforced."
    else "Fresh refit predictions: reference discrepancies reported, not forced to zero.", file.path(out, 
        "RESULT_SOURCE.txt"))
    invisible(out)
}

run_fresh_report(repo, result_dir("10_report"))

expected <- read.csv("config/manuscript_expected_values.csv", stringsAsFactors = FALSE)
stopifnot(nrow(expected) == 113L)
report_dir <- result_dir("10_report")
files <- c(Table2 = "Table2_internal_metrics.csv", Table3 = "Table3_internal_contrasts.csv", Table5_external_expression.csv = "Table5_external_expression.csv", Table6_external_reduced.csv = "Table6_external_reduced.csv", Table7 = "Table7_external_per_drug.csv")
checks <- expected
checks$recalculated <- NA_real_
for (i in seq_len(nrow(expected))) {
    tab <- read.csv(file.path(report_dir, files[[expected$table[i]]]), stringsAsFactors = FALSE, check.names = FALSE)
    key <- if (expected$table[i] == "Table3") "contrast" else if (expected$table[i] == "Table7") "drug" else "model"
    index <- which(tab[[key]] == expected$row[i])
    stopifnot(length(index) == 1L, expected$metric[i] %in% names(tab))
    checks$recalculated[i] <- tab[[expected$metric[i]]][index]
}
checks$passed <- is.finite(checks$recalculated) & abs(checks$recalculated - checks$manuscript) < 0.00005 + 1e-12
write.csv(checks, file.path(report_dir, "manuscript_numeric_checks.csv"), row.names = FALSE)
print(table(checks$passed))
if (!all(checks$passed)) stop("Some recalculated manuscript values differ. Inspect manuscript_numeric_checks.csv.")
message("PASS: 113 published numeric entries match at four-decimal rounding.")
