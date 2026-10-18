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

run_biology <- function(repo, out, stability_path = file.path(repo, "data/feature_stability.csv")) {
    suppressPackageStartupMessages({
        library(digest)
        library(jsonlite)
        library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        library(org.Hs.eg.db)
        library(GO.db)
        library(AnnotationDbi)
    })
    stopifnot(as.character(packageVersion("org.Hs.eg.db")) == "3.22.0", as.character(packageVersion("GO.db")) == 
        "3.22.0", as.character(packageVersion("IlluminaHumanMethylation450kanno.ilmn12.hg19")) == 
        "0.6.1")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    wc <- function(x, n) write.csv(x, file.path(out, n), row.names = FALSE, na = "NA")
    pp <- file.path(repo, "prepared_inputs/protocol_config.rds")
    p <- readRDS(pp)
    stopifnot(digest(file = file.path(repo, "data/internal_response_pair_trace.csv"), algo = "sha256") == 
        p$config$input_sha256[["trace"]])
    op <- file.path(repo, "prepared_inputs/omics_verified_headers.rds")
    stopifnot(digest(file = op, algo = "sha256") == p$config$input_sha256[["omics"]])
    o <- readRDS(op)
    dev <- sort(unique(p$analysis_pairs$cell_name[p$cell_pair_split$split == "development"]))
    stopifnot(length(dev) == 587)
    ap <- p$analysis_pairs[p$cell_pair_split$split == "development", ]
    meta <- unique(ap[, c("cell_name", "primary_site")])
    stopifnot(!anyDuplicated(meta$cell_name))
    tissue <- setNames(meta$primary_site, meta$cell_name)
    tissue[is.na(tissue) | tissue == ""] <- "UNKNOWN"
    x <- o$expression[dev, , drop = FALSE]
    m <- o$methylation[dev, , drop = FALSE]
    aa <- get("Other", asNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19"))
    ii <- match(colnames(m), rownames(aa))
    a <- as.data.frame(aa[ii[!is.na(ii)], c("UCSC_RefGene_Name", "UCSC_RefGene_Group")])
    ll <- get("Locations", asNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19"))
    loc <- as.data.frame(ll[match(rownames(a), rownames(ll)), ])
    message("Loaded annotation subset")
    mapping <- do.call(rbind, lapply(colnames(m), function(cg) {
        if (!cg %in% rownames(a)) 
            return(NULL)
        g <- strsplit(a[cg, "UCSC_RefGene_Name"], ";", fixed = TRUE)[[1]]
        rg <- strsplit(a[cg, "UCSC_RefGene_Group"], ";", fixed = TRUE)[[1]]
        if (!length(g)) 
            return(NULL)
        data.frame(CpG = cg, gene = g, region = if (length(rg) == length(g)) 
            rg
        else NA_character_)
    }))
    mapping <- unique(mapping[nzchar(mapping$gene), ])
    mapping$chr <- loc[mapping$CpG, "chr"]
    mapping$position_hg19 <- loc[mapping$CpG, "pos"]
    mapping$expression_measured <- mapping$gene %in% colnames(x)
    wc(mapping, "probe_gene_annotation.csv")
    wc(data.frame(CpG = colnames(m), annotation_matched = colnames(m) %in% rownames(a)), "probe_annotation_coverage.csv")
    wc(data.frame(cell = dev, tissue = unname(tissue[dev])), "biology_development_cells.csv")
    assoc <- function(u, v, ti) {
        ok <- is.finite(u) & is.finite(v) & !is.na(ti)
        u <- u[ok]
        v <- v[ok]
        ti <- droplevels(factor(ti[ok]))
        n <- length(u)
        D <- if (nlevels(ti) > 1) 
            model.matrix(~ti)
        else matrix(1, n, 1)
        q <- qr(D)
        df <- n - q$rank - 1
        r <- pv <- raw <- NA_real_
        if (n > 3 && sd(u) > 0 && sd(v) > 0 && df > 2) {
            ru <- rank(u)
            rv <- rank(v)
            raw <- cor(ru, rv)
            e1 <- qr.resid(q, ru)
            e2 <- qr.resid(q, rv)
            if (sum(e1^2) > 1e-12 && sum(e2^2) > 1e-12) {
                r <- sum(e1 * e2)/sqrt(sum(e1^2) * sum(e2^2))
                r <- max(-1, min(1, r))
                pv <- 2 * pt(-abs(r * sqrt(df/max(1e-15, 1 - r * r))), df)
            }
        }
        data.frame(n = n, tissue_design_rank = q$rank, df = df, spearman = raw, tissue_partial_rank_r = r, 
            p_value = pv)
    }
    cis <- unique(mapping[mapping$expression_measured, c("CpG", "gene")])
    message("cis tests: ", nrow(cis))
    D <- model.matrix(~factor(tissue[dev]))
    q <- qr(D)
    xr <- qr.resid(q, apply(x, 2, rank))
    mr <- qr.resid(q, apply(m, 2, rank))
    xr <- sweep(xr, 2, sqrt(colSums(xr^2)), "/")
    mr <- sweep(mr, 2, sqrt(colSums(mr^2)), "/")
    rr <- vapply(seq_len(nrow(cis)), function(i) sum(xr[, cis$gene[i]] * mr[, cis$CpG[i]]), 0)
    rr <- pmax(-1, pmin(1, rr))
    df <- length(dev) - q$rank - 1
    cis$n <- length(dev)
    cis$df <- df
    cis$tissue_partial_rank_r <- rr
    cis$p_value <- 2 * pt(-abs(rr * sqrt(df/pmax(1e-15, 1 - rr^2))), df)
    cis$q_BH <- p.adjust(cis$p_value, "BH")
    wc(cis, "cis_associations_all.csv")
    casepairs <- read.csv(file.path(repo, "data/internal_response_pair_trace.csv"), stringsAsFactors = FALSE)
    casepairs <- casepairs[casepairs$cell_name %in% dev & is.finite(casepairs$response_value) & 
        casepairs$equal_ln_6dp_candidates > 0, ]
    cases <- data.frame(gene = c("MGMT", "SLFN11"), drug = c("Temozolomide", "Cisplatin"))
    cr <- list()
    cd <- list()
    availability <- list()
    for (i in seq_len(nrow(cases))) {
        g <- cases$gene[i]
        drug <- cases$drug[i]
        pr <- sort(unique(mapping$CpG[mapping$gene == g]))
        z <- casepairs[casepairs$drug_key == drug, ]
        availability[[i]] <- data.frame(gene = g, drug = drug, expression_present = g %in% colnames(x), 
            measured_annotated_probes = length(pr), observed_development_pairs = nrow(z))
        if (!g %in% colnames(x) || !nrow(z)) 
            next
        stopifnot(!anyDuplicated(z$cell_name))
        cells <- z$cell_name
        response <- z$response_value
        expr <- x[cells, g]
        tt <- unname(tissue[cells])
        cr[[paste(g, "expression")]] <- cbind(gene = g, drug = drug, CpG = NA_character_, association = "expression_response", 
            assoc(expr, response, tt))
        cd[[paste(g, "expression")]] <- data.frame(gene = g, drug = drug, CpG = NA_character_, cell = cells, 
            tissue = tt, expression = expr, methylation = NA_real_, response = response)
        for (cg in pr) {
            cr[[paste(g, cg, "response")]] <- cbind(gene = g, drug = drug, CpG = cg, association = "methylation_response", 
                assoc(m[cells, cg], response, tt))
            cr[[paste(g, cg, "expression")]] <- cbind(gene = g, drug = drug, CpG = cg, association = "methylation_expression", 
                assoc(m[dev, cg], x[dev, g], tissue[dev]))
            cd[[paste(g, cg)]] <- data.frame(gene = g, drug = drug, CpG = cg, cell = cells, tissue = tt, 
                expression = expr, methylation = m[cells, cg], response = response)
        }
    }
    cr <- do.call(rbind, cr)
    cr$q_BH_case_family <- p.adjust(cr$p_value, "BH")
    wc(cr, "case_study_associations.csv")
    wc(do.call(rbind, cd), "case_study_data.csv")
    wc(do.call(rbind, availability), "case_study_availability.csv")
    message("Methylation enrichment")
    genes <- sort(unique(mapping$gene))
    valid <- intersect(genes, keys(org.Hs.eg.db, keytype = "SYMBOL"))
    gm <- AnnotationDbi::select(org.Hs.eg.db, valid, "ENTREZID", "SYMBOL")
    gm <- gm[!is.na(gm$ENTREZID), ]
    counts <- table(gm$SYMBOL)
    gm <- gm[counts[gm$SYMBOL] == 1, ]
    mp <- merge(mapping, gm, by.x = "gene", by.y = "SYMBOL")
    mp <- unique(mp[, c("CpG", "ENTREZID")])
    ann <- AnnotationDbi::select(org.Hs.eg.db, unique(mp$ENTREZID), c("GOALL", "ONTOLOGYALL"), "ENTREZID")
    ann <- unique(ann[!is.na(ann$GOALL) & ann$ONTOLOGYALL == "BP", c("ENTREZID", "GOALL")])
    universe <- sort(unique(ann$ENTREZID))
    terms <- lapply(split(ann$ENTREZID, ann$GOALL), unique)
    terms <- terms[lengths(terms) >= 10 & lengths(terms) <= 500]
    tn <- AnnotationDbi::select(GO.db, names(terms), "TERM", "GOID")
    st <- read.csv(stability_path)
    np <- vapply(universe, function(e) length(unique(mp$CpG[mp$ENTREZID == e])), 1L)
    en <- list()
    qc <- list()
    for (gr in c("methylation_joint", "methylation_specific")) {
        probes <- st$feature[st$group == gr & st$stable_top5pct]
        query <- intersect(unique(mp$ENTREZID[mp$CpG %in% probes]), universe)
        qc[[gr]] <- data.frame(group = gr, selected_probes = length(probes), query_genes = length(query), 
            universe_genes = length(universe))
        sel <- as.integer(universe %in% query)
        n <- length(query)
        N <- length(universe)
        en[[gr]] <- do.call(rbind, lapply(names(terms), function(term) {
            mem <- as.integer(universe %in% terms[[term]])
            k <- sum(sel * mem)
            M <- sum(mem)
            pv <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
            fit <- suppressWarnings(glm(sel ~ mem + log1p(np), family = binomial()))
            cf <- summary(fit)$coefficients
            ok <- isTRUE(fit$converged) && "mem" %in% rownames(cf) && all(is.finite(cf["mem", ])) && 
                max(abs(coef(fit)), na.rm = TRUE) < 20
            data.frame(group = gr, GO_ID = term, term = tn$TERM[match(term, tn$GOID)], overlap = k, 
                query_size = n, term_size = M, universe_size = N, p_ORA = pv, probe_count_adjusted_logOR = if (ok) 
                  cf["mem", 1]
                else NA_real_, p_probe_count_adjusted = if (ok) 
                  cf["mem", 4]
                else NA_real_, logistic_valid = ok)
        }))
    }
    en <- do.call(rbind, en)
    en$q_ORA_BH <- p.adjust(en$p_ORA, "BH")
    en$q_probe_count_adjusted_BH <- p.adjust(en$p_probe_count_adjusted, "BH")
    wc(en, "methylation_GO_BP_all.csv")
    wc(do.call(rbind, qc), "methylation_enrichment_coverage.csv")
    wc(data.frame(entrez = universe, measured_probe_count = np), "methylation_enrichment_universe.csv")
    png(file.path(out, "case_studies.png"), width = 2200, height = 1300, res = 160)
    par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
    for (i in seq_len(nrow(cases))) {
        z <- cd[[paste(cases$gene[i], "expression")]]
        plot(z$expression, z$response, pch = 16, col = adjustcolor("#176B87", alpha.f = 0.45), xlab = paste(cases$gene[i], 
            "expression (published scale)"), ylab = "ln(IC50 / 1 micromolar)", main = paste(cases$gene[i], 
            "-", cases$drug[i]))
        abline(lm(response ~ expression, z), col = "#AD4130", lwd = 2)
    }
    dev.off()
    check <- assoc(m[dev, cis$CpG[1]], x[dev, cis$gene[1]], tissue[dev])
    stopifnot(abs(check$tissue_partial_rank_r - cis$tissue_partial_rank_r[1]) < 1e-10)
    ref <- read.csv(file.path(repo, "data/cis_associations_all.csv"))
    ix <- match(paste(cis$CpG, cis$gene), paste(ref$CpG, ref$gene))
    stopifnot(!anyNA(ix), max(abs(cis$tissue_partial_rank_r - ref$tissue_partial_rank_r[ix])) < 
        1e-10, max(abs(cis$q_BH - ref$q_BH[ix])) < 1e-10)
    ref <- read.csv(file.path(repo, "data/methylation_GO_all.csv"))
    ix <- match(paste(en$group, en$GO_ID), paste(ref$group, ref$GO_ID))
    stopifnot(!anyNA(ix), max(abs(en$q_ORA_BH - ref$q_ORA_BH[ix])) < 1e-10)
    stopifnot(identical(en$logistic_valid, ref$logistic_valid[ix]))
    for (n in c("probe_count_adjusted_logOR", "p_probe_count_adjusted", "q_probe_count_adjusted_BH")) {
        stopifnot(identical(is.na(en[[n]]), is.na(ref[[n]][ix])))
        k <- is.finite(en[[n]])
        stopifnot(max(abs(en[[n]][k] - ref[[n]][ix][k])) < 1e-08)
    }
    ref <- read.csv(file.path(repo, "data/case_study_associations.csv"))
    ix <- match(paste(cr$gene, cr$association), paste(ref$gene, ref$association))
    stopifnot(!anyNA(ix), max(abs(cr$tissue_partial_rank_r - ref$tissue_partial_rank_r[ix])) < 1e-10, 
        max(abs(cr$q_BH_case_family - ref$q_BH_case_family[ix])) < 1e-10)
    writeLines("PASS: cis associations, case studies and methylation ORA match recorded values", 
        file.path(out, "verification.txt"))
    capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
}

run_expression_enrichment <- function(repo, out, stability_path = file.path(repo, "data/feature_stability.csv")) {
    suppressPackageStartupMessages({
        library(org.Hs.eg.db)
        library(GO.db)
        library(AnnotationDbi)
    })
    stopifnot(as.character(packageVersion("org.Hs.eg.db")) == "3.22.0", as.character(packageVersion("GO.db")) == 
        "3.22.0")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    wc <- function(x, n) write.csv(x, file.path(out, n), row.names = FALSE, na = "")
    check <- function(n, b) {
        if (!isTRUE(b)) 
            stop(n)
    }
    f <- read.csv(stability_path)
    ex <- f[f$group == "expression_joint", ]
    ids <- ex$feature
    check("8046_unique_expression_features", length(ids) == 8046 && !anyDuplicated(ids))
    exactkeys <- intersect(ids, keys(org.Hs.eg.db, keytype = "SYMBOL"))
    direct <- AnnotationDbi::select(org.Hs.eg.db, keys = exactkeys, keytype = "SYMBOL", columns = "ENTREZID")
    missing <- setdiff(ids, exactkeys)
    ak <- intersect(missing, keys(org.Hs.eg.db, keytype = "ALIAS"))
    alias <- AnnotationDbi::select(org.Hs.eg.db, keys = ak, keytype = "ALIAS", columns = "ENTREZID")
    map <- do.call(rbind, lapply(ids, function(id) {
        v <- if (id %in% exactkeys) 
            direct$ENTREZID[direct$SYMBOL == id]
        else alias$ENTREZID[alias$ALIAS == id]
        v <- sort(unique(v[!is.na(v)]))
        data.frame(feature = id, status = if (length(v) == 0) 
            "unmapped"
        else if (length(v) > 1) 
            "ambiguous_excluded"
        else if (id %in% exactkeys) 
            "exact_symbol"
        else "unique_alias", candidate_entrez = paste(v, collapse = ";"), entrez = if (length(v) == 
            1) 
            v
        else NA_character_)
    }))
    entrez <- unique(na.omit(map$entrez))
    syms <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = c("SYMBOL", 
        "GENENAME"))
    check("one_current_symbol_per_entrez", !anyDuplicated(syms$ENTREZID))
    map$current_symbol <- syms$SYMBOL[match(map$entrez, syms$ENTREZID)]
    map$gene_name <- syms$GENENAME[match(map$entrez, syms$ENTREZID)]
    ann <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = c("GOALL", 
        "ONTOLOGYALL", "EVIDENCEALL"))
    ann <- ann[!is.na(ann$GOALL) & ann$ONTOLOGYALL == "BP", ]
    wc(ann, "measured_gene_BP_annotations.csv")
    uni <- sort(unique(ann$ENTREZID))
    map$in_BP_universe <- !is.na(map$entrez) & map$entrez %in% uni
    wc(map, "gene_identifier_mapping.csv")
    wc(data.frame(entrez = uni, symbol = syms$SYMBOL[match(uni, syms$ENTREZID)]), "background_genes.csv")
    termsets <- lapply(split(ann$ENTREZID, ann$GOALL), unique)
    termsets <- termsets[lengths(termsets) >= 10 & lengths(termsets) <= 500]
    terms <- AnnotationDbi::select(GO.db, keys = names(termsets), keytype = "GOID", columns = c("TERM", 
        "ONTOLOGY"))
    check("all_terms_BP", all(terms$ONTOLOGY == "BP") && !anyDuplicated(terms$GOID))
    wc(data.frame(GO_ID = names(termsets), background_count = lengths(termsets), description = terms$TERM[match(names(termsets), 
        terms$GOID)]), "eligible_terms.csv")
    N <- length(uni)
    allres <- list()
    coverage <- list()
    queries <- list()
    for (t in c(5, 1, 10)) {
        tr <- list()
        for (g in c("expression_joint", "expression_specific")) {
            rows <- f[f$group == g & f[[paste0("stable_top", t, "pct")]], ]
            mp <- map[match(rows$feature, map$feature), ]
            q <- sort(intersect(na.omit(mp$entrez), uni))
            n <- length(q)
            check(paste(g, t, "query_subset"), n > 0 && all(q %in% uni))
            queries[[paste(g, t)]] <- data.frame(group = g, threshold_pct = t, entrez = q, symbol = syms$SYMBOL[match(q, 
                syms$ENTREZID)])
            coverage[[paste(g, t)]] <- data.frame(group = g, threshold_pct = t, selected_features = nrow(rows), 
                mapped_features = sum(!is.na(mp$entrez)), ambiguous = sum(mp$status == "ambiguous_excluded"), 
                unmapped = sum(mp$status == "unmapped"), unique_mapped_genes = length(unique(na.omit(mp$entrez))), 
                BP_query_genes = n, BP_background_genes = N)
            overlaps <- lapply(termsets, intersect, y = q)
            k <- lengths(overlaps)
            M <- lengths(termsets)
            pv <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
            tr[[g]] <- data.frame(group = g, threshold_pct = t, GO_ID = names(termsets), description = terms$TERM[match(names(termsets), 
                terms$GOID)], overlap = k, query_size = n, term_background_size = M, background_size = N, 
                fold_enrichment = (k/n)/(M/N), p_value = pv, overlap_entrez = vapply(overlaps, paste, 
                  collapse = ";", FUN.VALUE = ""), overlap_symbols = vapply(overlaps, function(v) paste(syms$SYMBOL[match(v, 
                  syms$ENTREZID)], collapse = ";"), ""))
        }
        r <- do.call(rbind, tr)
        r$p_adjust_BH_family <- p.adjust(r$p_value, "BH")
        r$analysis <- if (t == 5) 
            "primary"
        else "sensitivity"
        r <- r[order(r$p_adjust_BH_family, r$p_value, r$GO_ID), ]
        allres[[as.character(t)]] <- r
    }
    res <- do.call(rbind, allres)
    wc(res, "GO_BP_all_tests.csv")
    primary <- res[res$threshold_pct == 5, ]
    wc(primary, "GO_BP_primary_all_tests.csv")
    wc(primary[primary$p_adjust_BH_family < 0.05, ], "GO_BP_primary_FDR05.csv")
    wc(do.call(rbind, coverage), "mapping_coverage.csv")
    wc(do.call(rbind, queries), "query_gene_lists.csv")
    z <- primary[order(primary$p_value), ]
    ix <- unique(c(seq_len(min(10, nrow(z))), nrow(z)))
    for (i in ix) {
        d <- z[i, ]
        tab <- matrix(c(d$overlap, d$query_size - d$overlap, d$term_background_size - d$overlap, 
            d$background_size - d$query_size - d$term_background_size + d$overlap), 2)
        fp <- fisher.test(tab, alternative = "greater")$p.value
        check(paste0("fisher_agreement_", i), abs(log(fp) - log(d$p_value)) < 1e-08)
    }
    check("zero_overlap_retained", all(res$p_value[res$overlap == 0] == 1))
    check("six_complete_test_sets", nrow(res) == 6 * length(termsets))
    ref <- read.csv(file.path(repo, "data/expression_GO_all.csv"))
    ix <- match(paste(ref$group, ref$threshold_pct, ref$GO_ID), paste(res$group, res$threshold_pct, 
        res$GO_ID))
    stopifnot(!anyNA(ix), max(abs(res$p_adjust_BH_family[ix] - ref$p_adjust_BH_family)) < 1e-10)
    writeLines("PASS: expression GO enrichment matches recorded values", file.path(out, "verification.txt"))
    capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
}

run_feature_stability <- function(repo, cv_dir, out) {
    files <- file.path(cv_dir, paste0("fold", 1:5, "_O2PLS_transform.rds"))
    if (!all(file.exists(files))) 
        stop("Run 03_train_and_tune_models.R first.")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    groups <- c(expression_joint = "W.", methylation_joint = "C.", expression_specific = "P_Yosc.", 
        methylation_specific = "P_Xosc.")
    models <- lapply(files, readRDS)
    results <- lapply(names(groups), function(group) {
        matrices <- lapply(models, function(m) m$fit[[groups[[group]]]])
        ids <- rownames(matrices[[1]])
        stopifnot(length(ids) > 0L, all(vapply(matrices, function(m) identical(rownames(m), ids), 
            logical(1))))
        leverage <- sapply(matrices, function(m) {
            q <- qr(m)
            stopifnot(q$rank == ncol(m))
            rowSums(qr.Q(q)^2)
        })
        ranks <- apply(leverage, 2, function(h) {
            r <- integer(length(h))
            r[order(-h, ids)] <- seq_along(h)
            r
        })
        result <- data.frame(group = group, feature = ids, mean_leverage = rowMeans(leverage), mean_rank = rowMeans(ranks), 
            sd_rank = apply(ranks, 1, sd), min_rank = apply(ranks, 1, min), max_rank = apply(ranks, 
                1, max))
        for (fold in 1:5) {
            result[[paste0("rank_fold", fold)]] <- ranks[, fold]
            result[[paste0("leverage_fold", fold)]] <- leverage[, fold]
        }
        for (percent in c(1, 5, 10)) {
            count <- rowSums(ranks <= ceiling(percent/100 * length(ids)))
            result[[paste0("count_top", percent, "pct")]] <- count
            result[[paste0("stable_top", percent, "pct")]] <- count >= 4L
        }
        result[order(result$mean_rank, result$feature), ]
    })
    result <- do.call(rbind, results)
    write.csv(result, file.path(out, "feature_stability.csv"), row.names = FALSE)
    original <- read.csv(file.path(repo, "data/feature_stability.csv"))
    index <- match(paste(result$group, result$feature), paste(original$group, original$feature))
    stopifnot(!anyNA(index))
    comparison <- data.frame(group = result$group, feature = result$feature, rerun_count_top5pct = result$count_top5pct, 
        recorded_count_top5pct = original$count_top5pct[index])
    write.csv(comparison, file.path(out, "selection_comparison.csv"), row.names = FALSE)
    invisible(result)
}

run_feature_stability(repo, result_dir("02_cv"), result_dir("08_features"))

run_biology(repo, result_dir("09_biology"),
  stability_path = file.path(result_dir("08_features"), "feature_stability.csv"))
run_expression_enrichment(repo, result_dir("09_expression_GO"),
  stability_path = file.path(result_dir("08_features"), "feature_stability.csv"))
