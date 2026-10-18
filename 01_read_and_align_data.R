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
required <- c("digest", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
expected <- jsonlite::fromJSON("config/raw_input_hashes.json")
for (p in names(expected)) {
    stopifnot(file.exists(p), digest::digest(file = p, algo = "sha256") == expected[[p]])
}
message("PASS: four original files and source metadata verified.")
read_table <- function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
objects <- list()
for (ds in c("expression", "methylation")) {
    filename <- if (ds == "expression") "Gene_expression.csv" else "Methylation.csv"
    raw <- read_table(file.path("raw_data", filename))
    map <- read_table(file.path("data/source_metadata", paste0(ds, "_feature_name_map.csv")))
    stopifnot(identical(map$current_header, names(raw)[-1L]), identical(map$numeric_column, seq_len(ncol(raw) - 1L)))
    mat <- as.matrix(raw[, -1L, drop = FALSE])
    rownames(mat) <- as.character(raw[[1L]])
    colnames(mat) <- map$source_feature_id
    stopifnot(is.numeric(mat), !anyDuplicated(rownames(mat)))
    objects[[ds]] <- mat
}
expected <- jsonlite::fromJSON("config/reconstructed_object_hashes.json")
checks <- data.frame(object = names(objects), passed = vapply(names(objects), function(n) digest::digest(objects[[n]], algo = "sha256") == expected[[n]], logical(1)))
stopifnot(all(checks$passed))
dir.create("results/01_data", recursive = TRUE, showWarnings = FALSE)
saveRDS(objects, "results/01_data/omics_from_csv.rds")
write.csv(checks, "results/01_data/omics_checks.csv", row.names = FALSE)
print(checks)
