#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(GO.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: classify_tf_go.R INPUT.cys OUTPUT.csv")
}

input_cys <- normalizePath(args[[1]], mustWork = TRUE)
output_csv <- args[[2]]

entries <- unzip(input_cys, list = TRUE)$Name
node_table <- entries[grepl(
  "tables/87-M\\.\\+musculus\\+%281%29/SHARED_ATTRS-org\\.cytoscape\\.model\\.CyNode",
  entries
)][[1]]

raw <- read.csv(
  unz(input_cys, node_table),
  skip = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  colClasses = "character"
)
nodes <- raw[-seq_len(3), , drop = FALSE]
symbols <- sort(unique(nodes$`gene name`[nodes$`node type` == "query"]))

tf_root <- "GO:0003700"
tf_go_ids <- unique(c(tf_root, GOMFOFFSPRING[[tf_root]]))
annotations <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = symbols,
  keytype = "SYMBOL",
  columns = c("ENTREZID", "GO", "ONTOLOGY", "EVIDENCE")
)

tf_hits <- annotations[
  annotations$ONTOLOGY == "MF" &
    annotations$GO %in% tf_go_ids &
    !is.na(annotations$GO),
  ,
  drop = FALSE
]
tf_hits$TERM <- Term(GOTERM[tf_hits$GO])

collapse_values <- function(values) {
  values <- sort(unique(values[!is.na(values) & nzchar(values)]))
  paste(values, collapse = "; ")
}

results <- lapply(symbols, function(symbol) {
  hit <- tf_hits[tf_hits$SYMBOL == symbol, , drop = FALSE]
  data.frame(
    gene_name = symbol,
    node_type = if (nrow(hit) > 0) "TF" else "DEG",
    entrez_id = collapse_values(hit$ENTREZID),
    tf_database_type = "Gene Ontology - Molecular Function",
    tf_database_version = paste0(
      "org.Mm.eg.db ", packageVersion("org.Mm.eg.db"),
      "; GO.db ", packageVersion("GO.db")
    ),
    tf_reference_rule = paste0(
      tf_root,
      " DNA-binding transcription factor activity and descendants"
    ),
    tf_go_terms = collapse_values(hit$GO),
    tf_go_term_names = collapse_values(hit$TERM),
    tf_evidence_codes = collapse_values(hit$EVIDENCE),
    classification_date = "2026-08-06",
    stringsAsFactors = FALSE
  )
})

result <- do.call(rbind, results)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(result, output_csv, row.names = FALSE, fileEncoding = "UTF-8")

tf_symbols <- result$gene_name[result$node_type == "TF"]
cat("Query genes:", nrow(result), "\n")
cat("Strict GO-MF TFs:", length(tf_symbols), "\n")
cat("TF symbols:", paste(tf_symbols, collapse = ", "), "\n")
