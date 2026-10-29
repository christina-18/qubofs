# dedup_resequenced_libraries.R
#
# Identify duplicate cells from resequenced 10x libraries.
# Within each cohort, donor, and compartment, sample pairs
# sharing >=50% of 16-bp cell barcodes are treated as
# resequenced libraries.
# For duplicated barcodes, retain the cell with the
# highest nCount_RNA (ties resolved by sample ID).
#
# Sourced by 01_clean_cohorts.R.
dedup_keep_cells <- function(md, cohort_col = "prj", donor_col = "pid", sid_col = "sid",
                             tissue_col = "tissue", thr = 0.5, log_file = "dedup_log.csv") {
  ## md = so@meta.data; returns a logical vector (TRUE = keep), named by cell
  md$.cell <- rownames(md)
  # 16-bp barcode core, wherever it sits in the cell name (prefix or suffix style)
  md$.core <- regmatches(md$.cell, regexpr("[ACGT]{16}", md$.cell))
  stopifnot(length(md$.core) == nrow(md))   # every cell name must contain a 16-bp barcode
  keep <- rep(TRUE, nrow(md)); names(keep) <- md$.cell
  log <- list()
  grp <- split(md, list(md[[cohort_col]], md[[donor_col]], md[[tissue_col]]), drop = TRUE)
  for (g in grp) {
    sids <- sort(unique(g[[sid_col]]))
    if (length(sids) < 2) next
    cores <- lapply(sids, function(s) unique(g$.core[g[[sid_col]] == s])); names(cores) <- sids
    parent <- setNames(sids, sids)
    find <- function(x) { while (parent[[x]] != x) x <- parent[[x]]; x }
    for (i in seq_along(sids)) for (j in seq_along(sids)) if (i < j) {
      a <- sids[i]; b <- sids[j]
      f <- length(intersect(cores[[a]], cores[[b]])) / min(length(cores[[a]]), length(cores[[b]]))
      log[[length(log) + 1]] <- data.frame(cohort = g[[cohort_col]][1], donor = g[[donor_col]][1],
                                           sid_a = a, sid_b = b, shared_frac = round(f, 4))
      if (f >= thr) parent[[find(b)]] <- find(a)
    }
    lib <- vapply(g[[sid_col]], find, character(1))
    for (L in unique(lib)) {
      h <- g[lib == L, ]
      if (length(unique(h[[sid_col]])) < 2) next
      h <- h[order(h$.core, -h$nCount_RNA, h[[sid_col]]), ]
      keep[h$.cell[duplicated(h$.core)]] <- FALSE
    }
  }
  if (length(log)) write.csv(do.call(rbind, log), log_file, row.names = FALSE)
  message(sprintf("dedup: %d duplicate cells flagged of %d", sum(!keep), length(keep)))
  keep
}

## convenience wrapper (makes a subset copy; prefer dedup_keep_cells for big objects)
dedup_resequenced <- function(so, ...) {
  keep <- dedup_keep_cells(so@meta.data, ...)
  subset(so, cells = names(keep)[keep])
}
