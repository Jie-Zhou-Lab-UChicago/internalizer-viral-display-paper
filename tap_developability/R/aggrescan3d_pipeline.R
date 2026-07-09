## =============================================================================
## Aggrescan3D (A3D) COMPUTE pipeline
## Reads per-VHH A3D output (per-residue `score`) + the A3D structure and writes
## a master table of total / mean / per-CDR / CDR-vicinity aggregation scores.
## Group comparison & figures are in analysis/aggregation_phage_vs_viral.R
## =============================================================================
suppressPackageStartupMessages({ library(bio3d); library(stringr) })

## ------------------------------- CONFIG --------------------------------------
A3D_DIR     <- "data/VHH_AD3"                  # Aggrescan3D results; one sub-folder per ID (Ab_viral_1, ...)
CSV_PATH    <- "data/VHH_to_analyze.csv"       # columns: <index>, ID, seq
CHAIN       <- "A"
VIC_CUTOFF  <- 4.0                             # Å, CDR-vicinity neighbour cutoff
PDB_PATTERN <- "input\\.pdb$"                  # the A3D structure (log showed folded/input.pdb) - VERIFY
A3D_CSV_PAT <- "\\.csv$"                        # the A3D per-residue csv (e.g. A3D.csv)
OUT_TABLE   <- "Aggregation3D_analysis.csv"
OUT_DIR     <- "results"

## ---- helpers: score keyed by residue number --------------------------------
score_by_resno <- function(csv, chain = NULL) {
  nm <- tolower(names(csv))
  rescol   <- names(csv)[which(nm %in% c("residue","resid","resno","res_num","position"))[1]]
  scorecol <- names(csv)[which(nm %in% c("score","a3d","a3d_score"))[1]]
  if (is.na(scorecol)) stop("no score column found in A3D csv")
  if (!is.null(chain) && "chain" %in% nm) {
    cc <- names(csv)[which(nm == "chain")[1]]; csv <- csv[csv[[cc]] == chain, ]
  }
  if (is.na(rescol)) { warning("A3D csv has no residue-number column; using row position")
    return(setNames(csv[[scorecol]], seq_len(nrow(csv)))) }
  setNames(csv[[scorecol]], as.integer(csv[[rescol]]))
}
sum_scores <- function(sb, res_vec) sum(sb[as.character(res_vec)], na.rm = TRUE)

## ======================= BUILD TABLE + CDRs ==================================
c <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
c <- c[, -1]                                    # drop leading index -> ID, seq

folders  <- list.files(A3D_DIR, pattern = "^Ab_")
c <- c[c$ID %in% folders, ]
c <- c[match(folders[folders %in% c$ID], c$ID), ]
file.name <- file.path(A3D_DIR, c$ID)
stopifnot(identical(basename(file.name), c$ID)) # aborts if ever misaligned

c$cdr1 <- substr(c$seq, 26, 35)
c$cdr2 <- substr(c$seq, 49, 59)
c$cdr3 <- vapply(c$seq, function(s) {
  m <- str_locate(s, "WGQGTLVTVSS")[1]
  if (is.na(m) || (m-1) < 97) NA_character_ else substr(s, 97, m-1)
}, character(1))
c$cdr_length <- nchar(c$cdr1) + nchar(c$cdr2) + ifelse(is.na(c$cdr3), NA, nchar(c$cdr3))
c$cdr1_range <- "26-35"; c$cdr2_range <- "49-59"; c$cdr3_range <- paste0("97-", nchar(c$seq)-11)
for (col in c("total_agg_score","mean_agg_score","cdr1_agg_score","cdr2_agg_score",
              "cdr3_agg_score","vicinity_res","vicinity_agg_score","error")) c[[col]] <- NA

## ======================= PER-STRUCTURE LOOP ==================================
for (i in seq_len(length(file.name))) {
  c$error[i] <- tryCatch({
    csv.file <- list.files(file.name[i], pattern = A3D_CSV_PAT, full.names = TRUE)[1]
    if (is.na(csv.file)) stop("no A3D csv")
    a3d <- read.csv(csv.file, stringsAsFactors = FALSE)
    sb  <- score_by_resno(a3d, chain = CHAIN)

    cdr1 <- 26:35; cdr2 <- 49:59; cdr3 <- 97:(nchar(c$seq[i]) - 11)
    cdr_residues <- sort(unique(c(cdr1, cdr2, cdr3)))

    c$total_agg_score[i] <- sum(sb, na.rm = TRUE)
    c$mean_agg_score[i]  <- mean(sb, na.rm = TRUE)
    c$cdr1_agg_score[i]  <- sum_scores(sb, cdr1)
    c$cdr2_agg_score[i]  <- sum_scores(sb, cdr2)
    c$cdr3_agg_score[i]  <- sum_scores(sb, cdr3)

    pdb.file <- list.files(file.name[i], pattern = PDB_PATTERN, full.names = TRUE)[1]
    if (is.na(pdb.file)) stop("no A3D pdb")
    pdb <- read.pdb(pdb.file)
    cdr_sel <- atom.select(pdb, "protein", chain = CHAIN, resno = cdr_residues, verbose = FALSE)
    if (length(cdr_sel$atom) == 0L) { c$vicinity_res[i] <- ""; return("no_cdr_atoms") }
    d <- dm(pdb$xyz, mask.lower = FALSE)
    within_cols <- which(apply(d[cdr_sel$atom, , drop = FALSE], 2,
                               function(x) any(x > 0 & x < VIC_CUTOFF, na.rm = TRUE)))
    vic_res <- sort(unique(as.integer(pdb$atom$resno[within_cols])))
    c$vicinity_res[i]       <- paste(vic_res, collapse = ",")
    c$vicinity_agg_score[i] <- sum_scores(sb, vic_res)
    NA_character_
  }, error = function(e) as.character(conditionMessage(e)))
  if (i %% 100 == 0) message(sprintf("... %d / %d", i, length(file.name)))
}

dir.create(OUT_DIR, showWarnings = FALSE)
write.csv(c, file.path(OUT_DIR, OUT_TABLE), row.names = FALSE)
message(sprintf("Done. %d VHHs -> %s (%d errors).",
                nrow(c), file.path(OUT_DIR, OUT_TABLE), sum(!is.na(c$error))))
