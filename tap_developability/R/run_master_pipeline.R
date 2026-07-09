## =============================================================================
## VHH developability master table: CDRs, vicinity, surface exposure, PSH, PPC/PNC
## Standalone & runnable. Requires: bio3d, stringr, and the FreeSASA .rsa files.
## =============================================================================
suppressPackageStartupMessages({
  library(bio3d)
  library(stringr)
})

## ------------------------------- CONFIG --------------------------------------
CSV_PATH   <- "data/VHH_to_analyze.csv"     # columns: <index>, ID, seq
STRUCT_DIR <- "data/VHH_model_str"          # predicted VHH structures + FreeSASA .rsa; one sub-folder per ID
CHAIN      <- "A"                            # VHH chain id in the PDBs
VIC_CUTOFF <- 4.0                           # Å, CDR-vicinity neighbour cutoff
EXPO_CUT   <- 15                            # % relative side-chain SASA -> "surface exposed"
SCORE_SET  <- "cdr"                         # "cdr" = score exposed CDR residues (original behaviour)
                                            # "vicinity" = score exposed CDR + framework neighbours (more TAP-like)
OUT_CSV    <- "results/analysis.csv"

## ======================= 1. HELPER FUNCTIONS =================================

## ---- PSH (patches of surface hydrophobicity) --------------------------------
psh_from_residue_string <- function(pdb, residue_str, chain = "A",
                                    scale = 3,            # 1=Wimley-White, 2=Eisenberg-McLachlan, 3=Kyte-Doolittle (TAP)
                                    cutoff = 7.5, salt_cut = 3.2,
                                    normalize_to_1_2 = TRUE) {
  ww_interfacial <- c(A=0.17,R=2.58,N=0.42,D=1.23,C=-0.24,Q=0.58,E=2.02,G=0.74,H=0.11,
                      I=-0.31,L=-0.56,K=2.71,M=-0.23,F=-1.13,P=0.45,S=0.13,T=0.25,W=-1.85,Y=-0.94,V=-0.07)
  em_consensus   <- c(A=0.62,R=-2.53,N=-0.78,D=-0.90,C=0.29,Q=-0.85,E=-0.74,G=0.48,H=-0.40,
                      I=1.38,L=1.06,K=-1.50,M=0.64,F=1.19,P=0.12,S=-0.18,T=-0.05,W=0.81,Y=0.26,V=1.08)
  kyte_doolittle <- c(I=4.5,V=4.2,L=3.8,F=2.8,C=2.5,M=1.9,A=1.8,G=-0.4,T=-0.7,S=-0.8,
                      W=-0.9,Y=-1.3,P=-1.6,H=-3.2,E=-3.5,Q=-3.5,D=-3.5,N=-3.5,K=-3.9,R=-4.5)
  scale_raw <- switch(as.character(scale),
                      "1" = -ww_interfacial, "2" = em_consensus, "3" = kyte_doolittle,
                      stop("scale must be 1, 2, or 3"))
  if (normalize_to_1_2) {
    rng <- range(scale_raw, na.rm = TRUE); H12 <- 1 + (scale_raw - rng[1])/(rng[2]-rng[1])
  } else H12 <- scale_raw
  aa3to1 <- c(ALA="A",ARG="R",ASN="N",ASP="D",CYS="C",GLN="Q",GLU="E",GLY="G",HIS="H",ILE="I",
              LEU="L",LYS="K",MET="M",PHE="F",PRO="P",SER="S",THR="T",TRP="W",TYR="Y",VAL="V")
  get_H <- function(r3, gly=FALSE) { a<-aa3to1[r3]; if(is.na(a)) return(NA_real_); if(gly) H12["G"] else H12[a] }
  resnos <- unique(na.omit(as.integer(strsplit(gsub("\\s+","",residue_str), ",")[[1]])))
  if (length(resnos) < 2) return(0.0)
  at_all <- pdb$atom; heavy <- !grepl("^H", at_all$elety)
  at_all <- at_all[heavy,,drop=FALSE]; coords <- matrix(pdb$xyz,ncol=3,byrow=TRUE)[heavy,,drop=FALSE]
  keep <- at_all$chain==chain & at_all$resno %in% resnos
  at <- at_all[keep,,drop=FALSE]; if(!nrow(at)) return(0.0); coords <- coords[keep,,drop=FALSE]
  key <- paste(at$chain, at$resno, sep=":")
  resname_lookup <- tapply(at$resid, key, function(x) x[1])
  idx_by_res <- split(seq_len(nrow(at)), key); keys <- names(idx_by_res)
  if (length(keys) < 2) return(0.0)
  is_posN <- (at$resid=="LYS"&at$elety=="NZ")|(at$resid=="ARG"&at$elety%in%c("NE","NH1","NH2"))
  is_negO <- (at$resid=="ASP"&at$elety%in%c("OD1","OD2"))|(at$resid=="GLU"&at$elety%in%c("OE1","OE2"))
  pos_idx<-which(is_posN); neg_idx<-which(is_negO); salt_keys<-character(0)
  if (length(pos_idx)&&length(neg_idx)) for (ii in pos_idx) {
    d2 <- colSums((t(coords[neg_idx,,drop=FALSE]) - coords[ii,])^2)
    if (length(d2) && sqrt(min(d2)) <= salt_cut) {
      jj <- neg_idx[which.min(d2)]
      salt_keys <- unique(c(salt_keys, paste(at$chain[ii],at$resno[ii],sep=":"),
                            paste(at$chain[jj],at$resno[jj],sep=":")))
    }
  }
  Hmap <- setNames(mapply(function(r3,k) get_H(r3, k %in% salt_keys),
                          unname(resname_lookup[keys]), keys), keys)
  PSH <- 0.0; n <- length(keys)
  for (i in 1:(n-1)) { ci<-coords[idx_by_res[[keys[i]]],,drop=FALSE]
    for (j in (i+1):n) { cj<-coords[idx_by_res[[keys[j]]],,drop=FALSE]
      d2<-as.matrix(dist(rbind(ci,cj)))^2; nci<-nrow(ci)
      r2min<-min(d2[seq_len(nci), nci+seq_len(nrow(cj)), drop=FALSE])
      if (is.finite(r2min) && sqrt(r2min) < cutoff) PSH <- PSH + (Hmap[[keys[i]]]*Hmap[[keys[j]]])/r2min
    } }
  PSH
}

## ---- PPC / PNC (patches of positive / negative charge) ----------------------
ppc_pnc_from_residue_string <- function(pdb, residue_str, chain = "A",
                                        cutoff = 7.5, salt_cut = 3.2, his_charge = 0.1) {
  empty <- list(PPC=0.0, PNC=0.0)
  resnos <- unique(na.omit(as.integer(strsplit(gsub("\\s+","",residue_str), ",")[[1]])))
  if (length(resnos) < 2) return(empty)
  at_all <- pdb$atom; heavy <- !grepl("^H", at_all$elety)
  at_all <- at_all[heavy,,drop=FALSE]; coords <- matrix(pdb$xyz,ncol=3,byrow=TRUE)[heavy,,drop=FALSE]
  keep <- at_all$chain==chain & at_all$resno %in% resnos
  at <- at_all[keep,,drop=FALSE]; if(!nrow(at)) return(empty); coords <- coords[keep,,drop=FALSE]
  key <- paste(at$chain, at$resno, sep=":")
  resname_lookup <- tapply(at$resid, key, function(x) x[1])
  idx_by_res <- split(seq_len(nrow(at)), key); keys <- names(idx_by_res)
  if (length(keys) < 2) return(empty)
  is_posN <- (at$resid=="LYS"&at$elety=="NZ")|(at$resid=="ARG"&at$elety%in%c("NE","NH1","NH2"))
  is_negO <- (at$resid=="ASP"&at$elety%in%c("OD1","OD2"))|(at$resid=="GLU"&at$elety%in%c("OE1","OE2"))
  pos_idx<-which(is_posN); neg_idx<-which(is_negO); salt_keys<-character(0)
  if (length(pos_idx)&&length(neg_idx)) for (ii in pos_idx) {
    d2 <- colSums((t(coords[neg_idx,,drop=FALSE]) - coords[ii,])^2)
    if (length(d2) && sqrt(min(d2)) <= salt_cut) {
      jj <- neg_idx[which.min(d2)]
      salt_keys <- unique(c(salt_keys, paste(at$chain[ii],at$resno[ii],sep=":"),
                            paste(at$chain[jj],at$resno[jj],sep=":")))
    }
  }
  rn <- unname(resname_lookup[keys])
  q  <- ifelse(rn %in% c("LYS","ARG"), 1, ifelse(rn=="HIS", his_charge, ifelse(rn %in% c("ASP","GLU"), -1, 0)))
  q[keys %in% salt_keys] <- 0
  Qmap <- setNames(abs(q), keys); signv <- setNames(ifelse(q>0,"pos",ifelse(q<0,"neg","zero")), keys)
  sum_patch <- function(which_sign) {
    kv <- keys[signv[keys]==which_sign]; kv <- kv[Qmap[kv] > 0]
    if (length(kv) < 2) return(0.0)
    total<-0.0; n<-length(kv)
    for (i in 1:(n-1)) { ci<-coords[idx_by_res[[kv[i]]],,drop=FALSE]; Qi<-Qmap[[kv[i]]]
      for (j in (i+1):n) { cj<-coords[idx_by_res[[kv[j]]],,drop=FALSE]; Qj<-Qmap[[kv[j]]]
        d2<-as.matrix(dist(rbind(ci,cj)))^2; nci<-nrow(ci)
        r2min<-min(d2[seq_len(nci), nci+seq_len(nrow(cj)), drop=FALSE])
        if (is.finite(r2min) && sqrt(r2min) < cutoff) total <- total + (Qi*Qj)/r2min
      } }
    total
  }
  list(PPC = sum_patch("pos"), PNC = sum_patch("neg"))
}

## ---- FreeSASA .rsa reader ---------------------------------------------------
read_freesasa_rsa <- function(rsa_file) {
  stopifnot(file.exists(rsa_file))
  lines <- readLines(rsa_file, warn = FALSE)
  res_lines <- grep("^RES", lines, value = TRUE)
  res_lines <- gsub("\\bN/A\\b", "NA", res_lines)
  df <- read.table(text = res_lines, stringsAsFactors = FALSE)
  colnames(df) <- c("tag","resname","chain","resno","all_abs","all_rel","side_abs","side_rel",
                    "main_abs","main_rel","nonpol_abs","nonpol_rel","polar_abs","polar_rel")
  df$resno <- suppressWarnings(as.integer(df$resno))
  num_cols <- colnames(df)[5:14]; df[num_cols] <- lapply(df[num_cols], as.numeric)
  df[, -1]
}

## ======================= 2. BUILD TABLE + CDRs ===============================
df <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
df <- df[, -1]                                     # drop leading index column -> ID, seq

## CDR extraction (fixed scaffold; positions are SEQUENCE positions)
df$cdr1 <- substr(df$seq, 26, 35)
df$cdr2 <- substr(df$seq, 49, 59)
df$cdr3 <- vapply(df$seq, function(s) {
  m <- str_locate(s, "WGQGTLVTVSS")[1]
  if (is.na(m) || (m-1) < 97) NA_character_ else substr(s, 97, m-1)
}, character(1))
df$cdr_length  <- nchar(df$cdr1) + nchar(df$cdr2) + ifelse(is.na(df$cdr3), NA, nchar(df$cdr3))
df$cdr1_range  <- "26-35"
df$cdr2_range  <- "49-59"
df$cdr3_range  <- paste0("97-", nchar(df$seq) - 11)

## result columns
for (col in c("vicinity_res","vicinity_res_noncdr","vicinity_res_cdr","sur_expo_15","ser_cdr_15",
              "score_set","psh_15_ww","psh_15_em","psh_15_kd","ppc","pnc","error"))
  df[[col]] <- NA_character_

## ---- align structure folders to df BY ID (bulletproof) ----------------------
have <- dir.exists(file.path(STRUCT_DIR, df$ID))
if (!all(have)) message(sprintf("Note: %d/%d VHHs have no structure folder; skipping.", sum(!have), length(have)))
df <- df[have, ]
file.name <- file.path(STRUCT_DIR, df$ID)
stopifnot(identical(basename(file.name), df$ID))

## ======================= 3. PER-STRUCTURE LOOP ===============================
for (i in seq_len(nrow(df))) {
  df$error[i] <- tryCatch({

    pdb.files <- list.files(file.name[i], pattern = "\\.pdb$", full.names = TRUE)
    if (!length(pdb.files)) stop("no .pdb")
    pdb <- read.pdb(pdb.files[1])

    ## one-time numbering sanity check (resno must equal sequence position)
    if (i == 1) {
      ca <- pdb$atom[pdb$atom$elety=="CA" & pdb$atom$chain==CHAIN, ]
      if (max(ca$resno) != nchar(df$seq[i]))
        warning(sprintf("resno max (%d) != seq length (%d) - check PDB numbering!",
                        max(ca$resno), nchar(df$seq[i])))
    }

    ## CDR residues (sequence-position based)
    cdr1 <- 26:35; cdr2 <- 49:59
    cdr3 <- 97:(nchar(df$seq[i]) - 11)
    cdr_residues <- sort(unique(c(cdr1, cdr2, cdr3)))

    cdr_sel <- atom.select(pdb, "protein", chain = CHAIN, resno = cdr_residues, verbose = FALSE)
    if (length(cdr_sel$atom) == 0L) {
      df$vicinity_res[i] <- df$vicinity_res_noncdr[i] <- df$vicinity_res_cdr[i] <- ""
      return("no_cdr_atoms")
    }

    ## CDR vicinity: any atom within VIC_CUTOFF of a CDR atom (full symmetric matrix!)
    d <- dm(pdb$xyz, mask.lower = FALSE)
    within_cols <- which(apply(d[cdr_sel$atom, , drop = FALSE], 2,
                               function(x) any(x > 0 & x < VIC_CUTOFF, na.rm = TRUE)))
    vic_keys <- unique(paste(pdb$atom$chain[within_cols], pdb$atom$resno[within_cols], sep=":"))
    cdr_keys <- unique(paste(pdb$atom$chain[cdr_sel$atom], pdb$atom$resno[cdr_sel$atom], sep=":"))
    nbr_keys <- setdiff(vic_keys, cdr_keys)

    strip <- function(x) gsub(paste0(CHAIN, ":"), "", x)
    df$vicinity_res[i]        <- strip(paste(vic_keys, collapse=","))
    df$vicinity_res_noncdr[i] <- strip(paste(nbr_keys, collapse=","))
    df$vicinity_res_cdr[i]    <- strip(paste(cdr_keys, collapse=","))

    ## surface-exposed residues (FreeSASA rel side-chain SASA > EXPO_CUT %)
    rsa.files <- list.files(file.name[i], pattern = "\\.rsa$", full.names = TRUE)
    if (!length(rsa.files)) stop("no .rsa")
    rsa <- read_freesasa_rsa(rsa.files[1])
    exp_res <- sort(unique(na.omit(rsa$resno[rsa$chain == CHAIN & !is.na(rsa$side_rel) & rsa$side_rel > EXPO_CUT])))
    df$sur_expo_15[i] <- paste(exp_res, collapse=",")

    ## residue set to score: exposed  ∩  (CDR only | full vicinity)
    base_set <- if (SCORE_SET == "vicinity") strip(paste(vic_keys, collapse=",")) else df$vicinity_res_cdr[i]
    base_num <- as.integer(unlist(strsplit(base_set, ",")))
    score_res <- sort(intersect(base_num, exp_res))
    score_str <- paste(score_res, collapse=",")
    df$ser_cdr_15[i] <- score_str
    df$score_set[i]  <- SCORE_SET

    ## PSH (3 scales) + PPC/PNC on the scored set
    df$psh_15_ww[i] <- psh_from_residue_string(pdb, score_str, chain=CHAIN, scale=1)
    df$psh_15_em[i] <- psh_from_residue_string(pdb, score_str, chain=CHAIN, scale=2)
    df$psh_15_kd[i] <- psh_from_residue_string(pdb, score_str, chain=CHAIN, scale=3)  # TAP scale
    patch <- ppc_pnc_from_residue_string(pdb, score_str, chain=CHAIN,
                                         cutoff=7.5, salt_cut=3.2, his_charge=0.1)
    df$ppc[i] <- patch$PPC
    df$pnc[i] <- patch$PNC

    NA_character_                                   # no error
  },
  error = function(e) as.character(conditionMessage(e)))

  if (i %% 100 == 0) message(sprintf("... %d / %d", i, nrow(df)))
}

## ======================= 4. WRITE ============================================
## numeric columns come back as character via the NA init; coerce the scores:
for (col in c("cdr_length","psh_15_ww","psh_15_em","psh_15_kd","ppc","pnc"))
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))

write.csv(df, OUT_CSV, row.names = FALSE)
message(sprintf("Done. %d rows written to %s (%d errors).",
                nrow(df), OUT_CSV, sum(!is.na(df$error))))
