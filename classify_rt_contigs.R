#!/usr/bin/env Rscript

# Rule-based classification of reverse-transcriptase (RT)-containing contigs from GFF3
# Usage: Rscript classify_rt_contigs.R seq_genomes.gff3 [output_prefix]
# Outputs:
#   <prefix>_rt_contig_classification.tsv
#   <prefix>_rt_feature_evidence.tsv
#   <prefix>_rt_class_summary.tsv
#
#
#
# Design principles:
# 1) Require RT evidence to enter the classification set.
# 2) Use element-specific annotations (CDD/Pfam/InterPro/PROSITE names) as strong evidence.
# 3) Use combinations of domains as supporting evidence.
# 4) Keep ambiguous cases explicit instead of forcing a class.




args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript classify_rt_contigs.R <input.gff3> [output_prefix]")
infile <- args[1]
prefix <- if (length(args) >= 2) args[2] else sub("\\.gff3?$", "", basename(infile), ignore.case = TRUE)


############################################################################################
################# helpers 
############################################################################################

url_decode <- function(x) utils::URLdecode(x)
extract_attr <- function(attr, key = "Name") {
  m <- regexec(paste0("(?:^|;)", key, "=([^;]*)"), attr, perl = TRUE)
  z <- regmatches(attr, m)
  out <- vapply(z, function(v) if (length(v) >= 2) v[2] else NA_character_, character(1))
  url_decode(out)
}
rx <- function(x, pattern) grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
any_rx <- function(x, pattern) any(rx(x, pattern), na.rm = TRUE)
collapse_unique <- function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse = " | ")


############################################################################################
################# read GFF3 
############################################################################################

gff <- read.delim(infile, comment.char = "#", header = FALSE, sep = "\t",
                  quote = "", stringsAsFactors = FALSE, fill = TRUE,
                  col.names = c("seqid","source","type","start","end","score","strand","phase","attributes"))
if (ncol(gff) < 9) stop("Input does not look like a valid 9-column GFF3 file")
gff$name <- extract_attr(gff$attributes, "Name")
gff$name[is.na(gff$name)] <- ""
gff$feature_text <- paste(gff$source, gff$type, gff$name)



############################################################################################
################## evidence patterns
############################################################################################

# Generic RT evidence
p_rt_generic <- paste(
  "reverse transcript", "RNA-directed DNA polymerase", "RNA dependent DNA polymerase",
  "RVT[_ -]?1", "RT_POL", "RT_RNaseH", "RT_G2_intron", "RT_nLTR", "RT_LTR",
  "endonuclease-reverse transcriptase", "endonuclease/reverse transcriptase",
  sep = "|"
)

# Strong subtype-specific evidence
p_nonltr_strong <- paste(
  "RT_nLTR", "Rnase_HI_RT_non_LTR", "non[- ]?LTR",
  "mobile element jockey", "LINE[- ]?like", "LINE retrotranspos", "R1 retrotranspos",
  sep = "|"
)
p_ltr_strong <- paste(
  "RT_LTR", "RNase_HI_RT_Ty3", "RNase_HI_RT_DIRS", "Ty3", "gypsy", "copia",
  "LTR retrotranspos", "retrotransposon.*gag-pol", "gag-pol polyprotein",
  sep = "|"
)
p_retrovirus_strong <- paste(
  "endogenous retrovirus", "retrovirus-related", "retroviral", "HIV", "simian retrovirus",
  "GAG RETROVIRAL", "gag protein", "Gag polyprotein", "viral envelope protein", "ENV polyprotein",
  sep = "|"
)
p_groupII_strong <- paste(
  "RT_G2_intron", "group II intron", "group 2 intron", "intron maturase",
  "reverse transcriptase/maturase", "NUCLEAR INTRON MATURASE",
  sep = "|"
)
p_groupII_support <- "GsI-IIC RT/DNA"
p_telomerase <- "telomerase|TERT|telomere recombination"
p_retron <- "retron|msDNA|multicopy single-stranded DNA"
p_dgr <- "diversity-generating retroelement|DGR"
p_crispr_rt <- "CRISPR.*reverse transcript|Cas6-RT-Cas1|RT-Cas1|CRISPR-associated.*RT"

# Supporting domains
p_integrase <- "integrase|\\brve\\b"
p_rnaseh <- "RNase[_ -]?H|ribonuclease H"
p_retro_protease <- "retroviral.*protease|retroviral-type|retropepsin|acid protease"
p_gag_capsid <- "\\bgag\\b|capsid.*retro|retrovirus capsid|matrix protein"
p_env <- "\\benv\\b|envelope protein|ENV polyprotein"
p_endo <- "endonuclease|Exo_endo_phos|DNase I-like"
p_zf_retro <- "retrovirus zinc finger|gag.*zinc|integrase.*zinc"


############################################################################################
################# classify each RT-containing contig 
############################################################################################
seqs <- unique(gff$seqid)
rows <- vector("list", length(seqs))
evidence_rows <- list(); ei <- 1L

for (i in seq_along(seqs)) {
  s <- seqs[i]
  d <- gff[gff$seqid == s, , drop = FALSE]
  txt <- d$feature_text
  nm <- d$name

  has_rt <- any_rx(txt, p_rt_generic)
  if (!has_rt) next

  # Evidence flags
  f_nonltr <- any_rx(txt, p_nonltr_strong)
  f_ltr <- any_rx(txt, p_ltr_strong)
  f_retrovirus <- any_rx(txt, p_retrovirus_strong)
  f_groupII <- any_rx(txt, p_groupII_strong)
  f_groupII_support <- any_rx(txt, p_groupII_support)
  f_telomerase <- any_rx(txt, p_telomerase)
  f_retron <- any_rx(txt, p_retron)
  f_dgr <- any_rx(txt, p_dgr)
  f_crispr <- any_rx(txt, p_crispr_rt)
  f_integrase <- any_rx(txt, p_integrase)
  f_rnaseh <- any_rx(txt, p_rnaseh)
  f_protease <- any_rx(txt, p_retro_protease)
  f_gag <- any_rx(txt, p_gag_capsid)
  f_env <- any_rx(txt, p_env)
  f_endo <- any_rx(txt, p_endo)
  f_zf <- any_rx(txt, p_zf_retro)

  # Weighted scores. Strong sequence-family annotations dominate architecture-only evidence.
  score_nonltr <- 5*f_nonltr + 1.5*f_endo + 0.5*f_rnaseh
  score_ltr <- 5*f_ltr + 2*f_integrase + 1.5*f_rnaseh + 1*f_protease + 1*f_gag + 0.5*f_zf
  score_retrovirus <- 6*f_retrovirus + 2*f_gag + 2*f_env + 1.5*f_integrase + 1*f_protease + 1*f_rnaseh
  score_groupII <- 7*f_groupII + 2*f_groupII_support
  score_cellular <- 7*f_telomerase
  score_retron <- 7*f_retron
  score_dgr <- 7*f_dgr
  score_crispr <- 7*f_crispr

  scores <- c(
    "Retrotransposon-like: non-LTR" = score_nonltr,
    "Retrotransposon-like: LTR" = score_ltr,
    "Retrovirus-related" = score_retrovirus,
    "Group II intron-like" = score_groupII,
    "Cellular RT-containing: telomerase-like" = score_cellular,
    "Other retroelement: retron-like" = score_retron,
    "Other retroelement: DGR-like" = score_dgr,
    "Other retroelement: CRISPR-RT-like" = score_crispr
  )

  ord <- order(scores, decreasing = TRUE)
  top <- scores[ord[1]]
  second <- scores[ord[2]]
  detailed <- names(scores)[ord[1]]

  # Explicit rules for ambiguous/weak cases
  if (top < 4) {
    detailed <- "Unclassified RT-containing"
    confidence <- "low"
  } else if (top >= 7 && (top - second) >= 2) {
    confidence <- "high"
  } else if (top >= 5 && (top - second) >= 1) {
    confidence <- "moderate"
  } else {
    confidence <- "ambiguous"
    detailed <- paste0("Ambiguous: ", names(scores)[ord[1]], " vs ", names(scores)[ord[2]])
  }

  # Prefer explicit retrovirus evidence over generic LTR architecture when both are present.
  if (f_retrovirus && score_retrovirus >= 6) {
    detailed <- "Retrovirus-related"
    confidence <- if (f_gag || f_env) "high" else "moderate"
  }
  # Group II and telomerase annotations are highly diagnostic.
  if (f_groupII) { detailed <- "Group II intron-like"; confidence <- "high" }
  if (f_telomerase && !f_groupII) { detailed <- "Cellular RT-containing: telomerase-like"; confidence <- "high" }

  # Broad class requested by user
  broad <- if (grepl("^Retrotransposon-like", detailed)) {
    "Retrotransposon-like"
  } else if (grepl("^Retrovirus-related", detailed)) {
    "Retrovirus-related"
  } else if (grepl("^Group II intron-like", detailed)) {
    "Group II intron-like"
  } else if (grepl("^Cellular RT-containing", detailed)) {
    "Cellular RT-containing"
  } else if (grepl("^Other retroelement", detailed)) {
    "Other retroelement"
  } else if (grepl("^Ambiguous", detailed)) {
    "Ambiguous RT-containing"
  } else {
    "Unclassified RT-containing"
  }

  # Ordered compact architecture using informative names only
  dd <- d[order(d$start, d$end), ]
  informative <- !rx(dd$name, "^unknown$|^Extracted region from|^DNA/RNA polymerases$|^Reverse transcriptase \\(RT\\) catalytic domain profile\\.$")
  arch <- collapse_unique(dd$name[informative])

  evidence <- c(
    if (f_nonltr) "non-LTR RT signature",
    if (f_ltr) "LTR/Ty3-like RT signature",
    if (f_retrovirus) "retrovirus-related annotation",
    if (f_groupII) "group II intron/maturase signature",
    if (f_groupII_support) "group II intron RT structural similarity (supporting only)",
    if (f_telomerase) "telomerase/TERT signature",
    if (f_retron) "retron signature",
    if (f_dgr) "DGR signature",
    if (f_crispr) "CRISPR-RT signature",
    if (f_integrase) "integrase",
    if (f_rnaseh) "RNase H",
    if (f_protease) "retroviral-type protease",
    if (f_gag) "gag/capsid",
    if (f_env) "env/envelope",
    if (f_endo) "endonuclease",
    if (f_zf) "retroelement zinc-finger"
  )

  contig_len <- suppressWarnings(max(c(d$end[d$type == "extracted region"], d$end), na.rm = TRUE))
  n_2a <- sum(rx(txt, "2A_peptidase|2A peptidase|2A_peptide|2A peptide"))

  rows[[i]] <- data.frame(
    seqid = s,
    contig_length = contig_len,
    n_2A_features = n_2a,
    broad_class = broad,
    detailed_class = detailed,
    confidence = confidence,
    top_score = unname(top),
    second_score = unname(second),
    evidence = paste(evidence, collapse = "; "),
    architecture = arch,
    stringsAsFactors = FALSE
  )

  # feature-level evidence table
  ev_idx <- rx(txt, paste(p_rt_generic, p_nonltr_strong, p_ltr_strong, p_retrovirus_strong,
                         p_groupII_strong, p_groupII_support, p_telomerase, p_retron, p_dgr, p_crispr_rt,
                         p_integrase, p_rnaseh, p_retro_protease, p_gag_capsid, p_env, p_endo,
                         sep = "|"))
  if (any(ev_idx)) {
    ee <- d[ev_idx, c("seqid","source","type","start","end","strand","name")]
    ee$broad_class <- broad
    ee$detailed_class <- detailed
    ee$confidence <- confidence
    evidence_rows[[ei]] <- ee
    ei <- ei + 1L
  }
}

res <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
if (is.null(res) || nrow(res) == 0) stop("No RT-containing contigs detected with the current evidence patterns")
ev <- if (length(evidence_rows)) do.call(rbind, evidence_rows) else data.frame()

############################################################################################
################# Summary
############################################################################################

summary <- as.data.frame(table(res$broad_class, res$detailed_class, res$confidence), stringsAsFactors = FALSE)
summary <- summary[summary$Freq > 0, ]
names(summary) <- c("broad_class","detailed_class","confidence","n_contigs")
summary <- summary[order(summary$broad_class, -summary$n_contigs), ]

write.table(res, paste0(prefix, "_rt_contig_classification.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ev, paste0(prefix, "_rt_feature_evidence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, paste0(prefix, "_rt_class_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


