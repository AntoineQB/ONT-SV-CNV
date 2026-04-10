# =============================================================================
#   Breakpoint (BND) Annotation with GENCODE v19 (hg19)
#
#   Annotates structural variant breakpoints from a VCF file with overlapping
#   or nearby genes from the GENCODE v19 GTF annotation.
#   Compatible with both Delly (CHR2 + POS2) and Sniffles2 (CHR2 + END).
#
# =============================================================================

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(dplyr)
  library(rtracklayer)
  # txdbmaker requis depuis GenomicFeatures >= 1.61.1
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    library(txdbmaker)
  }
})

# --- 1. Retrieve paths from Snakemake ---
vcf_file   <- snakemake@input[[1]]
gtf_file   <- snakemake@input[[2]]
output_csv <- snakemake@output[[1]]

cat("Reading VCF file:", vcf_file, "\n")
cat("Reading GTF file:", gtf_file, "\n")

# --- 2. Read VCF ---
vcf <- readVcf(vcf_file)

info_df <- info(vcf)
chrom   <- as.character(seqnames(rowRanges(vcf)))
pos1    <- start(rowRanges(vcf))
sv_id   <- names(rowRanges(vcf))

# --- 3. Extract Chr2 / Pos2 (compatible with Sniffles2 and Delly) ---
#
#   Sniffles2 : CHR2 in INFO, position in END (not POS2)
#   Delly     : CHR2 in INFO, position in POS2
#
#   For intra-chromosomal events (INV/DUP), CHR2 may be absent;
#   in that case, we fall back to Chr1.

extract_chr2 <- function(info_df, chrom) {
  if ("CHR2" %in% colnames(info_df)) {
    chr2 <- as.character(info_df$CHR2)
    chr2 <- ifelse(is.na(chr2) | chr2 == "", chrom, chr2)
  } else {
    chr2 <- chrom
  }
  chr2
}

extract_pos2 <- function(info_df, pos1) {
  if ("POS2" %in% colnames(info_df)) {
    pos2 <- as.integer(info_df$POS2)
  } else if ("END" %in% colnames(info_df)) {
    pos2 <- as.integer(info_df$END)
  } else {
    pos2 <- rep(NA_integer_, length(pos1))
  }
  pos2 <- ifelse(is.na(pos2), pos1, pos2)
  pos2
}

chr2 <- extract_chr2(info_df, chrom)
pos2 <- extract_pos2(info_df, pos1)

cat("Available INFO columns:", paste(colnames(info_df), collapse = ", "), "\n")
cat("Number of variants read:", length(sv_id), "\n")

# --- 4. Build SV data frame ---
sv_df <- data.frame(
  SV_ID = sv_id,
  Chr1  = chrom,
  Pos1  = pos1,
  Chr2  = chr2,
  Pos2  = pos2,
  stringsAsFactors = FALSE
)

# --- 5. Load GENCODE v19 gene annotations ---
cat("Loading GENCODE GTF annotation...\n")
gtf       <- rtracklayer::import(gtf_file)
genes_gtf <- gtf[gtf$type == "gene"]

gene_ranges <- GRanges(
  seqnames  = seqnames(genes_gtf),
  ranges    = ranges(genes_gtf),
  strand    = strand(genes_gtf),
  gene_id   = genes_gtf$gene_id,
  gene_name = genes_gtf$gene_name
)

# --- 6. Find overlapping or nearby genes ---
find_gene <- function(chr, pos, window = 50000) {
  if (is.na(chr) || is.na(pos)) return(NA_character_)
  
  # First, look for direct overlaps
  gr   <- GRanges(seqnames = chr, ranges = IRanges(pos, pos))
  hits <- findOverlaps(gr, gene_ranges)
  
  # If no direct overlap, extend the search window
  if (length(hits) == 0) {
    gr_window <- GRanges(
      seqnames = chr,
      ranges   = IRanges(max(1, pos - window), pos + window)
    )
    hits <- findOverlaps(gr_window, gene_ranges)
  }
  
  if (length(hits) == 0) return(NA_character_)
  
  gene_names <- unique(mcols(gene_ranges)$gene_name[subjectHits(hits)])
  paste(unique(na.omit(gene_names)), collapse = ", ")
}

# --- 7. Annotate each breakpoint ---
cat("Annotating breakpoints...\n")
sv_df$Gene1 <- mapply(find_gene, sv_df$Chr1, sv_df$Pos1)
sv_df$Gene2 <- mapply(find_gene, sv_df$Chr2, sv_df$Pos2)

sv_df$Gene1[is.na(sv_df$Gene1)] <- "intergenic"
sv_df$Gene2[is.na(sv_df$Gene2)] <- "intergenic"

# Remove rows with missing coordinates
sv_df <- sv_df[!is.na(sv_df$Chr1) & !is.na(sv_df$Pos1) &
               !is.na(sv_df$Chr2) & !is.na(sv_df$Pos2), ]

cat(nrow(sv_df), "variants annotated and saved to:", output_csv, "\n")

# --- 8. Save output ---
write.csv(sv_df, output_csv, row.names = FALSE)
