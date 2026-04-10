# =============================================================================
#   CNV Annotation with GENCODE v19 (hg19)
#
#   Annotates copy number variants from a Spectre VCF with overlapping genes
#   from the GENCODE v19 annotation. Prioritizes protein-coding genes.
#
#   Inputs/Outputs are controlled by Snakemake.
# =============================================================================

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(dplyr)
  library(rtracklayer)
})

# --- 1. Snakemake parameters ---
vcf_file   <- snakemake@input[["vcf"]]
gtf_file   <- snakemake@input[["gtf"]]
output_csv <- snakemake@output[["csv"]]
sample_id  <- snakemake@params[["sample"]]

cat("Reading VCF:", vcf_file, "\n")
cat("Reading GTF:", gtf_file, "\n")

# --- 2. Read Spectre VCF ---
vcf     <- readVcf(vcf_file)
info_df <- info(vcf)

chrom <- as.character(seqnames(rowRanges(vcf)))
start <- start(rowRanges(vcf))
sv_id <- names(rowRanges(vcf))

# Extract END coordinate from INFO field (Spectre stores it there)
if ("END" %in% colnames(info_df)) {
  end <- as.numeric(info_df$END)
} else if ("SVLEN" %in% colnames(info_df)) {
  end <- start + abs(as.numeric(info_df$SVLEN)) - 1
} else {
  end <- start
  warning("Neither END nor SVLEN found in VCF — end coordinates will be incorrect")
}

# CNV type
cnv_type <- as.character(info_df$SVTYPE)
if ("CNVTYPE" %in% colnames(info_df)) {
  cnv_type <- as.character(info_df$CNVTYPE)
}

# Copy number
copy_number <- rep(NA_real_, nrow(info_df))
if ("CN" %in% colnames(info_df)) {
  copy_number <- as.numeric(info_df$CN)
}

# CNV length
cnv_length <- end - start + 1

cnv_df <- data.frame(
  SV_ID       = sv_id,
  Chr         = chrom,
  Start       = start,
  End         = end,
  Length_bp   = cnv_length,
  CNV_Type    = cnv_type,
  Copy_Number = copy_number,
  stringsAsFactors = FALSE
)

cat(nrow(cnv_df), "CNVs detected before annotation\n")
cat("  Gains :", sum(grepl("gain|DUP", cnv_df$CNV_Type, ignore.case = TRUE)), "\n")
cat("  Losses:", sum(grepl("loss|DEL", cnv_df$CNV_Type, ignore.case = TRUE)), "\n")
cat("  LOH   :", sum(grepl("loh|LOH",  cnv_df$CNV_Type, ignore.case = TRUE)), "\n")

# --- 3. Load GENCODE v19 ---
cat("Loading GENCODE v19 annotation...\n")
gtf        <- rtracklayer::import(gtf_file)
genes_gtf  <- gtf[gtf$type == "gene"]
gene_ranges <- GRanges(
  seqnames  = seqnames(genes_gtf),
  ranges    = ranges(genes_gtf),
  strand    = strand(genes_gtf),
  gene_id   = genes_gtf$gene_id,
  gene_name = genes_gtf$gene_name,
  gene_type = genes_gtf$gene_type
)

# --- 4. Find genes overlapping each CNV ---
find_overlapping_genes <- function(chr, start, end, max_genes = 20) {
  if (is.na(chr) || is.na(start) || is.na(end)) {
    return(list(names = "intergenic", count = 0))
  }
  
  gr   <- GRanges(seqnames = chr, ranges = IRanges(start, end))
  hits <- findOverlaps(gr, gene_ranges)
  
  if (length(hits) == 0) return(list(names = "intergenic", count = 0))
  
  matched    <- gene_ranges[subjectHits(hits)]
  gene_names <- unique(na.omit(mcols(matched)$gene_name))
  
  # Prioritize protein-coding genes
  protein_coding <- gene_names[
    mcols(matched)$gene_type[
      match(gene_names, mcols(matched)$gene_name)
    ] == "protein_coding"
  ]
  protein_coding <- unique(na.omit(protein_coding))
  
  # Limit display length
  if (length(protein_coding) > max_genes) {
    genes_label <- paste0(
      paste(protein_coding[1:max_genes], collapse = ", "),
      " [+", length(protein_coding) - max_genes, " more]"
    )
  } else if (length(protein_coding) > 0) {
    genes_label <- paste(protein_coding, collapse = ", ")
  } else {
    genes_label <- paste(
      gene_names[1:min(5, length(gene_names))], collapse = ", "
    )
  }
  
  list(names = genes_label, count = length(gene_names))
}

# --- 5. Apply annotation ---
cat("Annotating CNVs...\n")
annotation_results <- mapply(
  find_overlapping_genes,
  cnv_df$Chr, cnv_df$Start, cnv_df$End,
  SIMPLIFY = FALSE
)

cnv_df$Genes_Overlapping   <- sapply(annotation_results, function(x) x$names)
cnv_df$N_Genes_Overlapping <- sapply(annotation_results, function(x) x$count)

# --- 6. Add length in Mb for readability ---
cnv_df$Length_Mb <- round(cnv_df$Length_bp / 1e6, 3)

# --- 7. Sort by chromosome and position ---
chr_order <- paste0("chr", c(1:22, "X", "Y"))
cnv_df$Chr <- factor(cnv_df$Chr, levels = chr_order)
cnv_df     <- cnv_df[order(cnv_df$Chr, cnv_df$Start), ]
cnv_df$Chr <- as.character(cnv_df$Chr)

# --- 8. Save ---
write.csv(cnv_df, output_csv, row.names = FALSE)
cat("Annotation complete:", nrow(cnv_df), "CNVs annotated\n")
cat("Results saved to:", output_csv, "\n")
