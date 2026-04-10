# =============================================================================
#   Per-Chromosome CNV Plot — Publication-Ready
#
#   Generates a single-chromosome log2 ratio plot with gene annotations.
#   Key sarcoma-related genes are labeled when they overlap a CNV segment.
#
#   Inputs/Outputs are controlled by Snakemake.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(R.utils)
  library(scales)
  library(ggrepel)
})

# --- Parameters from Snakemake ---
csv_file   <- snakemake@input[["csv"]]
bed_file   <- snakemake@input[["bed"]]
chr_target <- snakemake@params[["chr"]]
sample_id  <- snakemake@params[["sample"]]
pdf_file   <- snakemake@output[["pdf"]]
png_file   <- snakemake@output[["png"]]

# --- Genes of interest (hg19 coordinates) ---
# Key oncogenes and tumor suppressors relevant to sarcoma diagnostics
genes_of_interest <- data.frame(
  gene  = c("MDM2",    "CDK4",    "HMGA2",
            "MYC",     "MYCN",    "CCND1",
            "EGFR",    "ERBB2",   "CDK6",
            "RB1",     "TP53",    "CDKN2A",
            "BRCA1",   "BRCA2",   "ATM",     "TFE3"),
  chr   = c("chr12",   "chr12",   "chr12",
            "chr8",    "chr2",    "chr11",
            "chr7",    "chr17",   "chr7",
            "chr13",   "chr17",   "chr9",
            "chr17",   "chr13",   "chr11",   "chrX"),
  start = c(69202420,  58141510,  65996361,
            128748315, 16080683,  69165054,
            55086725,  37844167,  92234235,
            48303747,  7571720,   21967751,
            41196312,  32889617,  108093558, 48886237),
  end   = c(69239060,  58149796,  66301449,
            128753680, 16087129,  69178423,
            55275031,  37886679,  92837983,
            49030552,  7590868,   22009312,
            41277500,  32973809,  108239826, 48901002),
  stringsAsFactors = FALSE
)

# Filter genes on the target chromosome
genes_chr     <- genes_of_interest[genes_of_interest$chr == chr_target, ]
genes_chr$mid <- (genes_chr$start + genes_chr$end) / 2

# --- Load coverage data ---
bed <- read.table(
  gzfile(bed_file), header = FALSE, sep = "\t",
  col.names = c("Chr", "Start", "End", "Coverage")
)
if (!any(grepl("^chr", bed$Chr))) bed$Chr <- paste0("chr", bed$Chr)
bed_chr <- bed[bed$Chr == chr_target, ]
if (nrow(bed_chr) == 0) stop(paste("Chromosome", chr_target, "not found in coverage data"))

# Compute log2 ratio using genome-wide median
median_cov     <- median(bed$Coverage[bed$Coverage > 0], na.rm = TRUE)
bed_chr$Log2R  <- log2((bed_chr$Coverage + 0.01) / (median_cov + 0.01))
bed_chr$Log2R  <- pmax(pmin(bed_chr$Log2R, 5), -5)
bed_chr$MidPos <- (bed_chr$Start + bed_chr$End) / 2

# --- Load CNV calls ---
cnv_df <- read.csv(csv_file, stringsAsFactors = FALSE)
if (!any(grepl("^chr", cnv_df$Chr))) cnv_df$Chr <- paste0("chr", cnv_df$Chr)
cnv_chr <- cnv_df[cnv_df$Chr == chr_target, ]

cnv_chr$Color <- dplyr::case_when(
  grepl("DUP|gain", cnv_chr$CNV_Type, ignore.case = TRUE) ~ "#C0392B",
  grepl("DEL|loss", cnv_chr$CNV_Type, ignore.case = TRUE) ~ "#2471A3",
  TRUE ~ "#666666"
)
cnv_chr$Log2R_expected <- pmax(
  pmin(log2(pmax(cnv_chr$Copy_Number, 0.1) / 2), 5), -5
)

# --- Match genes of interest to CNV segments ---
if (nrow(genes_chr) > 0 && nrow(cnv_chr) > 0) {
  genes_chr$in_cnv <- sapply(genes_chr$gene, function(g) {
    any(grepl(paste0("\\b", g, "\\b"), cnv_chr$Genes_Overlapping))
  })
  genes_chr$cnv_y <- sapply(genes_chr$gene, function(g) {
    idx <- which(grepl(paste0("\\b", g, "\\b"), cnv_chr$Genes_Overlapping))
    if (length(idx) > 0) cnv_chr$Log2R_expected[idx[1]] else NA
  })
  genes_chr$plot_x <- sapply(genes_chr$gene, function(g) {
    idx <- which(grepl(paste0("\\b", g, "\\b"), cnv_chr$Genes_Overlapping))
    if (length(idx) > 0) {
      (cnv_chr$Start[idx[1]] + cnv_chr$End[idx[1]]) / 2
    } else {
      genes_chr$mid[genes_chr$gene == g]
    }
  })
} else if (nrow(genes_chr) > 0) {
  genes_chr$in_cnv <- rep(FALSE, nrow(genes_chr))
  genes_chr$cnv_y  <- rep(NA,    nrow(genes_chr))
  genes_chr$plot_x <- genes_chr$mid
} else {
  genes_chr$in_cnv <- logical(0)
  genes_chr$cnv_y  <- numeric(0)
  genes_chr$plot_x <- numeric(0)
}

# Only annotate genes that overlap a CNV
genes_annotate <- genes_chr[!is.na(genes_chr$in_cnv) & genes_chr$in_cnv, ]

cat("INFO:", chr_target,
    "| CNVs:", nrow(cnv_chr),
    "| Genes annotated:", nrow(genes_annotate), "\n")
if (nrow(genes_annotate) > 0) {
  cat("Genes found:", paste(genes_annotate$gene, collapse = ", "), "\n")
}

# --- Build the plot ---
make_plot <- function() {
  
  p <- ggplot() +
    
    # Reference line at 0
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5,
               linetype = "dashed") +
    
    # Coverage points
    geom_point(
      data = bed_chr,
      aes(x = MidPos, y = Log2R),
      color = "black", size = 0.15, alpha = 0.4, shape = 16
    ) +
    
    # CNV segments
    {if (nrow(cnv_chr) > 0)
      geom_rect(
        data = cnv_chr,
        aes(xmin = Start, xmax = End,
            ymin = Log2R_expected - 0.08,
            ymax = Log2R_expected + 0.08,
            fill = Color),
        alpha = 0.9
      )
    } +
    scale_fill_identity() +
    scale_color_identity() +
    
    # Gene labels with repel to avoid overlaps
    {if (nrow(genes_annotate) > 0)
      geom_label_repel(
        data              = genes_annotate,
        aes(x = plot_x, y = cnv_y, label = gene),
        fill              = "white",
        color             = "#111111",
        size              = 4.5,
        fontface          = "bold.italic",
        label.padding     = unit(0.3, "lines"),
        label.r           = unit(0.1, "lines"),
        linewidth         = 0.4,
        nudge_y           = 7.5 - genes_annotate$cnv_y,
        direction         = "x",
        segment.color     = "black",
        segment.size      = 0.5,
        segment.linetype  = "solid",
        min.segment.length = 0.1,
        box.padding       = 0.4,
        force             = 4,
        max.overlaps      = Inf,
        ylim              = c(NA, 6.5)
      )
    } +
    
    # Axes
    scale_x_continuous(
      labels = function(x) paste0(round(x / 1e6, 0), " Mb"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(-5, 6),
      breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
      labels = c("-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", "+5"),
      oob    = scales::squish
    ) +
    
    labs(
      title = paste0(sample_id, " — ", chr_target),
      x     = "Genomic position",
      y     = "Log2 copy number ratio"
    ) +
    
    theme_classic(base_size = 16) +
    theme(
      plot.title        = element_text(face = "bold", size = 16, hjust = 0,
                                       margin = margin(b = 10)),
      axis.title.x      = element_text(size = 15, color = "black",
                                       margin = margin(t = 8), face = "bold"),
      axis.title.y      = element_text(size = 15, color = "black",
                                       margin = margin(r = 8), face = "bold"),
      axis.text.x       = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y       = element_text(size = 14, color = "black", face = "bold"),
      axis.line         = element_line(color = "black", linewidth = 0.7),
      axis.ticks        = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid        = element_blank(),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(12, 20, 10, 12)
    )
  
  return(p)
}

# --- Save outputs ---
p <- make_plot()

cat("Saving PDF...\n")
ggsave(pdf_file, plot = p, width = 16, height = 7, device = "pdf")

cat("Saving PNG...\n")
ggsave(png_file, plot = p, width = 16, height = 7, dpi = 300, device = "png")

cat("Done:", chr_target, "-", sample_id, "\n")
