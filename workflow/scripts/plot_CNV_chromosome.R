# =============================================================================
#   Per-Chromosome CNV Plot
#
#   If no CNVs are detected, the script exits cleanly after generating a simple
#   placeholder or coverage-only figure.
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
png_file   <- snakemake@output[["png"]]

# --- Genes of interest (hg19 coordinates) ---
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

# =============================================================================
# Helper functions
# =============================================================================

make_placeholder_plot <- function(title_text, message_text) {
  ggplot() +
    annotate("text", x = 0, y = 0.15, label = title_text,
             fontface = "bold", size = 7) +
    annotate("text", x = 0, y = -0.05, label = message_text,
             size = 5.5, color = "grey30") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

save_png <- function(plot_obj, png_path, w = 16, h = 7) {
  ggsave(png_path, plot = plot_obj, width = w, height = h, dpi = 300, device = "png")
}

safe_read_cnv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame())
  }
  
  out <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE),
    error = function(e) {
      message("WARNING: Could not read CNV CSV file. Treating as empty.")
      data.frame()
    }
  )
  
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame())
  }
  
  out
}

# =============================================================================
# Load coverage data
# =============================================================================

bed <- read.table(
  gzfile(bed_file), header = FALSE, sep = "\t",
  col.names = c("Chr", "Start", "End", "Coverage")
)

bed$Chr <- as.character(bed$Chr)
if (!any(grepl("^chr", bed$Chr))) {
  bed$Chr <- paste0("chr", bed$Chr)
}

bed_chr <- bed[bed$Chr == chr_target, , drop = FALSE]

if (nrow(bed_chr) == 0) {
  p <- make_placeholder_plot(
    paste0(sample_id, " — ", chr_target),
    paste("Chromosome", chr_target, "was not found in the coverage file.")
  )
  save_png(p, png_file)
  message("INFO: Chromosome not found in coverage data. Placeholder figure created.")
  quit(save = "no", status = 0)
}

median_cov <- median(bed$Coverage[bed$Coverage > 0], na.rm = TRUE)

if (!is.finite(median_cov)) {
  p <- make_placeholder_plot(
    paste0(sample_id, " — ", chr_target),
    "Coverage values are not valid for copy-number plotting."
  )
  save_png(p, png_file)
  message("INFO: Invalid coverage values. Placeholder figure created.")
  quit(save = "no", status = 0)
}

bed_chr$Log2R  <- log2((bed_chr$Coverage + 0.01) / (median_cov + 0.01))
bed_chr$Log2R  <- pmax(pmin(bed_chr$Log2R, 5), -5)
bed_chr$MidPos <- (bed_chr$Start + bed_chr$End) / 2

# =============================================================================
# Load CNV calls safely
# =============================================================================

cnv_df <- safe_read_cnv(csv_file)

if (nrow(cnv_df) == 0) {
  p <- ggplot(bed_chr, aes(x = MidPos, y = Log2R)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_point(color = "black", size = 0.15, alpha = 0.4, shape = 16) +
    scale_x_continuous(
      labels = function(x) paste0(round(x / 1e6, 0), " Mb"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(-5, 5),
      breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
      labels = c("-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", "+5"),
      oob = scales::squish
    ) +
    labs(
      title = paste0(sample_id, " — ", chr_target),
      subtitle = "No CNVs were detected",
      x = "Genomic position",
      y = "Log2 copy number ratio"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = 13, color = "grey30", margin = margin(b = 10)),
      axis.title.x = element_text(size = 15, color = "black", margin = margin(t = 8), face = "bold"),
      axis.title.y = element_text(size = 15, color = "black", margin = margin(r = 8), face = "bold"),
      axis.text.x = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 14, color = "black", face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 20, 10, 12)
    )
  
  save_png(p, png_file)
  message("INFO: No CNVs were detected. Coverage-only figure created.")
  quit(save = "no", status = 0)
}

required_cols <- c("Chr", "Start", "End", "CNV_Type", "Copy_Number", "Genes_Overlapping")
missing_cols <- setdiff(required_cols, colnames(cnv_df))

if (length(missing_cols) > 0) {
  p <- make_placeholder_plot(
    paste0(sample_id, " — ", chr_target),
    paste0("CNV table is missing required columns: ", paste(missing_cols, collapse = ", "))
  )
  save_png(p, png_file)
  message("INFO: CNV CSV is malformed. Placeholder figure created.")
  quit(save = "no", status = 0)
}

cnv_df$Chr <- as.character(cnv_df$Chr)
if (!any(grepl("^chr", cnv_df$Chr))) {
  cnv_df$Chr <- paste0("chr", cnv_df$Chr)
}

cnv_df$Start <- suppressWarnings(as.numeric(cnv_df$Start))
cnv_df$End <- suppressWarnings(as.numeric(cnv_df$End))
cnv_df$Copy_Number <- suppressWarnings(as.numeric(cnv_df$Copy_Number))
cnv_df$Genes_Overlapping <- as.character(cnv_df$Genes_Overlapping)
cnv_df$CNV_Type <- as.character(cnv_df$CNV_Type)

cnv_chr <- cnv_df[cnv_df$Chr == chr_target, , drop = FALSE]

if (nrow(cnv_chr) == 0) {
  p <- ggplot(bed_chr, aes(x = MidPos, y = Log2R)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_point(color = "black", size = 0.15, alpha = 0.4, shape = 16) +
    scale_x_continuous(
      labels = function(x) paste0(round(x / 1e6, 0), " Mb"),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(-5, 5),
      breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
      labels = c("-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", "+5"),
      oob = scales::squish
    ) +
    labs(
      title = paste0(sample_id, " — ", chr_target),
      subtitle = "No CNVs were detected on this chromosome",
      x = "Genomic position",
      y = "Log2 copy number ratio"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = 13, color = "grey30", margin = margin(b = 10)),
      axis.title.x = element_text(size = 15, color = "black", margin = margin(t = 8), face = "bold"),
      axis.title.y = element_text(size = 15, color = "black", margin = margin(r = 8), face = "bold"),
      axis.text.x = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 14, color = "black", face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 20, 10, 12)
    )
  
  save_png(p, png_file)
  message("INFO: No CNVs were detected on this chromosome. Coverage-only figure created.")
  quit(save = "no", status = 0)
}

# =============================================================================
# Continue normal plotting
# =============================================================================

cnv_chr$Color <- dplyr::case_when(
  grepl("DUP|gain", cnv_chr$CNV_Type, ignore.case = TRUE) ~ "#C0392B",
  grepl("DEL|loss", cnv_chr$CNV_Type, ignore.case = TRUE) ~ "#2471A3",
  TRUE ~ "#666666"
)

cnv_chr$Log2R_expected <- pmax(
  pmin(log2(pmax(cnv_chr$Copy_Number, 0.1) / 2), 5), -5
)

genes_chr <- genes_of_interest[genes_of_interest$chr == chr_target, , drop = FALSE]
if (nrow(genes_chr) > 0) {
  genes_chr$mid <- (genes_chr$start + genes_chr$end) / 2
}

if (nrow(genes_chr) > 0) {
  genes_chr$in_cnv <- sapply(genes_chr$gene, function(g) {
    any(grepl(paste0("\\b", g, "\\b"), cnv_chr$Genes_Overlapping))
  })
  
  genes_chr$cnv_y <- sapply(genes_chr$gene, function(g) {
    idx <- which(grepl(paste0("\\b", g, "\\b"), cnv_chr$Genes_Overlapping))
    if (length(idx) > 0) cnv_chr$Log2R_expected[idx[1]] else NA_real_
  })
  
  genes_chr$plot_x <- sapply(genes_chr$gene, function(g) {
    idx <- which(grepl(paste0("\\b", g, "\\b"), cnv_chr$Genes_Overlapping))
    if (length(idx) > 0) {
      (cnv_chr$Start[idx[1]] + cnv_chr$End[idx[1]]) / 2
    } else {
      genes_chr$mid[genes_chr$gene == g]
    }
  })
  
  genes_annotate <- genes_chr[!is.na(genes_chr$in_cnv) & genes_chr$in_cnv, , drop = FALSE]
} else {
  genes_annotate <- data.frame()
}

p <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_point(
    data = bed_chr,
    aes(x = MidPos, y = Log2R),
    color = "black", size = 0.15, alpha = 0.4, shape = 16
  ) +
  geom_rect(
    data = cnv_chr,
    aes(
      xmin = Start, xmax = End,
      ymin = Log2R_expected - 0.08,
      ymax = Log2R_expected + 0.08,
      fill = Color
    ),
    alpha = 0.9
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    labels = function(x) paste0(round(x / 1e6, 0), " Mb"),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(-5, 6),
    breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
    labels = c("-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", "+5"),
    oob = scales::squish
  ) +
  labs(
    title = paste0(sample_id, " — ", chr_target),
    x = "Genomic position",
    y = "Log2 copy number ratio"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(b = 10)),
    axis.title.x = element_text(size = 15, color = "black", margin = margin(t = 8), face = "bold"),
    axis.title.y = element_text(size = 15, color = "black", margin = margin(r = 8), face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 20, 10, 12)
  )

if (nrow(genes_annotate) > 0) {
  p <- p + geom_label_repel(
    data = genes_annotate,
    aes(x = plot_x, y = cnv_y, label = gene),
    fill = "white",
    color = "#111111",
    size = 4.5,
    fontface = "bold.italic",
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.1, "lines"),
    linewidth = 0.4,
    nudge_y = 7.5 - genes_annotate$cnv_y,
    direction = "x",
    segment.color = "black",
    segment.size = 0.5,
    segment.linetype = "solid",
    min.segment.length = 0.1,
    box.padding = 0.4,
    force = 4,
    max.overlaps = Inf,
    ylim = c(NA, 6.5)
  )
}

save_png(p, png_file)
message("INFO: Chromosome plot completed successfully.")
quit(save = "no", status = 0)