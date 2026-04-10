# =============================================================================
#   Genome-Wide CNV Plot
#
#   If no CNVs are detected, the script exits cleanly after generating a
#   coverage-only genome-wide figure with a clear message.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(R.utils)
})

# --- Parameters ---
csv_file      <- snakemake@input[["csv"]]
bed_file      <- snakemake@input[["bed"]]
png_file      <- snakemake@output[["png"]]
sample_id     <- snakemake@params[["sample"]]
highlight_chr <- snakemake@params[["highlight_chr"]]

if (is.null(highlight_chr) || highlight_chr == "NA" || highlight_chr == "") {
  highlight_chr <- NA
}

# =============================================================================
# Helper functions
# =============================================================================

save_png <- function(plot_obj, png_path, w = 20, h = 6) {
  ggsave(png_path, plot = plot_obj, width = w, height = h, dpi = 300, device = "png")
}

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

chr_order <- paste0("chr", c(1:22, "X", "Y"))
bed <- bed[bed$Chr %in% chr_order, , drop = FALSE]

if (nrow(bed) == 0) {
  p <- make_placeholder_plot(
    paste0("Genome-wide CNV profile — ", sample_id),
    "No valid chromosomes were found in the coverage file."
  )
  save_png(p, png_file)
  message("INFO: No valid coverage data found. Placeholder figure created.")
  quit(save = "no", status = 0)
}

bed$Chr <- factor(bed$Chr, levels = chr_order)

median_cov <- median(bed$Coverage[bed$Coverage > 0], na.rm = TRUE)

if (!is.finite(median_cov)) {
  p <- make_placeholder_plot(
    paste0("Genome-wide CNV profile — ", sample_id),
    "Coverage values are not valid for copy-number plotting."
  )
  save_png(p, png_file)
  message("INFO: Invalid coverage values. Placeholder figure created.")
  quit(save = "no", status = 0)
}

bed$Log2R <- log2((bed$Coverage + 0.01) / (median_cov + 0.01))
bed$Log2R <- pmax(pmin(bed$Log2R, 3), -3)

# =============================================================================
# Prepare genome-wide coordinates
# =============================================================================

chr_sizes <- bed %>%
  group_by(Chr) %>%
  summarise(Size = max(End), .groups = "drop") %>%
  arrange(Chr)

chr_sizes$Offset <- cumsum(c(0, head(chr_sizes$Size, -1)))
chr_sizes$MidPos <- chr_sizes$Offset + chr_sizes$Size / 2

bed2 <- bed %>%
  left_join(chr_sizes[, c("Chr", "Offset")], by = "Chr") %>%
  mutate(GenomicPos = (Start + End) / 2 + Offset)

chr_bands <- chr_sizes %>%
  mutate(
    xmin = Offset,
    xmax = Offset + Size,
    fill = ifelse(seq_len(n()) %% 2 == 0, "#C7C7C7", "#FFFFFF")
  )

step <- max(1, round(nrow(bed2) / 80000))
bed_plot <- bed2[seq(1, nrow(bed2), by = step), , drop = FALSE]

# =============================================================================
# Load CNV calls safely
# =============================================================================

cnv_df <- safe_read_cnv(csv_file)

if (nrow(cnv_df) == 0) {
  p <- ggplot() +
    geom_rect(
      data = chr_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -3.5, ymax = 3.5, fill = fill),
      alpha = 0.6, inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_point(
      data = bed_plot,
      aes(x = GenomicPos, y = Log2R),
      size = 0.15, alpha = 0.3, color = "black"
    ) +
    scale_x_continuous(
      breaks = chr_sizes$MidPos,
      labels = gsub("chr", "", chr_sizes$Chr),
      expand = c(0.005, 0)
    ) +
    scale_y_continuous(
      limits = c(-3.5, 3.5),
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      labels = c("-3", "-2", "-1", "0", "+1", "+2", "+3")
    ) +
    labs(
      title = paste0("Genome-wide CNV profile — ", sample_id),
      subtitle = "No CNVs were detected",
      x = "Chromosome",
      y = "Log2 copy number ratio"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 14, color = "#555555", margin = margin(b = 10)),
      axis.title.x = element_text(size = 15, color = "black", margin = margin(t = 8), face = "bold"),
      axis.title.y = element_text(size = 15, color = "black", margin = margin(r = 8), face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 20, 10, 12)
    )
  
  save_png(p, png_file)
  message("INFO: No CNVs were detected. Coverage-only genome-wide figure created.")
  quit(save = "no", status = 0)
}

required_cols <- c("Chr", "Start", "End", "CNV_Type", "Copy_Number")
missing_cols <- setdiff(required_cols, colnames(cnv_df))

if (length(missing_cols) > 0) {
  p <- make_placeholder_plot(
    paste0("Genome-wide CNV profile — ", sample_id),
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
cnv_df$CNV_Type <- as.character(cnv_df$CNV_Type)

cnv_df <- cnv_df[cnv_df$Chr %in% chr_order, , drop = FALSE]

if (nrow(cnv_df) == 0) {
  p <- ggplot() +
    geom_rect(
      data = chr_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -3.5, ymax = 3.5, fill = fill),
      alpha = 0.6, inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_point(
      data = bed_plot,
      aes(x = GenomicPos, y = Log2R),
      size = 0.15, alpha = 0.3, color = "black"
    ) +
    scale_x_continuous(
      breaks = chr_sizes$MidPos,
      labels = gsub("chr", "", chr_sizes$Chr),
      expand = c(0.005, 0)
    ) +
    scale_y_continuous(
      limits = c(-3.5, 3.5),
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      labels = c("-3", "-2", "-1", "0", "+1", "+2", "+3")
    ) +
    labs(
      title = paste0("Genome-wide CNV profile — ", sample_id),
      subtitle = "No valid CNVs were detected",
      x = "Chromosome",
      y = "Log2 copy number ratio"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 14, color = "#555555", margin = margin(b = 10)),
      axis.title.x = element_text(size = 15, color = "black", margin = margin(t = 8), face = "bold"),
      axis.title.y = element_text(size = 15, color = "black", margin = margin(r = 8), face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 20, 10, 12)
    )
  
  save_png(p, png_file)
  message("INFO: No valid CNVs were detected. Coverage-only genome-wide figure created.")
  quit(save = "no", status = 0)
}

cnv_df$Chr <- factor(cnv_df$Chr, levels = chr_order)

cnv_df$Color <- dplyr::case_when(
  grepl("gain|DUP", cnv_df$CNV_Type, ignore.case = TRUE) ~ "#C0392B",
  grepl("loss|DEL", cnv_df$CNV_Type, ignore.case = TRUE) ~ "#2471A3",
  grepl("loh|LOH",  cnv_df$CNV_Type, ignore.case = TRUE) ~ "#E67E22",
  TRUE ~ "#666666"
)

cnv2 <- cnv_df %>%
  left_join(chr_sizes[, c("Chr", "Offset")], by = "Chr") %>%
  mutate(
    GenomicStart = Start + Offset,
    GenomicEnd   = End + Offset,
    Log2R_exp    = pmax(pmin(log2(pmax(Copy_Number, 0.1) / 2), 2.8), -2.8)
  )

n_gain <- sum(grepl("gain|DUP", cnv_df$CNV_Type, ignore.case = TRUE))
n_loss <- sum(grepl("loss|DEL", cnv_df$CNV_Type, ignore.case = TRUE))
n_loh  <- sum(grepl("loh|LOH",  cnv_df$CNV_Type, ignore.case = TRUE))

p <- ggplot() +
  geom_rect(
    data = chr_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -3.5, ymax = 3.5, fill = fill),
    alpha = 0.6, inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 1, color = "#C0392B", linewidth = 0.3, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -1, color = "#2471A3", linewidth = 0.3, linetype = "dashed", alpha = 0.5) +
  geom_point(
    data = bed_plot,
    aes(x = GenomicPos, y = Log2R),
    size = 0.15, alpha = 0.3, color = "black"
  ) +
  geom_rect(
    data = cnv2,
    aes(
      xmin = GenomicStart, xmax = GenomicEnd,
      ymin = Log2R_exp - 0.12, ymax = Log2R_exp + 0.12,
      fill = Color
    ),
    alpha = 0.95, inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = chr_sizes$MidPos,
    labels = gsub("chr", "", chr_sizes$Chr),
    expand = c(0.005, 0)
  ) +
  scale_y_continuous(
    limits = c(-3.5, 3.5),
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c("-3", "-2", "-1", "0", "+1", "+2", "+3")
  ) +
  labs(
    title = paste0("Genome-wide CNV profile — ", sample_id),
    subtitle = paste0(
      nrow(cnv_df), " CNVs | ",
      n_gain, " gains | ",
      n_loss, " losses | ",
      n_loh, " LOH"
    ),
    x = "Chromosome",
    y = "Log2 copy number ratio"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 14, color = "#555555", margin = margin(b = 10)),
    axis.title.x = element_text(size = 15, color = "black", margin = margin(t = 8), face = "bold"),
    axis.title.y = element_text(size = 15, color = "black", margin = margin(r = 8), face = "bold"),
    axis.text.x = element_text(size = 13, color = "black", face = "bold"),
    axis.text.y = element_text(size = 13, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 20, 10, 12)
  )

if (!is.na(highlight_chr)) {
  hchr <- ifelse(grepl("^chr", highlight_chr), highlight_chr, paste0("chr", highlight_chr))
  hdata <- chr_sizes[chr_sizes$Chr == hchr, , drop = FALSE]
  if (nrow(hdata) > 0) {
    p <- p + annotate(
      "rect",
      xmin = hdata$Offset,
      xmax = hdata$Offset + hdata$Size,
      ymin = -3.5,
      ymax = 3.5,
      fill = "#FF5733",
      alpha = 0.08
    )
  }
}

save_png(p, png_file)
message("INFO: Genome-wide plot completed successfully.")
quit(save = "no", status = 0)