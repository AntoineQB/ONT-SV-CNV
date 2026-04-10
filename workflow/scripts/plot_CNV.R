# =============================================================================
#   Genome-Wide CNV Plot — Publication-Ready
#
#   Generates a genome-wide log2 ratio plot from mosdepth coverage bins,
#   with CNV segments overlaid. Optionally highlights a chromosome of interest.
#
#   Inputs/Outputs are controlled by Snakemake.
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
pdf_file      <- snakemake@output[["pdf"]]
png_file      <- snakemake@output[["png"]]
sample_id     <- snakemake@params[["sample"]]
highlight_chr <- snakemake@params[["highlight_chr"]]

if (is.null(highlight_chr) || highlight_chr == "NA" || highlight_chr == "") {
  highlight_chr <- NA
}

# --- Load coverage data ---
bed <- read.table(
  gzfile(bed_file), header = FALSE, sep = "\t",
  col.names = c("Chr", "Start", "End", "Coverage")
)
if (!any(grepl("^chr", bed$Chr))) bed$Chr <- paste0("chr", bed$Chr)

chr_order <- paste0("chr", c(1:22, "X", "Y"))
bed       <- bed[bed$Chr %in% chr_order, ]
bed$Chr   <- factor(bed$Chr, levels = chr_order)

# Compute log2 ratio relative to median coverage
median_cov <- median(bed$Coverage[bed$Coverage > 0], na.rm = TRUE)
bed$Log2R  <- log2((bed$Coverage + 0.01) / (median_cov + 0.01))
bed$Log2R  <- pmax(pmin(bed$Log2R, 3), -3)

# --- Load CNV calls ---
cnv_df <- read.csv(csv_file, stringsAsFactors = FALSE)
if (!any(grepl("^chr", cnv_df$Chr))) cnv_df$Chr <- paste0("chr", cnv_df$Chr)
cnv_df$Chr <- factor(cnv_df$Chr, levels = chr_order)

cnv_df$Color <- dplyr::case_when(
  grepl("gain|DUP", cnv_df$CNV_Type, ignore.case = TRUE) ~ "#C0392B",
  grepl("loss|DEL", cnv_df$CNV_Type, ignore.case = TRUE) ~ "#2471A3",
  grepl("loh|LOH",  cnv_df$CNV_Type, ignore.case = TRUE) ~ "#E67E22",
  TRUE ~ "#666666"
)

# --- Build the plot ---
make_cnv_plot <- function() {
  
  # Compute genomic offsets for linear genome layout
  chr_sizes <- bed %>%
    group_by(Chr) %>%
    summarise(Size = max(End), .groups = "drop") %>%
    arrange(Chr)
  chr_sizes$Offset <- cumsum(c(0, head(chr_sizes$Size, -1)))
  chr_sizes$MidPos <- chr_sizes$Offset + chr_sizes$Size / 2
  
  bed2 <- bed %>%
    left_join(chr_sizes[, c("Chr", "Offset")], by = "Chr") %>%
    mutate(GenomicPos = (Start + End) / 2 + Offset)
  
  cnv2 <- cnv_df %>%
    left_join(chr_sizes[, c("Chr", "Offset")], by = "Chr") %>%
    mutate(
      GenomicStart = Start + Offset,
      GenomicEnd   = End + Offset,
      Log2R_exp    = pmax(pmin(log2(pmax(Copy_Number, 0.1) / 2), 2.8), -2.8)
    )
  
  # Alternating chromosome background bands
  chr_bands <- chr_sizes %>%
    mutate(
      xmin = Offset,
      xmax = Offset + Size,
      fill = ifelse(as.integer(Chr) %% 2 == 0, "#C7C7C7", "#FFFFFF")
    )
  
  # Subsample coverage points for rendering speed
  step     <- max(1, round(nrow(bed2) / 80000))
  bed_plot <- bed2[seq(1, nrow(bed2), by = step), ]
  
  # Count CNV types for subtitle
  n_gain <- sum(grepl("gain|DUP", cnv_df$CNV_Type, ignore.case = TRUE))
  n_loss <- sum(grepl("loss|DEL", cnv_df$CNV_Type, ignore.case = TRUE))
  n_loh  <- sum(grepl("loh|LOH",  cnv_df$CNV_Type, ignore.case = TRUE))
  
  p <- ggplot() +
    
    # Alternating chromosome bands
    geom_rect(
      data = chr_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -3.5, ymax = 3.5, fill = fill),
      alpha = 0.6, inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    
    # Reference lines
    geom_hline(yintercept = 0,  color = "black",   linewidth = 0.5) +
    geom_hline(yintercept =  1, color = "#C0392B", linewidth = 0.3,
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -1, color = "#2471A3", linewidth = 0.3,
               linetype = "dashed", alpha = 0.5) +
    
    # Coverage points
    geom_point(
      data = bed_plot,
      aes(x = GenomicPos, y = Log2R),
      size = 0.15, alpha = 0.3, color = "black"
    ) +
    
    # CNV segments
    geom_rect(
      data = cnv2,
      aes(xmin = GenomicStart, xmax = GenomicEnd,
          ymin = Log2R_exp - 0.12, ymax = Log2R_exp + 0.12,
          fill = Color),
      alpha = 0.95, inherit.aes = FALSE
    ) +
    
    # Optional chromosome highlight
    {if (!is.na(highlight_chr)) {
      hchr  <- ifelse(grepl("^chr", highlight_chr), highlight_chr,
                      paste0("chr", highlight_chr))
      hdata <- chr_sizes[chr_sizes$Chr == hchr, ]
      if (nrow(hdata) > 0)
        annotate("rect",
                 xmin = hdata$Offset, xmax = hdata$Offset + hdata$Size,
                 ymin = -3.5, ymax = 3.5,
                 fill = "#FF5733", alpha = 0.08)
    }} +
    
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
      title    = paste0("Genome-wide CNV profile — ", sample_id),
      subtitle = paste0(nrow(cnv_df), " CNVs | ",
                        n_gain, " gains | ",
                        n_loss, " losses | ",
                        n_loh,  " LOH"),
      x = "Chromosome",
      y = "Log2 copy number ratio"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title        = element_text(face = "bold", size = 16, hjust = 0,
                                       margin = margin(b = 4)),
      plot.subtitle     = element_text(size = 14, color = "#555555",
                                       margin = margin(b = 10)),
      axis.title.x      = element_text(size = 15, color = "black",
                                       margin = margin(t = 8), face = "bold"),
      axis.title.y      = element_text(size = 15, color = "black",
                                       margin = margin(r = 8), face = "bold"),
      axis.text.x       = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y       = element_text(size = 13, color = "black", face = "bold"),
      axis.line         = element_line(color = "black", linewidth = 0.6),
      axis.ticks        = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid        = element_blank(),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(12, 20, 10, 12)
    )
  
  return(p)
}

# --- Save outputs ---
cat("Saving PDF...\n")
p <- make_cnv_plot()
ggsave(pdf_file, plot = p, width = 20, height = 6, device = "pdf")

cat("Saving PNG...\n")
ggsave(png_file, plot = p, width = 20, height = 6, dpi = 300, device = "png")

cat("Done:", sample_id, "\n")
