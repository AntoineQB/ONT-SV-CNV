# =============================================================================
#   Circos Plot — Publication-Ready
#
#   Generates a circos plot of structural variant breakpoints (BND/TRA/INV/DUP)
#   with optional highlighting of links between two chromosomes of interest.
#
#   Inputs/Outputs are controlled by Snakemake.
# =============================================================================

suppressPackageStartupMessages({
  library(circlize)
  library(RColorBrewer)
})

# --- 1. Parameters from Snakemake ---
sample       <- snakemake@wildcards[["sample"]]
sample_label <- snakemake@params[["sample_label"]]

highlight_chr1 <- snakemake@params[["highlight_chr1"]]
highlight_chr2 <- snakemake@params[["highlight_chr2"]]

# Set to NA if not provided or invalid
if (is.null(highlight_chr1) || highlight_chr1 == "" || highlight_chr1 == "NA" ||
    is.null(highlight_chr2) || highlight_chr2 == "" || highlight_chr2 == "NA") {
  highlight_chr1 <- NA
  highlight_chr2 <- NA
}

cat("Sample:", sample, "| Label:", sample_label, "\n")
if (!is.na(highlight_chr1) & !is.na(highlight_chr2)) {
  cat("Highlighted chromosomes:", highlight_chr1, "and", highlight_chr2, "\n")
} else {
  cat("No chromosome highlighting — all links drawn uniformly\n")
}

# --- 2. Read annotated CSV ---
input_file <- snakemake@input[["csv"]]
cat("Reading:", input_file, "\n")
sv_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# --- 3. Harmonize column names ---
if (!all(c("chrom1", "pos1", "chrom2", "pos2") %in% names(sv_data))) {
  sv_data$chrom1 <- sv_data$Chr1
  sv_data$pos1   <- as.numeric(sv_data$Pos1)
  sv_data$chrom2 <- sv_data$Chr2
  sv_data$pos2   <- as.numeric(sv_data$Pos2)
}

# --- 4. Ensure "chr" prefix on chromosome names ---
if (!any(grepl("^chr", sv_data$chrom1))) {
  sv_data$chrom1 <- paste0("chr", sv_data$chrom1)
  sv_data$chrom2 <- paste0("chr", sv_data$chrom2)
}

if (!is.na(highlight_chr1)) {
  highlight_chr1 <- ifelse(
    grepl("^chr", highlight_chr1), highlight_chr1, paste0("chr", highlight_chr1)
  )
}
if (!is.na(highlight_chr2)) {
  highlight_chr2 <- ifelse(
    grepl("^chr", highlight_chr2), highlight_chr2, paste0("chr", highlight_chr2)
  )
}

# --- 5. Output files ---
pdf_file <- snakemake@output[["pdf"]]
png_file <- snakemake@output[["png"]]

# --- 6. Circos plotting function ---
plot_circos_function <- function(output_type = c("pdf", "png")) {
  output_type <- match.arg(output_type)
  
  if (output_type == "pdf") {
    pdf(pdf_file, width = 10, height = 10)
  } else {
    png(png_file, width = 3000, height = 3000, res = 400)
  }
  
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  circos.clear()
  circos.initializeWithIdeogram(species = "hg19")
  
  # --- 7. Identify links to highlight ---
  sv_data$highlight <- FALSE
  if (!is.na(highlight_chr1) & !is.na(highlight_chr2)) {
    sv_data$highlight <- (
      (sv_data$chrom1 == highlight_chr1 & sv_data$chrom2 == highlight_chr2) |
      (sv_data$chrom1 == highlight_chr2 & sv_data$chrom2 == highlight_chr1)
    )
    n_highlight <- sum(sv_data$highlight)
    if (n_highlight == 0) {
      cat("No links found between", highlight_chr1, "and", highlight_chr2,
          "- drawing all links uniformly.\n")
    } else {
      cat(n_highlight, "links found between", highlight_chr1, "and",
          highlight_chr2, "\n")
    }
  }
  
  # --- 8. Draw all links ---
  for (i in seq_len(nrow(sv_data))) {
    chr1 <- sv_data$chrom1[i]
    p1   <- sv_data$pos1[i]
    chr2 <- sv_data$chrom2[i]
    p2   <- sv_data$pos2[i]
    
    if (!is.na(highlight_chr1) & !is.na(highlight_chr2) & sv_data$highlight[i]) {
      col    <- "#FF5733A0"
      border <- "#FF5733"
      lwd    <- 3
    } else {
      col    <- "#0072B250"
      border <- NA
      lwd    <- 1
    }
    
    try(
      circos.link(chr1, p1, chr2, p2, col = col, border = border, lwd = lwd),
      silent = TRUE
    )
  }
  
  dev.off()
}

# --- 9. Generate outputs ---
cat("Generating PDF...\n")
plot_circos_function("pdf")

cat("Generating PNG...\n")
plot_circos_function("png")

cat("Circos plots generated for", sample_label, "(", sample, ")\n")
