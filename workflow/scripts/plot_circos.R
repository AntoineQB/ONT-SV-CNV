# =============================================================================
#   Circos Plot
#
#   Generates a circos plot of structural variant breakpoints
#   with optional chromosome-pair highlighting.
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

# --- 2. Output file ---
png_file <- snakemake@output[["png"]]

# --- 3. Read annotated CSV ---
input_file <- snakemake@input[["csv"]]
cat("Reading:", input_file, "\n")
sv_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# --- 4. Handle empty input gracefully ---
if (nrow(sv_data) == 0) {
  cat("No variants in input — generating empty circos plot.\n")
  png(png_file, width = 3000, height = 3000, res = 400)
  par(mar = c(0.5, 0.5, 2, 0.5))
  circos.clear()
  circos.initializeWithIdeogram(species = "hg19")
  title(paste0(sample_label, " — No structural variants detected"), cex.main = 1.2)
  dev.off()
  cat("Empty circos plot generated for", sample_label, "(", sample, ")\n")
  quit(save = "no", status = 0)
}

# --- 5. Harmonize column names ---
if (!all(c("chrom1", "pos1", "chrom2", "pos2") %in% names(sv_data))) {
  sv_data$chrom1 <- sv_data$Chr1
  sv_data$pos1   <- as.numeric(sv_data$Pos1)
  sv_data$chrom2 <- sv_data$Chr2
  sv_data$pos2   <- as.numeric(sv_data$Pos2)
}

# --- 6. Ensure "chr" prefix on chromosome names ---
if (!any(grepl("^chr", sv_data$chrom1))) {
  sv_data$chrom1 <- paste0("chr", sv_data$chrom1)
  sv_data$chrom2 <- paste0("chr", sv_data$chrom2)
}

if (!is.na(highlight_chr1)) {
  highlight_chr1 <- ifelse(grepl("^chr", highlight_chr1), highlight_chr1, paste0("chr", highlight_chr1))
}
if (!is.na(highlight_chr2)) {
  highlight_chr2 <- ifelse(grepl("^chr", highlight_chr2), highlight_chr2, paste0("chr", highlight_chr2))
}

# --- 7. Plotting function ---
plot_circos_png <- function() {
  png(png_file, width = 3000, height = 3000, res = 400)
  
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  circos.clear()
  circos.initializeWithIdeogram(species = "hg19")
  
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
      cat(n_highlight, "links found between", highlight_chr1, "and", highlight_chr2, "\n")
    }
  }
  
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

# --- 8. Generate output ---
cat("Generating PNG...\n")
plot_circos_png()
cat("Circos plot generated for", sample_label, "(", sample, ")\n")