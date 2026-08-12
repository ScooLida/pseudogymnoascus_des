#!/usr/bin/env Rscript

# Create one depth-of-coverage plot per target gene.
# The x-axis shows nucleotide position, the y-axis shows read depth,
# and each line represents one sample that passed the filtering step.

# Configuration
OUT_DIR <- "./genes_report"
TARGET_SAMPLES_FILE <- file.path(OUT_DIR, "target_samples.txt")
PLOT_DIR <- file.path(OUT_DIR, "depth_plots")
TABLE_DIR <- file.path(OUT_DIR, "depth_tables")

GENES <- c("Pd_18S", "Pd_ITS", "Pd_28S", "Pd_MCM7", "Pd_TEF1alpha", "Pd_RPB2")
PLOT_WIDTH <- 12
PLOT_HEIGHT <- 7
PLOT_DPI <- 300

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required. Install it with install.packages('ggplot2').")
}

if (!file.exists(TARGET_SAMPLES_FILE)) {
  stop("Selected-sample list not found: ", TARGET_SAMPLES_FILE)
}

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

target_samples <- readLines(TARGET_SAMPLES_FILE, warn = FALSE)
target_samples <- trimws(target_samples)
target_samples <- unique(target_samples[nzchar(target_samples)])

if (length(target_samples) == 0) {
  stop("The selected-sample list is empty.")
}

read_gene_table <- function(gene) {
  input_file <- file.path(
    OUT_DIR,
    paste0(gene, "_fil_cov_", length(target_samples), ".tsv")
  )

  if (!file.exists(input_file)) {
    stop("Coverage table not found for ", gene, ": ", input_file)
  }

  data <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  required_columns <- c("Sample", "Position", "Depth")
  if (!all(required_columns %in% names(data))) {
    stop("Unexpected columns in ", input_file, ". Expected: Sample, Position, Depth")
  }

  data <- data[data$Sample %in% target_samples, required_columns]
  data$Position <- as.numeric(data$Position)
  data$Depth <- as.numeric(data$Depth)
  data <- data[!is.na(data$Position) & !is.na(data$Depth), ]
  data$Sample <- factor(data$Sample, levels = target_samples)
  data[order(data$Sample, data$Position), ]
}

for (gene in GENES) {
  depth_table <- read_gene_table(gene)

  if (nrow(depth_table) == 0) {
    warning("No depth data for ", gene, "; skipping.")
    next
  }

  write.table(
    depth_table,
    file = file.path(TABLE_DIR, paste0(gene, "_depth.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  plot <- ggplot2::ggplot(
    depth_table,
    ggplot2::aes(x = Position, y = Depth, color = Sample, group = Sample)
  ) +
    ggplot2::geom_line(linewidth = 0.4, na.rm = TRUE) +
    ggplot2::labs(
      title = paste("Read depth across", gene),
      subtitle = paste(length(target_samples), "selected samples"),
      x = "Nucleotide position",
      y = "Read depth",
      color = "Sample"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  ggplot2::ggsave(
    filename = file.path(PLOT_DIR, paste0(gene, "_depth.png")),
    plot = plot,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT,
    dpi = PLOT_DPI
  )
}

message("Depth plots saved to: ", PLOT_DIR)
message("Filtered depth tables saved to: ", TABLE_DIR)
