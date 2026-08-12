#!/usr/bin/env Rscript

# Create depth plots from per-gene tables in Sample/Position/Depth format.
# Each plot contains individual sample lines, the interquartile range,
# and the median depth across the selected samples.

INPUT_DIR <- "/home/lidacool/PycharmProjects/grib"
TARGET_SAMPLES_FILE <- file.path(INPUT_DIR, "target_samples.txt")
INPUT_PATTERN <- "^Pd_(18S|ITS|28S|MCM7|TEF1alpha|RPB2)_fil_cov_71\\.tsv$"
OUTPUT_DIR <- file.path(INPUT_DIR, "depth_plots")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 first: install.packages('ggplot2')")
}
if (!file.exists(TARGET_SAMPLES_FILE)) {
  stop("Target sample list not found: ", TARGET_SAMPLES_FILE)
}

target_samples <- unique(trimws(readLines(TARGET_SAMPLES_FILE, warn = FALSE)))
target_samples <- target_samples[nzchar(target_samples)]
if (length(target_samples) == 0) stop("The target sample list is empty.")

input_files <- list.files(INPUT_DIR, pattern = INPUT_PATTERN, full.names = TRUE)
if (length(input_files) != 6) {
  stop("Expected 6 gene tables, found ", length(input_files), ": ", INPUT_PATTERN)
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (input_file in sort(input_files)) {
  data <- read.delim(
    input_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  required_columns <- c("Sample", "Position", "Depth")
  if (!all(required_columns %in% names(data))) {
    stop("Unexpected columns in: ", input_file)
  }

  data <- data[data$Sample %in% target_samples, required_columns]
  data$Position <- as.numeric(data$Position)
  data$Depth <- as.numeric(data$Depth)
  data <- data[!is.na(data$Position) & !is.na(data$Depth), ]

  if (nrow(data) == 0) {
    warning("No target samples found in: ", input_file)
    next
  }

  data$Sample <- factor(data$Sample, levels = target_samples)
  data <- data[order(data$Sample, data$Position), ]

  quantiles <- by(data$Depth, data$Position, function(values) {
    c(
      Q25 = unname(quantile(values, 0.25, na.rm = TRUE)),
      Median = median(values, na.rm = TRUE),
      Q75 = unname(quantile(values, 0.75, na.rm = TRUE))
    )
  })
  summary_data <- data.frame(
    Position = as.numeric(names(quantiles)),
    do.call(rbind, quantiles),
    row.names = NULL
  )

  gene <- sub("_fil_cov_71\\.tsv$", "", basename(input_file))
  output_file <- file.path(OUTPUT_DIR, paste0(gene, "_depth.png"))

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = Position, y = Depth, group = Sample)
  ) +
    ggplot2::geom_ribbon(
      data = summary_data,
      ggplot2::aes(x = Position, ymin = Q25, ymax = Q75),
      inherit.aes = FALSE,
      fill = "steelblue",
      alpha = 0.25
    ) +
    ggplot2::geom_line(color = "steelblue3", linewidth = 0.45, alpha = 0.4) +
    ggplot2::geom_line(
      data = summary_data,
      ggplot2::aes(x = Position, y = Median),
      inherit.aes = FALSE,
      color = "navy",
      linewidth = 1.3
    ) +
    ggplot2::labs(
      title = paste(gene, "coverage depth across selected samples"),
      subtitle = paste(length(unique(data$Sample)), "samples"),
      x = "Nucleotide position along the gene",
      y = "Sequencing depth (number of reads)"
    ) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      axis.title = ggplot2::element_text(size = 17),
      axis.text = ggplot2::element_text(size = 14),
      legend.position = "none"
    )

  ggplot2::ggsave(output_file, plot, width = 16, height = 9, dpi = 300)
  message("Saved: ", output_file)
}
