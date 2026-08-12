#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(circlize))

usage <- function() {
  cat(paste0(
    "Usage:\n",
    "  Rscript plot_mitogenome_circos.R --repeats repeats.tsv --genes genes.tsv ",
    "--gc gc.tsv --links links.tsv --genome-length 376041 --output figure.pdf ",
    "[--chromosome-id Mito] [--start-degree 83] [--gap-degree 14]\n\n",
    "Input tables are tab-separated and supplied through command-line arguments.\n"
  ))
}

parse_args <- function(values) {
  result <- list(chromosome_id = "Mito", start_degree = 83, gap_degree = 14)
  i <- 1
  while (i <= length(values)) {
    key <- values[[i]]
    if (key %in% c("-h", "--help")) {
      usage()
      quit(status = 0)
    }
    if (!startsWith(key, "--") || i == length(values)) {
      stop(paste("Invalid argument:", key), call. = FALSE)
    }
    name <- gsub("-", "_", substring(key, 3))
    result[[name]] <- values[[i + 1]]
    i <- i + 2
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("repeats", "genes", "gc", "links", "genome_length", "output")
missing <- required[!required %in% names(args)]
if (length(missing) > 0) {
  usage()
  stop(paste("Missing required arguments:", paste(missing, collapse = ", ")), call. = FALSE)
}

genome_length <- as.numeric(args$genome_length)
start_degree <- as.numeric(args$start_degree)
gap_degree <- as.numeric(args$gap_degree)
chromosome_id <- args$chromosome_id
if (!is.finite(genome_length) || genome_length <= 0) {
  stop("--genome-length must be a positive number", call. = FALSE)
}

read_table <- function(path, columns) {
  data <- read.table(
    path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
    comment.char = "", quote = "", fill = FALSE
  )
  if (ncol(data) != length(columns)) {
    stop(sprintf("%s must contain exactly %d tab-separated columns", path, length(columns)), call. = FALSE)
  }
  colnames(data) <- columns
  data
}

repair_coordinates <- function(data, start_column, end_column) {
  starts <- data[[start_column]]
  data[[start_column]] <- pmin(starts, data[[end_column]])
  data[[end_column]] <- pmax(starts, data[[end_column]])
  data
}

normalize_color <- function(value, soft_yellow = "#F9E79F") {
  if (tolower(value) == "yellow" || toupper(value) == "#FFFF00" || value == "255,255,0") {
    return(soft_yellow)
  }
  if (grepl("^[0-9]+,[0-9]+,[0-9]+$", value)) {
    rgb_values <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
    if (length(rgb_values) != 3 || any(rgb_values < 0 | rgb_values > 255)) {
      stop(paste("Invalid RGB color:", value), call. = FALSE)
    }
    return(rgb(rgb_values[1], rgb_values[2], rgb_values[3], maxColorValue = 255))
  }
  value
}

repeats <- read_table(args$repeats, c("chr", "start", "end", "value", "color"))
genes <- read_table(args$genes, c("chr", "start", "end", "color"))
gc_data <- read_table(args$gc, c("chr", "start", "end", "value"))
links <- read_table(args$links, c("q_chr", "q_start", "q_end", "s_chr", "s_start", "s_end"))

repeats$chr <- chromosome_id
genes$chr <- chromosome_id
gc_data$chr <- chromosome_id
links$q_chr <- chromosome_id
links$s_chr <- chromosome_id
repeats <- repair_coordinates(repeats, "start", "end")
genes <- repair_coordinates(genes, "start", "end")
links <- repair_coordinates(links, "q_start", "q_end")
links <- repair_coordinates(links, "s_start", "s_end")
links <- links[!(links$q_start == links$s_start & links$q_end == links$s_end), , drop = FALSE]
repeats$color <- vapply(repeats$color, normalize_color, character(1))
genes$color <- vapply(genes$color, normalize_color, character(1))

link_color <- function(start, end) {
  length <- abs(end - start)
  color <- if (length >= 1000) {
    "magenta"
  } else if (length >= 300) {
    "green"
  } else if (length >= 200) {
    "blue"
  } else if (length >= 100) {
    "orange"
  } else if (length >= 50) {
    "#F9E79F"
  } else {
    "grey"
  }
  add_transparency(color, 0.5)
}
link_colors <- mapply(link_color, links$q_start, links$q_end)

open_device <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension == "pdf") {
    pdf(path, width = 7.2, height = 7.2, useDingbats = FALSE)
  } else if (extension == "svg") {
    svg(path, width = 7.2, height = 7.2)
  } else if (extension == "png") {
    png(path, width = 2400, height = 2400, res = 300)
  } else {
    stop("Output extension must be .pdf, .svg, or .png", call. = FALSE)
  }
}

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
open_device(args$output)
on.exit({ circos.clear(); dev.off() }, add = TRUE)

circos.clear()
circos.par(
  track.margin = c(0.01, 0.01), start.degree = start_degree,
  gap.after = gap_degree, cell.padding = c(0, 0, 0, 0)
)
reference <- data.frame(chr = chromosome_id, start = 0, end = genome_length)
circos.genomicInitialize(reference, plotType = NULL)

circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.col = "#F28E2B", bg.border = "black",
  panel.fun = function(x, y) {
    ticks <- seq(0, floor(genome_length / 10000) * 10000, by = 10000)
    labels <- vapply(ticks, function(value) {
      if (value %% 50000 == 0) {
        if (value == 0) "0" else paste0(value / 1000, "K")
      } else ""
    }, character(1))
    tick_lengths <- ifelse(ticks %% 50000 == 0, 0.4, 0.2)
    circos.axis(
      h = "top", major.at = ticks, labels = labels, minor.ticks = 0,
      labels.cex = 0.5, labels.facing = "downward", labels.niceFacing = TRUE,
      major.tick.length = tick_lengths
    )
  }
)

circos.genomicTrack(
  gc_data, track.height = 0.06, bg.col = "white", bg.border = "black",
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, type = "l", col = "#7F7F7F", lwd = 0.6, ...)
  }
)

circos.genomicTrack(
  genes, stack = TRUE, track.height = 0.05, bg.border = NA,
  panel.fun = function(region, value, ...) {
    level <- getI(...)
    circos.genomicRect(
      region, value, ybottom = level - 0.45, ytop = level + 0.45,
      col = value$color, border = NA, ...
    )
  }
)

circos.genomicTrack(
  repeats, stack = TRUE, track.height = 0.05, bg.border = NA,
  panel.fun = function(region, value, ...) {
    level <- getI(...)
    circos.genomicRect(
      region, value, ybottom = level - 0.45, ytop = level + 0.45,
      col = value$color, border = NA, ...
    )
  }
)

if (nrow(links) > 0) {
  circos.genomicLink(
    links[, 1:3], links[, 4:6], col = link_colors, border = NA, rou = 0.65
  )
}
legend(
  "center", legend = c(">=1000", "300-999", "200-299", "100-199", "50-99", "<50"),
  fill = c("magenta", "green", "blue", "orange", "#F9E79F", "grey"),
  border = NA, bty = "n", cex = 0.6, title = "Repeat (bp)"
)

message("Wrote ", args$output)
