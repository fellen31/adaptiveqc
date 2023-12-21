#!/usr/bin/env Rscript

library(arrow, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)
library(optparse, quietly = T, warn.conflicts = F)

parser <- OptionParser()
parser <- add_option(parser, c("-e", "--experiment"), type = "character",
                     help = "experiment",
                     metavar = "EXPERIMENT")
parser <- add_option(parser, c("-s", "--sample"), type = "character",
                     help = "sample",
                     metavar = "SAMPLE")
parser <- add_option(parser, c("-r", "--run_id"), type = "character",
                     help = "run_id",
                     metavar = "RUN_ID")
parser <- add_option(parser, c("-d", "--dir"), type = "character",
                     help = "Sometimes you have fastq_pass
                             or fastq_fail or fastq_skip",
                     metavar = "DIR")

parser <- add_option(parser, c("-b", "--barcode"), type = "character",
                     help = "barcode",
                     metavar = "BARCODE")

parser <- add_option(parser, c("-i", "--input"), type = "character",
                     help = "input",
                     metavar = "INPUT")

parser <- add_option(parser, c("-o", "--out"), type = "character",
                     help = "Hive output location",
                     metavar = "OUPUT_DIR", default = "./")

parser <- add_option(parser, c("-t", "--threads"), type = "integer",
                     help = paste0("Number of theads, default: available threads [", cpu_count() ,"]"),
                     metavar = "THREADS", default = cpu_count())
options <- parse_args(parser)

# Check required arguments
required_options <- c("experiment", "sample", "run_id", "dir", "barcode", "out")
for (opt in required_options) {
  if (is.null(options[[opt]])) {
    cat(paste0("Required option --", opt, " not found. --help to show more.\n"))
    quit('no', status = 1, runLast = FALSE)
  }
}

# Set number of threads to use from options (default: n available)
set_cpu_count(options$threads)

# Create directory, allow overwrite but show warning
dir.create(file.path(options$out), showWarnings = TRUE)

read_fqcrs_to_arrow <- function(file_path) {
  # Could replace this with fread(cmd = fqcrs ...) to avoid extra files on disk?
  read_tsv_arrow(file_path,
                 col_names = c("read_id", "length", "avg_qual"),
                 as_data_frame = TRUE)
}


file_pathx <- options$input
files <- list.files(file_pathx, recursive = TRUE, pattern = "*.txt.gz", full.names = TRUE)

out_pathx <- options$out

make_hive_from_fqcrs <- function(file_path, out_path) {

  # What to partition on?
  # UUIDs (Universally Unique Identifiers) & file. File needed since multiple files and otherwise might overwrite?

  split_path <- unlist(strsplit(file_path, "/"))
  barcode_dir <- options$barcode
  fastq_dir <- options$dir
  protocol_id_dir <- options$run_id
  sample_id_dir <- options$sample
  experiment_id_dir <- options$sample
  file <- split_path[length(split_path)]

  read_fqcrs_to_arrow(file_path) |>
    mutate(experiment = experiment_id_dir,
           sample = sample_id_dir,
           protocol = protocol_id_dir,
           dir = fastq_dir,
           barcode = barcode_dir,
           file = file) |>
    write_dataset(path = out_path,
                  format = "parquet")
}

lapply(files, make_hive_from_fqcrs, out_path = out_pathx)
