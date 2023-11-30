#!/usr/bin/env Rscript

library(arrow, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)
library(optparse, quietly = T, warn.conflicts = F)
library(MethylationUtils, quietly = T, warn.conflicts = F)

# Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("-s", "--summary"), type = "character",
                     help = "Input sequencing summary file from ONT seqeuncing",
                     metavar = "SEQUENCING_SUMMARY")
parser <- add_option(parser, c("-e", "--experiment"), type = "character",
                     help = "Experiment name",
                     metavar = "EXPERIMENT")

parser <- add_option(parser, c("-r", "--run_id"), type = "character",
                     help = "Run ID",
                     metavar = "RUN_ID")

parser <- add_option(parser, c("-o", "--out"), type="character",
                     help="Hive output location",
                     metavar="OUPUT_DIR", default = "./")

parser <- add_option(parser, c("-t", "--threads"), type="integer",
                     help=paste0("Number of theads, default: available threads [", cpu_count() ,"]"),
                     metavar="THREADS", default = cpu_count())
options <- parse_args(parser)

# Check required arguments
required_options <- c("summary", "experiment", "run_id", "out")
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

# Can/should it grab the experiment/sequencing names from
convert_sequencing_summary_to_parquet_hive <- function(file_path, out_path) {
  # Run ID should be the unique identifier
  # Experiment_id can contain different runs, sample_id could be the same
  # Is it necessary partition on different experiment_id as well?
  read_tsv_arrow(file_path,
                 as_data_frame = FALSE) |> write_dataset(path = out_path,
                                                         format = "parquet",
                                                         partition = c("end_reason"))
}

# Convert modkit bed to parquet hive
convert_sequencing_summary_to_parquet_hive(options$summary, options$out)
