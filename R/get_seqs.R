
#' Get Genomic Sequences
#'
#' @importFrom Rsamtools FaFile getSeq
#' @importFrom rtracklayer import
#'
#' @param zent_obj ZentTools object.
#' @param outdir Output directory for fasta file.
#' @param fixed_width Make the peaks a fixed width from the center of the peak.
#'
#' @export

get_seqs <- function(
  zent_obj,
  outdir = getwd(),
  fixed_width = NA
) {

  ## Validate inputs.
  if (!str_ends(outdir, "/")) outdir <- str_c(outdir, "/")

  ## Import the peak files.
  peak_files <- str_c(
    pull_setting(zent_obj, "peak_dir"),
    zent_obj@sample_sheet[["sample_name"]],
    "_peaks.narrowPeak"
  )
  names(peak_files) <- zent_obj@sample_sheet[["sample_name"]]

  narrowpeak_cols  <-  c(
    signal.value = "numeric",
    p.value.negLog10 = "numeric",
    q.value.negLog10 = "numeric",
    peak = "integer"
  )

  print_message("Importing the peak files.")
  peaks <- map(
    peak_files, import,
    format = "BED", extraCols = narrowpeak_cols
  )

  ## If specified, make the peaks a fixed width from the center.
  if (!is.na(fixed_width)) {
    NULL
  }

  ## Create reference to the genome assembly.
  assembly <- pull_setting(zent_obj, "genome_assembly")
  assembly <- FaFile(assembly)

}
