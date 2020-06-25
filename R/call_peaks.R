
#' Peak Calling
#'
#' @importFrom purrr pmap
#'
#' @param zent_obj ZentTools object.
#' @param outdir Output directory for peak files.
#' @param genome_size Effective genome size.
#' @param qvalue_cutoff q-value cutoff for peak calling.
#' @param broad_peaks Whether the peaks are broad or not.
#' @param broad_qvalue_cutoff If 'broad_peaks' is set,
#'  This controls the q-value cutoff of the broad peaks.
#'
#' @export

call_peaks <- function(
  zent_obj,
  outdir = getwd(),
  genome_size,
  qvalue_cutoff = 0.05,
  broad_peaks = FALSE,
  broad_qvalue_cutoff = 0.1
) {

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")
  paired_status <- as.logical(pull_setting(zent_obj, "paired"))

  ## Make sure the output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Create the peak calling command.
  commands <- pmap(zent_obj@sample_sheet, function(...) {
    args <- list(...)

    command <- str_c(
      "macs2", "callpeak",
      "-t", args$sample_bams,
      "-c", args$control_bams,
      "-n", args$sample_name,
      "--outdir", outdir,
      "-g", genome_size,
      "-q", qvalue_cutoff,
      sep = " "
    )

    if (broad_peaks) {
      command <- str_c(
        command,
        "--broad",
        "--broad-cutoff", broad_qvalue_cutoff,
        sep = " "
      )
    }

    if (paired_status) {
      command <- str_c(command, "-f", "BAMPE", sep = " ")
    } else {
      command <- str_c(command, "-f", "BAM", sep = " ")
    }

    return(command)
  })

  ## Run the commands.
  print_message("Calling peaks from the aligned reads.")
  walk(commands, system, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Save the peak directory.
  zent_obj <- set_settings(zent_obj, peak_dir = outdir)

  ## Return the zent object.
  return(zent_obj)
}

#' Annotate Peaks
#'
#' @importFrom ChIPseeker annotatePeak
#' @importFrom rtracklayer import
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom purrr iwalk
#'
#' @param zent_obj ZentTools object.
#' @param outdir Directory to output the annotated peaks.
#' @param genome_annotation Genome GTF annotation file.
#' @param promoter_downstream Bases downstream of TSS to define promoter.
#' @param promoter_upstream Bases upstream of TSS to define promoter.
#' @param feature_type Either 'gene' or 'transcript'.
#'
#' @export

annotate_peaks <- function(
  zent_obj,
  outdir = getwd(),
  genome_annotation,
  promoter_downstream = 1000,
  promoter_upstream = 1000,
  feature_type = "gene"
) {

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  ## Make sure the output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE) 

  ## Import the genome annotation file.
  genome_txdb <- makeTxDbFromGFF(genome_annotation)

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

  print_message("Annotating the called peaks.")
  iwalk(peak_files, function(x, y) {
    peaks <- import(x, format = "BED", extraCols = narrowpeak_cols)

    annotated_peaks <- annotatePeak(
      peaks,
      tssRegion = c(-promoter_upstream, promoter_downstream),
      TxDb = genome_txdb,
      level = feature_type
    )
    annotated_peaks <- as.data.table(annotated_peaks)

    fwrite(
      annotated_peaks, str_c(outdir, y, "_peaks_annotated.tsv"),
      sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
    )
  })

  ## Return the zent object.
  return(zent_obj)

}
