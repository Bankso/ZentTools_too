
#' Quality Control of Fastq Files
#'
#' @param zent_obj Zent object.
#' @param outdir Output directory for reports (-o).
#'
#' @export

fastqc <- function(
  zent_obj,
  outdir = getwd()
) {

  ## Ensure the output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  analysis_type <- zent_obj@settings[parameter == "analysis_type", value]
  paired_status <- as.logical(zent_obj@settings[parameter == "paired", value])

  ## Get file names.
  samples <- zent_obj@sample_sheet[["file_1"]]

  if (paired_status) {
    samples <- c(samples, zent_obj@sample_sheet[["file_2"]])
  }

  if (analysis_type %in% c("ChEC-seq", "ChIP-seq")) {
    samples <- c(samples, unique(zent_obj@sample_sheet[["control_file_1"]]))

    if (paired_status) {
      samples <- c(samples, unique(zent_obj@sample_sheet[["control_file_2"]]))
    }
  }

  ## Prepare the fastqc command.
  command <- str_c(
    "fastqc",
    "-o", outdir,
    "-t", zent_obj@settings[parameter == "ncores", value],
    str_c(samples, collapse = " "),
    sep = " "
  )

  ## Run the fastqc command.
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Return the zent object.
  return(zent_obj)

}
