
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
  analysis_type <- pull_setting(zent_obj, "analysis_type")
  paired_status <- as.logical(pull_setting(zent_obj, "paired"))

  ## Get file names.
  samples <- zent_obj@sample_sheet[["file_1"]]

  if (paired_status) {
    samples <- c(samples, zent_obj@sample_sheet[["file_2"]])
  }

  if (
    analysis_type %in% c("ChEC-seq", "ChIP-seq") &&
    any(!is.na(unique(zent_obj@sample_sheet[["control_file_1"]])))
  ) {
    controls_1 <- unique(zent_obj@sample_sheet[["control_file_1"]])
    controls_1 <- controls_1[!is.na(controls_1)]

    samples <- c(samples, controls_1)

    if (paired_status) {
      controls_2 <- unique(zent_obj@sample_sheet[["control_file_2"]])
      controls_2 <- controls_2[!is.na(controls_2)]

      samples <- c(samples, controls_2)
    }
  }

  ## Prepare the fastqc command.
  command <- str_c(
    "fastqc",
    "-o", outdir,
    "-t", pull_setting(zent_obj, "ncores"),
    str_c(samples, collapse = " "),
    sep = " "
  )

  ## Run the fastqc command.
  print_message("Running FASTQ quality control.")
  system(command)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Return the zent object.
  return(zent_obj)

}
