
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

  ## Get file names.
  paired_status <- as.logical(zent_obj@settings[parameter == "paired", value])

  samples <- str_c(zent_obj@sample_sheet[["file_1"]], collapse = " ")

  if (paired_status) {
    samples <- str_c(
      samples,
      str_c(zent_obj@sample_sheet[["file_2"]], collapse = " "),
      sep = " "
    )
  }

  ## Prepare the fastqc command.
  command <- str_c(
    "fastqc",
    "-o", outdir,
    "-t", zent_obj@settings[parameter == "ncores", value],
    samples,
    sep = " "
  )

  ## Run the fastqc command.
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Return the zent object.
  return(zent_obj)

}
