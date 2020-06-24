
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
  walk(commands, system, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Save the peak directory.
  zent_obj <- set_settings(zent_obj, peak_dir = outdir)

  ## Return the zent object.
  return(zent_obj)
}
