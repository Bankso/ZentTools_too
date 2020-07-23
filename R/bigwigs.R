
#' Generate Bigwigs
#'
#' @param zent_obj ZentTools object.
#' @param outdir Output directory.
#' @param bin_size Bin size for coverage summary.
#' @param normalize_using Either 'CPM' or 'RPGC'.
#' @param genome_size Effective genome size.
#' @param skip_non_covered Should regions without coverage be skipped.
#' @param min_fragment Minimum fragment length.
#' @param max_fragment Maximum fragment length.
#' @param extend_reads Distance to extend single-end reads.
#'   Set to NA to not extend reads.
#' @param scale_factors Takes a named vector, with the name being the sample name,
#'   and the value being the scale factor.
#'   If set will override 'normalzie_using' option.
#'
#' @export

make_bigwigs <- function(
  zent_obj,
  outdir = getwd(),
  bin_size = 10,
  normalize_using = NA,
  genome_size = NA,
  skip_non_covered = TRUE,
  min_fragment = NA,
  max_fragment = NA,
  extend_reads = 200,
  scale_factors = NA
) {

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")
  paired_status <- as.logical(pull_setting(zent_obj, "paired"))

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get bams.
  samples <- split(
    zent_obj@sample_sheet[, .(sample_name, sample_bams)],
    by = "sample_name",
    keep.by = FALSE
  )
  samples <- map(samples, as.character)

  if(any(!is.na(zent_obj@sample_sheet[["control_bams"]]))) {
    controls <- split(
      unique(zent_obj@sample_sheet[
        !is.na(control_bams),
        .(control_name, control_bams)
      ]),
      by = "control_name",
      keep.by = FALSE
    )
    controls <- map(controls, as.character)
    samples <- c(samples, controls)
  }

  ## Prepare command.
  commands <- imap(samples, function(x, y) {
    command <- str_c(
      "bamCoverage",
      "-b", x,
      "-o", str_c(outdir, y, ".bigwig"),
      "-of", "bigwig",
      "-bs", bin_size,
      "-p", pull_setting(zent_obj, "ncores"),
      sep = " "
    )

    if (all(is.na(scale_factors)) && !is.na(normalize_using)) {
      command <- str_c(command, "--normalizeUsing", normalize_using, sep = " ")
    }

    if (all(!is.na(scale_factors))) {
      command <- str_c(command, str_c("--scaleFactor", scale_factors[y], sep = " "), sep = " ")
    }

    if (!is.na(genome_size)) {
      command <- str_c(command, "--effectiveGenomeSize", genome_size, sep = " ")
    }

    if (skip_non_covered) {
      command <- str_c(command, "--skipNonCoveredRegions", sep = " ")
    }

    if (!is.na(min_fragment)) {
      command <- str_c(command, "--minFragmentLength", min_fragment, sep = " ")
    }
    if (!is.na(max_fragment)) {
      command <- str_c(command, "--maxFragmentLength", max_fragment, sep = " ")
    }

    if (!is.na(extend_reads) && paired_status) {
      command <- str_c(command, "-e", sep = " ")
    } else if (!is.na(extend_reads) && !paired_status) {
      command <- str_c(command, "-e", extend_reads, sep = " ")
    }

    return(command)
  })

  ## Run commands.
  print_message("Creating the BIGWIG coverage tracks.")
  walk(commands, system)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Return zent tools object.
  return(zent_obj)

}
