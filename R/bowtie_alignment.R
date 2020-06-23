
#' Bowtie2 Index Generation
#'
#' @param zent_obj Zent object.
#' @param outdir Output directory for STAR genome index.
#' @param genome_assembly Path to genome fasta file.
#' @param index_name The naming structure given to the index.
#'
#' @export

bowtie2_index <- function(
  zent_obj,
  genome_assembly,
  outdir = getwd(),
  index_name = "bowtie2_index"
) {

  ## Make sure output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Prepare Bowtie2 index command.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  command <- str_c(
    "bowtie2-build",
    "-f", genome_assembly,
    "--threads", zent_obj@settings[paremter == "ncores", value],
    str_c(outdir, index_name),
    sep = " "
  )

  ## Run the command.
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Store the genome directory.
  new_settings <- data.table(
    parameter = c("genome_dir", "genome_assembly"),
    value = c(str_c(outdir, index_name), genome_assembly)
  )
  settings <- copy(zent_obj@settings)
  settings <- rbindlist(list(settings, new_settings))

  zent_obj@settings <- settings

  ## Return the zent object.
  return(zent_obj)

}

#' Bowtie2 Alignment
#'
#' @importFrom purrr pmap
#'
#' @param zent_obj Zent object.
#' @param outdir Output directory for aligned reads.
#' @param alignment_mode Either 'end-to-end' or 'local'.
#' @param min_fragment Minimum fragment length (paired end).
#' @param max_fragment Maximum fragment length (paired end).
#'
#' @export

bowtie2_align <- function(
  zent_obj,
  outdir = getwd(),
  alignment_mode = "end-to-end",
  min_fragment = NA,
  max_fragment = NA
) {

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  ## Create output directory if it exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Prepare bowtie2 alignment command.
  paired_status <- as.logical(zent_obj@settings[parameter == "paired", value])

  commands <- pmap(zent_obj@sample_sheet, function(...) {
    args <- list(...)

    command <- str_c(
      "bowtie2",
      "-x", zent_obj@settings[paramemter == "genome_dir", value],
      "-S", str_c(outdir, args$sample_name, ".sam"),
      "--phred33",
      "--no-unal",
      "-p", zent_obj@settings[parameter == "ncores", value],
      sep = " "
    )

    if (paired_status) {
      command <- str_c(
        command,
        "-1", args$file_1,
        "-2", args$file_2,
        "--no-mixed",
        "--no-discordant",
        sep = " "
      )

      if (!is.na(min_fragment)) {
        command <- str_c(command, "-I", min_fragment, sep = " ")
      }
      if (!is.na(max_fragment)) {
        command <- str_c(command, "-X", max_fragment, sep = " ")
      }
    } else {
      command <- str_c(command, "-U", args$file_1, sep = " ")
    }

    return(command)
  })

  ## Run commands.
  walk(commands, function(x) {
    system(x, ignore.stdout = TRUE, ignore.stderr = TRUE)
  })

  ## Make coordinate sorted and indexed bams.

  ## Return the zent object.
  return(zent_obj)

}
