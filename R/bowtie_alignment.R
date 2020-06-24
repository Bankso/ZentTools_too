
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

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  ## Make sure output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Prepare Bowtie2 index command.
  command <- str_c(
    "bowtie2-build",
    "-f", genome_assembly,
    "--threads", zent_obj@settings[parameter == "ncores", value],
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
#' @importFrom purrr walk imap
#'
#' @param zent_obj Zent object.
#' @param outdir Output directory for aligned reads.
#' @param alignment_mode Either 'end-to-end' or 'local'.
#' @param min_fragment Minimum fragment length (paired end).
#' @param max_fragment Maximum fragment length (paired end).
#' @param max_memory Maximum memory per thread for samtools.
#'
#' @export

bowtie2_align <- function(
  zent_obj,
  outdir = getwd(),
  alignment_mode = "end-to-end",
  min_fragment = NA,
  max_fragment = NA,
  max_memory = "1G"
) {

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")
  paired_status <- as.logical(zent_obj@settings[parameter == "paired", value])

  ## Create output directory if it exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  if (paired_status) {
    samples <- split(
      zent_obj@sample_sheet[, .(sample_name, file_1, file_2)],
      by = "sample_name",
      keep.by = FALSE
    )
    samples <- map(samples, as.character)

    controls <- split(
      unique(zent_obj@sample_sheet[, .(control_name, control_file_1, control_file_2)]),
      by = "control_name",
      keep.by = FALSE
    )
    controls <- map(controls, as.character)
  } else {
    samples <- split(
      zent_obj@sample_sheet[, .(sample_name, file_1)],
      by = "sample_name",
      keep.by = FALSE
    )
    samples <- map(samples, as.character)

    controls <- split(
      unique(zent_obj@sample_sheet[, .(control_name, control_file_1)]),
      by = "control_name",
      keep.by = FALSE
    )
    controls <- map(controls, as.character)
  }

  samples <- c(samples, controls)

  ## Prepare bowtie2 alignment command.
  commands <- imap(samples, function(x, y) {
    command <- str_c(
      "bowtie2",
      "-x", zent_obj@settings[parameter == "genome_dir", value],
      "-S", str_c(outdir, y, ".sam"),
      "--phred33",
      "--no-unal",
      "-p", zent_obj@settings[parameter == "ncores", value],
      sep = " "
    )

    if (paired_status) {
      command <- str_c(
        command,
        "--no-mixed",
        "--no-discordant",
        "-1", x[1],
        "-2", x[2],
        sep = " "
      )

      if (!is.na(min_fragment)) {
        command <- str_c(command, "-I", min_fragment, sep = " ")
      }
      if (!is.na(max_fragment)) {
        command <- str_c(command, "-X", max_fragment, sep = " ")
      }
    } else {
      command <- str_c(command, "-U", x, sep = " ")
    }

    return(command)
  })

  ## Run the commands.
  walk(commands, system, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Make coordinate sorted and indexed bams.
  walk(names(samples), function(x) {
    command <- str_c(
      "samtools", "sort",
      "-m", max_memory,
      "-@", zent_obj@settings[parameter == "ncores", value],
      "-o", str_c(outdir, str_c(x, ".bam")),
      "-O", "BAM",
      str_c(outdir, str_c(x, ".sam")),
      sep = " "
    )
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

    command <- str_c(
      "samtools", "index",
      str_c(outdir, str_c(x, ".bam")),
      sep = " "
    )
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  })

  ## Add settings to zent object.
  new_settings <- data.table(alignment_dir = outdir)
  settings <- copy(zent_obj@settings)
  settings <- rbindlist(list(settings, new_settings))
  zent_obj@settings <- settings

  ## Add bam files to sample_sheet.
  sample_sheet <- copy(zent_obj@sample_sheet)
  sample_sheet[
    c("sample_bams", "control_bams") := list(
      str_c(outdir, sample_name, ".bam"),
      str_c(outdir, control_name, ".bam")
    )
  ]
  zent_obj@sample_sheet <- sample_sheet

  ## Return the zent object.
  return(zent_obj)

}
