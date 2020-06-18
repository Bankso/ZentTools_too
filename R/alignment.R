
#' STAR Genome Generation
#'
#' @param zent_obj Zent object.
#' @param outdir Output directory for STAR genome index (--genomeDir).
#' @param genome_assembly Path to genome fasta file (--genomeFastaFiles).
#' @param genome_annotation Path to genome gtf file (--sjdbGTFfile).
#' @param sa_bases Number of suffix array bases (--genomeSAindexNbases)
#'
#' @export

star_genome <- function(
  zent_obj,
  outdir = getwd(),
  genome_assembly,
  genome_annotation,
  sa_bases = NA
) {

  ## Make sure output directory exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Prepare STAR genome index command.
  command <- str_c(
    "STAR",
    "--runThreadN", zent_obj@settings[parameter == "ncores", value],
    "--runMode", "genomeGenerate",
    "--genomeDir", outdir,
    "--genomeFastaFiles", genome_assembly,
    "--sjdbGTFfile", genome_annotation,
    sep = " "
  )

  if (!is.na(sa_bases)) {
    command <- str_c(
      command,
      "--genomeSAindexNbases",
      sa_bases,
      sep = " "
    )
  }

  ## Run the command.
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Store the genome directory.
  new_settings <- list(data.table(parameter = "genome_dir", value = outdir))

  settings <- copy(zent_obj@settings)
  settings <- rbindlist(list(settings), new_settings)

  zent_obj@settings <- settings

  ## Return the zent object.
  return(zent_obj)

}

#' STAR Alignment
#'
#' @importFrom purrr iwalk map
#'
#' @param zent_obj Zent object.
#' @param outdir Output directory for aligned reads.
#'   For use with '--outFileNamePrefix'.
#'
#' @export

star_align <- function(
  zent_obj,
  outdir = getwd()
) {

  ## Create output directory if it exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get sample names.
  paired_status <- zent_obj@settings[parameter == "paired", value]

  if (paired_status) {
    samples <- split(zent_obj@sample_sheet, by = "sample_name", keep.by = FALSE)
    samples <- map(samples, function(x) {
      x <- str_c(as.character(x), collapse = " ")
      return(x)
    })
  } else {
    samples <- split(zent_obj@sample_sheet, by = "sample_name", keep.by = FALSE)
    samples <- map(samples, function(x) {
      x[, file_2 := NULL]
      x <- as.character(x)
      return(x)
    })
  }

  ## prepare STAR alignment command.
  command <- str_c(
    "STAR",
    "--runThreadN", zent_obj@settings[parameter == "ncols", value],
    "--genomeDir", zent_obj@settings[parameter == "genome_dir", value],
    "--outSAMtype", "BAM SortedByCoordinate"
  )

  ## Run the command.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  iwalk(samples, function(x, y) {
    command <- str_c(
      command,
      "--outFileNamePrefix", str_c(outdir, str_c(y, "_")),
      "--readFilesIn", x,
      sep = " "
    )
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  })

  ## Add the outdir directory for alignment to the settings.
  new_setting <- list(data.table(parameter = "alignment_dir", value = outdir))
  settings <- list(copy(zent_obj@settings))

  settings <- rbindlist(settings, new_setting)
  zent_obj@settings <- settings

  ## Return the zent object.
  return(zent_obj)

}
