
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
  new_settings <- data.table(
    parameter = c("genome_dir", "genome_annotation"),
    value = c(outdir, genome_annotation)
  )
  settings <- copy(zent_obj@settings)
  settings <- rbindlist(list(settings, new_settings))

  zent_obj@settings <- settings

  ## Return the zent object.
  return(zent_obj)

}

#' STAR Alignment
#'
#' @importFrom purrr walk map imap
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
  paired_status <- as.logical(zent_obj@settings[parameter == "paired", value])

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
    "--runThreadN", zent_obj@settings[parameter == "ncores", value],
    "--genomeDir", zent_obj@settings[parameter == "genome_dir", value],
    "--outSAMtype", "BAM SortedByCoordinate",
    sep = " "
  )

  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  command <- imap(samples, function(x, y) {
    command <- str_c(
      command,
      "--outFileNamePrefix", str_c(outdir, str_c(y, "_")),
      "--readFilesIn", x,
      sep = " "
    )
    return(command)
  })

  ## Run the command.
  walk(command, ~system(., ignore.stdout = TRUE, ignore.stderr = TRUE))

  ## Index the bam files.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  bam_files <- str_c(
    outdir, 
    zent_obj@sample_sheet[["sample_name"]],
    "_Aligned.sortedByCoord.out.bam"
  )

  zent_obj@sample_sheet[["bam_files"]] <- bam_files

  walk(bam_files, function(x) {
    command <- str_c("samtools", "index", x, sep = " ")
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  })

  ## Add the outdir directory for alignment to the settings.
  new_setting <- data.table(
    parameter = "alignment_dir",
    value = outdir
  )
  settings <- copy(zent_obj@settings)

  settings <- rbindlist(list(settings, new_setting))
  zent_obj@settings <- settings

  ## Return the zent object.
  return(zent_obj)

}
