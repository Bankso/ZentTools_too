#' STAR Genome Generation
#'
#' @param zent_obj Zent object.
#' @param outdir Output directory for STAR genome index.
#' @param genome_assembly Path to genome fasta file.
#' @param genome_annotation Path to genome gtf file.
#' @param sa_bases Number of suffix array bases (--genomeSAindexNbases)
#'
#' @export

star_index <- function(
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
    "--runThreadN", pull_setting(zent_obj, "ncores"),
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
  print_message("Generating the STAR genome index.")
  system(command)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Store the genome directory.
  zent_obj <- set_settings(
    zent_obj,
    genome_dir = outdir,
    genome_annotation = genome_annotation
  )

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

  ## Input check and getting settings.
  paired_status <- as.logical(pull_setting(zent_obj, "paired"))
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  ## Create output directory if it exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get sample names.
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
    "--runThreadN", pull_setting(zent_obj, "ncores"),
    "--genomeDir", pull_setting(zent_obj, "genome_dir"),
    "--outSAMtype", "BAM SortedByCoordinate",
    sep = " "
  )

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
  print_message("Aligning the FASTQ files using STAR.")
  walk(command, system)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Index the bam files.
  zent_obj <- add_bams(zent_obj, alignment_dir = outdir)

  print_message("Indexing the aligned BAM files.")
  walk(zent_obj@sample_sheet[["bam_files"]], function(x) {
    command <- str_c("samtools", "index", x, sep = " ")
    system(command)#, ignore.stdout = TRUE, ignore.stderr = TRUE)
  })

  ## Add the outdir directory for alignment to the settings.
  zent_obj <- set_settings(zent_obj, alignment_dir = outdir)

  ## Return the zent object.
  return(zent_obj)

}
