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
#' @param max_bam_sort_ram Maximum ammount of RAM per threads
#'   available for BAM sorting.
#'
#' @export

star_align <- function(
  zent_obj,
  outdir = getwd(),
  max_bam_sort_ram = "1G"
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
    "--outSAMtype", "BAM Unsorted",
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

  ## Sort and index the BAM files.
  unsorted_bams <- str_c(
    outdir,
    zent_obj@sample_sheet[["sample_name"]],
    "_Aligned.out.bam"
  )
  names(unsorted_bams) <- zent_obj@sample_sheet[["sample_name"]]

  iwalk(unsorted_bams, function(x, y) {
    command <- str_c(
      "samtools", "sort",
      "-m", max_bam_sort_ram,
      "-@", pull_setting(zent_obj, "ncores"),
      "-o", str_c(outdir, str_c(y, "_sorted.bam")),
      "-O", "BAM",
      x,
      sep = " "
    )
    system(command)

    command <- str_c(
      "samtools", "index", str_c(outdir, y, "_sorted.bam"),
      sep = " "
    )
    system(command)
  })

  ## Add the bam files to the sample sheet.
  zent_obj <- add_bams(zent_obj, alignment_dir = outdir)

  ## Add the outdir directory for alignment to the settings.
  zent_obj <- set_settings(zent_obj, alignment_dir = outdir)

  ## Return the zent object.
  return(zent_obj)

}
