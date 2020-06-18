
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
