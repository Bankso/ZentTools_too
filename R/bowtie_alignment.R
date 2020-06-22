
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
    str_c(outdir, index_name)
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
