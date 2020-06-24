
#' Set Settings.
#'
#' @param zent_obj ZentTools object.
#' @param analysis_type Type of experiment.
#'   Either 'RNA-seq', 'ChIP-seq', or 'ChEC-seq'.
#' @param paired Paired-end status of the reads.
#'   Either TRUE or FALSE.
#' @param ncores The number of cores/threads to use.
#' @param genome_dir The directory of the genome index.
#' @param genome_annotation The directory and file name of the
#'   the genome index GTF file.
#' @param genome_assembly The directory and file name of the
#'   the gnome assembly FASTA file.
#' @param alignment_dir The directory containing the aligned reads.
#'
#' @export

set_settings <- function(
  zent_obj,
  analysis_type = NA,
  paired = NA,
  ncores = NA,
  genome_dir = NA,
  genome_annotation = NA,
  genome_assembly = NA,
  alignment_dir = NA
) {

  settings <- copy(zent_obj@settings)

  if (!is.na(analysis_type)) {
    settings[parameter == "analysis_type", value := analysis_type]
  }

  if (!is.na(paired)) {
    settings[parameter == "paired", value := paired]
  }

  if (!is.na(ncores)) {
    settings[parameter == "ncores", value := ncores]
  }

  if (!is.na(genome_dir)) {
    settings[paremter == "genome_dir", value := genome_dir]
  }

  if (!is.na(genome_annotation)) {
    settings[parameter == "genome_annotation", value := genome_annotation]
  }

  if (!is.na(genome_assembly)) {
    settings[parameter == "genome_assembly", value := genome_assembly]
  }

  if (!is.na(alignment_dir)) {
    settings[parameter == "alignment_dir", value := alignment_dir]
  }

  zent_obj@settings <- settings
  return(zent_obj)
}

#' Pull Settings
#'
#' @param zent_obj ZentTools object.
#' @param setting Setting to pull.
#'
#' @export

pull_setting <- function(
  zent_obj,
  setting
) {

  setting_value <- zent_obj@settings[parameter == setting, value]
  return(setting_value)

}
