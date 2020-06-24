
#' Zent Class
#'
#' @import data.table
#'
#' @slot sample_sheet Sample sheet.
#' @slot settings Various run settings.
#'
#' @export

setClass(
  "zent_obj",
  representation(
    sample_sheet = "data.table",
    settings = "data.table"
  ),
  prototype(
    sample_sheet = data.table(),
    settings = data.table()
  )
)

#' Zent Object Constructor Function
#'
#' @import methods
#'
#' @param sample_sheet Either a data.frame or delimited file.
#' @param analysis_type Either 'RNA-seq' or 'ChIP-seq'.
#'   Only 'RNA-seq' is currently available.
#' @param sep If the sample sheet is a file, this specifies the delimiter.
#' @param paired Whether the run is paired (TRUE) or unpaired (FALSE)
#' @param ncores The number of cores/threads to use.
#'
#' @return A zent object.
#'
#' @export

zent_tools <- function(
  analysis_type,
  sample_sheet,
  sep = "\t",
  paired = FALSE,
  ncores = 1
) {

  ## Prepare the sample sheet.
  if (is(sample_sheet, "character")) {
    samples <- fread(sample_sheet, sep = sep)
  } else if (is(sample_sheet, "data.frame")) {
    samples <- as.data.table(sample_sheet)
  }

  ## prepare run settings.
  run_settings <- data.table(
     parameter = c(
       "analysis_type", "paired", "ncores", "genome_dir",
       "genome_annotation", "genome_assembly", "alignment_dir"
    ),
     value = c(analysis_type, paired, ncores, rep("", 4))
  )

  ## Create the zent object.
  zent_obj <- new(
    "zent_obj",
    sample_sheet = samples,
    settings = run_settings
  )

  return(zent_obj)

}
