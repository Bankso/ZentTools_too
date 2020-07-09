
#' Add BAMs to Sample Sheet
#'
#' @param zent_obj ZentTools object.
#' @param alignment_dir Directory of aligned reads.
#'
#' @export

add_bams <- function(
  zent_obj,
  alignment_dir
) {

  ## Grab some info from object and prepare inputs.
  analysis_type <- pull_setting(zent_obj, "analysis_type")
  if (!str_detect(alignment_dir, "/$")) {
    alignment_dir <- str_c(alignment_dir, "/")
  }
  sample_sheet <- copy(zent_obj@sample_sheet)

  ## Add BAMs if RNA-seq experiment.
  if (analysis_type == "RNA-seq") {
    sample_sheet[, bam_files := str_c(
      alignment_dir,
      sample_name,
      "_Aligned.sortedByCoord.out.bam"
    )]
  } else if (analysis_type %in% c("ChIP-seq", "ChEC-seq")) {
    sample_sheet[,
      sample_bams := str_c(alignment_dir, sample_name, ".bam")
    ]
    sample_sheet[
      !is.na(control_file_1),
      control_bams := str_c(alignment_dir, control_name, ".bam")
    ]
  }

  ## Add new sample sheet back to zent object.
  zent_obj@sample_sheet <- sample_sheet

  ## Return the zent object.
  return(zent_obj)
}
