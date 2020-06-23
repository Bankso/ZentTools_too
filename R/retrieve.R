
#' Retrieve Reads from SRA
#'
#' @importFrom stringr str_c str_detect
#'
#' @param zent_obj Zent object.
#' @param outdir Directory to donwload files to.
#'
#' @export

retrieve_reads <- function(
  zent_obj,
  outdir = getwd()
) {

  ## Input checks.
  analysis_type <- zent_obj@settings[parameter == "analysis_type", value]
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  ## Make sure outdir exists.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get vector of accession numbers.
  accessions <- zent_obj@sample_sheet[["file_1"]]

  if (analysis_type %in% c("ChIP-seq", "ChEC-seq")) {
    controls <- unique(zent_obj@sample_sheet[["control_file_1"]])
    accessions <- str_c(accessions, controls)
  }

  ## Prepare fasterq-dump command.
  paired_status <- as.logical(zent_obj@settings[parameter == "paired", value])

  command <- str_c(
    "fasterq-dump",
    "-O", outdir,
    "-e", zent_obj@settings[parameter == "ncores", value],
    sep = " "
  )

  if (paired_status) {
    command <- str_c(command, "--split-files", sep = " ")
  }

  command <- str_c(
    command,
    str_c(accessions, collapse = " "),
    sep = " "
  )

  ## Run the fasterq-dump command.
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Update the sample sheet.
  sample_sheet <- copy(zent_obj@sample_sheet)

  if (paired_status) {
    sample_sheet[, c("file_1", "file_2") := list(
      str_c(outdir, file_1, "_1.fastq"),
      str_c(outdir, file_1, "_2.fastq")
    )]

    if (analysis_type %in% c("ChIP-seq", "ChEC-seq")) {
      sample_sheet[, c("control_file_1", "control_file_2") := list(
        str_c(outdir, control_file_1, "_1.fastq"),
        str_c(outdir, control_file_2, "_2.fastq")
      )]
    }
  } else {
    sample_sheet[, file_1 := str_c(outdir, file_1, ".fastq")]

    if (analysis_type %in% c("ChIP-seq", "ChEC-seq")) {
      sample_sheet[,
        control_file_1 := str_c(outdir, control_file_1, ".fastq")
      ]
    }
  }

  zent_obj@sample_sheet <- sample_sheet

  ## Return the zent object.
  return(zent_obj)

}
