make_matrix <- function(
  zent_obj,
  outdir = getwd(),
  bw = bw_in,
  bed = bed_in,
  bin_size = 1,
  ref_point = 'center',
  min_fragment = 10,
  max_fragment = 700,
  extend_reads = 200,
  before_field = 1000,
  after_field = 1000,
  temp_dir = "./temp"
) {

  ## Input checks.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")
  paired_status <- as.logical(pull_setting(zent_obj, "paired"))
  analysis_type <- pull_setting(zent_obj, "analysis_type")

  ## Make output directory if it doesn't exist.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Prepare command.
  commands <- imap(samples, function(x, y) {
    command <- str_c(
      "computeMatrix",
	  "reference-point",
      "-R", bed,
      "-S", bw,
	  "--referencePoint", ref_point,
      "-o", "matrix",
      "-bs", bin_size,
      "-b", before_field,
      "-a", after_field,
      "-p", pull_setting(zent_obj, "ncores"),
      sep = " "
    )
    
    return(command)
  })

  ## Set temporary directory.
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  Sys.setenv(TMPDIR=temp_dir)

  ## Run commands.
  print_message("Creating the coverage matrix.")
  walk(commands, system)#, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ## Return zent tools object.
  return(zent_obj)

}