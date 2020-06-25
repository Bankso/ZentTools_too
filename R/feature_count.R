
#' Count Features
#'
#' @importFrom Rsubread featureCounts
#' @importFrom stringr str_replace
#'
#' @param zent_obj Zent object.
#' @param outdir Output directory for counts.
#' @param count_feature Reads will be counted against this
#'   genomic feature type, such as the default 'exon' (GTF.featureType).
#' @param aggregate_feature Features such as 'exons' will have
#'   their counts aggregated based on this feature type, such as
#'   the default 'gene_id' (GTF.attrType).
#' @param multiple_overlap Whether to allow reads to overlap with more
#'   than one feature.
#' @param minimum_overlap The minimum number of bases that the
#'   read needs to overlap to be counted.
#' @param largest_overlap Whether the read will be assigned
#'   to the feature it has the most overlapping bases with.
#' @param read_extension_3prime The number of bases to extend the read
#'    downstream from the 3' end.
#' @param multi_mapping Whether to count multi-mapping reads.
#' @param multi_mapping_frac Whether multi-mapping reads are
#'   counted as fractions.
#' @param strand_specific Either 0 (unstranded), 1 (stranded),
#'   or 2 (reverse stranded).
#'
#' @export

count_features <- function(
  zent_obj,
  outdir = getwd(),
  count_feature = "exon",
  aggregate_feature = "gene_id",
  multiple_overlap = FALSE,
  minimum_overlap = 10,
  largest_overlap = TRUE,
  read_extension_3prime = 0,
  multi_mapping = FALSE,
  multi_mapping_frac = FALSE,
  strand_specific = 0
) {

  ## Create output directory.
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  ## Get bam files.
  bam_files <- zent_obj@sample_sheet[["bam_files"]]  

  ## Run Rsubread featureCounts.
  print_message("Feature counting from aligned reads.")
  counts <- featureCounts(
    files = bam_files,
    annot.ext = zent_obj@settings[parameter == "genome_annotation", value],
    isGTFAnnotationFile = TRUE,
    GTF.featureType = count_feature,
    GTF.attrType = aggregate_feature,
    useMetaFeatures = TRUE,
    allowMultiOverlap = multiple_overlap,
    minOverlap = minimum_overlap,
    largestOverlap = largest_overlap,
    countMultiMappingReads = multi_mapping,
    fraction = multi_mapping_frac,
    strandSpecific = strand_specific,
    isPairedEnd = as.logical(zent_obj@settings[parameter == "paired", value]),
    nthreads = as.numeric(zent_obj@settings[parameter == "ncores", value])
  )

  ## Prepare the feature counts for export.
  feature_counts <- as.data.table(
    counts$counts,
    keep.rownames = aggregate_feature
  )

  colnames(feature_counts) <- str_replace(
    colnames(feature_counts),
    "_Aligned.sortedByCoord.out.bam",
    ""
  )

  ## Export the feature counts.
  if (!str_detect(outdir, "/$")) outdir <- str_c(outdir, "/")

  fwrite(
    feature_counts, str_c(outdir, "feature_counts.tsv"), sep = "\t",
    col.names = TRUE, row.names = FALSE, quote = FALSE
  )

  ## Prepare the assignment stats for export.
  assigned <- as.data.table(counts$stat)
  
  colnames(assigned) <- str_replace(
    colnames(assigned),
    "_Aligned.sortedByCoord.out.bam",
    ""
  )

  ## Export the assignment stats.
  fwrite(
    assigned, str_c(outdir, "assignment_stats.tsv"), sep = "\t",
    col.names = TRUE, row.names = FALSE, quote = FALSE
  )

  ## Return the zent object.
  return(zent_obj)

}
