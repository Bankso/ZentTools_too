# ZentTools

This repository houses ZentTools, a collection of pipelines for RNA-seq,
  ChIP-seq (not available yet), and ChEC-seq (not available yet) data analysis.
All required software is housed within a Singularity container to aid
  reproducibility and extensibility.
The code itself is written similarly to an R package and calls outside software
  when required.

## Preparing Software

### Installing Singularity

Before running the pipeline, Singularity must be installed.
Detailed instructions can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps).

A Singularity container packages all software required to run a given pipeline, facilitating the pipeline's use across diverse systems, standardizing software versions to enhance reproducibility, and preventing pipelines from breaking after software updates.

## Example Workflow

For this example workflow, we will be running the RNA-seq pipeline on six budding yeast samples.
Three of the samples are untreated and the other three are treated with the thiol-reactive oxidizing agent diamide, which results in a large-scale transcriptional stress response.

### Obtaining Genome Assembly and Annotation

For this example, it is recommended to create a clean working directory
  that contains all reference files and will house the resulting outputs.

Before starting the analysis, you first need to download a yeast genome assembly
  and an annotation file to your working directory.
Here, we obtain the budding yeast S288C assembly and annotation from ENSEMBL.

```
mkdir -p genome && cd genome

wget ftp://ftp.ensembl.org/pub/release-100/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz

wget ftp://ftp.ensembl.org/pub/release-100/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz

cd ..
```

### Starting Singularity

Once the assembly and annotation have been obtained, the Singularity container can be downloaded and the internal R version started.

```
singularity pull library://rpolicastro/default/zent_tools:rnaseq_v0.1
singularity exec -eCB `pwd` -H `pwd` zent_tools_rnaseq_v0.1.sif R
```

### Creating the ZentTools Object

To help keep track of samples, the workflow requires a sample sheet.
The sheet should have three columns: 'sample_name', 'file_1', and 'file_2'.
'sample_name' is the name of the sample, 'file_1' is either an NCBI SRA accession number
  or the R1 read from the sequencer (including the path to the file location),
  and 'file_2' is the R2 read in the case of paired end reads (also including
  the path to the file).
If there is no R2 file, the 'file_2' column can be filled with blank values or NAs. 

The sample sheet can be either a file read into R or a data.frame constructed within R (demonstrated here).

```
samples <- data.frame(
  sample_name = c(
    sprintf("S288C_WT_%s", seq_len(3)),
    sprintf("S288C_Diamide_%s", seq_len(3))
  ),
  file_1 = c(
   "SRR10759407", "SRR10759408", "SRR10759409",
   "SRR10759413", "SRR10759414", "SRR10759415"
  ),
  file_2 = rep(NA, 6)
)
```

Once a sample sheet has been created, the ZentTools object is constructed.

```
library("ZentTools")

zent <- zent_tools(
  analysis_type = "RNA-seq", sample_sheet = samples,
  paired = TRUE, ncores = 8
)
```

### Download Reads from SRA

Since we are using files from the SRA in this example,
  they need to be downloaded using sra-tools.
This step can be skipped if not using files from the SRA.

```
zent <- retrieve_reads(zent, outdir = "./sequences")
```

### FASTQ Quality Control

For all sequencing runs, it is recommended to run a quality control
  analysis on the FASTQ files to check for any problems in
  sample preparation or sequncing.
Quality control is done using the FastQC software.

```
zent <- fastqc(zent, outdir = "./fastqc_reports")
```

### Aligning Reads

The STAR aligner remains the gold standard for spliced-read alignment,
  so is utilized in this pipeline.
Before aligning the FASTQ files, a genome index must first be generated
  using the previously downloaded genome assembly and annotation files.

```
zent <- star_genome(
  zent, outdir = "./genome/index",
  genome_assembly = "./genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
  genome_annotation = "./genome/Saccharomyces_cerevisiae.R64-1-1.100.gtf"
)
```

Once the index is created, you are ready to align the reads to the genome.

```
zent <- star_align(zent, outdir = "./aligned")
```

### Feature Counting

Finally, after the reads are aligned, a feature count matrix can be calculated.
The pipeline uses featureCounts from Rsubread to assign the aligned reads
  to the gene (or transcipt) that they overlap.
The number of reads per gene (or transcript) is then counted.

This example data were generated from libraries prepared using the dUTP method; thus, the data are strand-specific, with reverse strandedness. That is specified by setting strand_specific to 2.
More details on strand specificity can be found in the documentation for the function.

```
zent <- count_features(zent, outdir = "./counts", strand_specific = 2)
```

## Submitting to HPC

This pipeline can easily be run on an HPC.
First, put all of the R commands into a file with an '.R' extension, such as 'run_pipeline.R'.
You then need to make a short bash script to run everything.
An example is provided here.

```
#!/bin/bash

# Change to the working directory (replace with the path to yours).
cd ~/workdir

# Ensure that Singularity is available
# For Indiana HPC users, it would be 'module load singularity'.

# Use the Singularity container to run your script.
singularity exec -eCB `pwd` -H `pwd` Rscript run_pipeline.R
```

## Notes on Directories

If any directories you need to access are not in your working directory, the singularity exec command must be slightly modified.
The `-B` part of the singularity command **B**inds directories to the container
  so they can be accessed.
You can add whatever directories you want by giving it a comma separated
  list of directories (without spaces).
