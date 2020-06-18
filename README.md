# ZentTools

## Quick Start

Download a yeast genome assembly and annotation file.

```
mkdir -p genome && cd genome

wget ftp://ftp.ensembl.org/pub/release-100/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz

wget ftp://ftp.ensembl.org/pub/release-100/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz

cd ..
```

Pull the singularity container and start R.

```
singularity pull library://rpolicastro/default/zent_tools:rnaseq_v0.1

singularity exec -eCB `pwd` -H `pwd` zent_tools_rnaseq_v0.1.sif R
```

Make a sample sheet.

```
samples <- data.frame(
  sample_name = c(
    sprintf("S288C_WT_%s", seq_len(3)),
    sprintf("S288C_Dimaide_%s", seq_len(3))
  ),
  file_1 = c(
   "SRR10759407", "SRR10759408", "SRR10759409",
   "SRR10759413", "SRR10759414", "SRR10759415"
  ),
  file_2 = rep(NA, 6)
)
```

Create the ZentTools object.

```
zent <- zent_tools(
  analysis_type = "RNA-seq", sample_sheet = samples,
  paired = TRUE, ncores = 8
)
```

Since we are using SRA accession numbers for the example analysis,
we need to download the fastq files first.

```
zent <- retrieve_reads(zent, outdir = "./sequences")
```

Run FastQC for quality control.

```
zent <- fastqc(zent, outdir = "./fastqc_reports")
```

Generate the genome index for STAR alignment.

```
zent <- star_genome(
  zent, outdir = "./genome/index",
  genome_assembly = "./genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
  genome_annotation = "./genome/Saccharomyces_cerevisiae.R64-1-1.100.gtf"
)
```

Align the reads using STAR.

```
zent <- star_align(zent, outdir = "./aligned")
```
Feature counting of the assigned reads.

```
zent <- count_features(zent, outdir = "./counts", strand_specific = 2)
```
