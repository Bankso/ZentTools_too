

Prepare a clean working directory and download the yeast genome assembly and annotation.

```
mkdir -p genome && cd genome

wget ftp://ftp.ensembl.org/pub/release-100/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz

wget ftp://ftp.ensembl.org/pub/release-100/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz

cd ..
```

Start singularity.

```
singularity pull --arch amd64 library://zentlab/default/zent_tools:software_v0.2
singularity exec -eCB `pwd` -H `pwd` zent_tools_software_v0.2.sif R
```

Create a sample sheet.

```
samples <- data.frame(
  sample_name = sprintf("Reb1_30s_%s", seq_len(3)),
  file_1 = sprintf("SRR19478%s", seq(19, 21, 1)),
  file_2 = rep(NA, 3),
  control_name = rep("Free_MNase_30s", 3),
  control_file_1 = rep("SRR1947777", 3),
  control_file_2 = rep(NA, 3)
)
```

Make the ZentTools object.

```
library("ZentTools")

zent <- zent_tools(
  analysis_type = "ChEC-seq", sample_sheet = samples,
  paired = TRUE, ncores = 8
)
```

Retrieve the reads from the SRA.

```
zent <- retrieve_reads(zent, outdir = "./sequences")
```

Quality control of the reads.

```
zent <- fastqc(zent, outdir = "./fastqc_reports")
```

Make a bowtie2 genome index.

```
zent <- bowtie2_index(
  zent, outdir = "./genome/index",
  genome_assembly = "./genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
  index_name = "S288C_index"
)
```

Align the reads using bowtie2.

```
zent <- bowtie2_align(zent, outdir = "./aligned")
```

Make bigwigs.

```
zent <- make_bigwigs(
  zent, outdir = "./bigwigs", bin_size = 25,
  normalize_using = "CPM"
)
```

Call peaks.

```
zent <- call_peaks(zent, outdir = "./peaks", genome_size = "12e6")
```

Annotate the peaks.

```
zent <- annotate_peaks(
  zent, outdir = "./peaks",
  genome_annotation = "./genome/Saccharomyces_cerevisiae.R64-1-1.100.gtf"
)
```
