STEP1 README

TITLE: Read alignment
PURPOSE: Align exome sequencing reads to a reference genome, sort them, 
and index them to create a bam file.
SUMMARY: Two fastq.gz files are taken in - read 1 and read 2. These are trimmed
using trim_galore. The reference genome (in fasta format) is indexed using the
"bwa index" command. The trimmed reads files are then aligned to the indexed 
reference genome to make a sam file. Using samtools, the output is sorted and 
turned into a bam file, which is then indexed using samtools. A unique 
read group is then added to the header using gatk, and the file reindexed with 
samtools.

PACKAGES:
    bamtools
    samtools
    bowtie
    GATK/4.1.3.0
    bwa
    trimgalore

INPUT FILES: 
    1: read1 sequencing file (fastq.gz)
    2: read2 sequencing file (fastq.gz)
    3: reference genome (fasta)

INPUT DIRECTORIES:
    1: exome_data: directory contaning reads files.
    2: aligned_dir: output directory

PARAMETERS:
    1: num_cores: Number of cpu cores. Recommended max of 8 cores.

MAIN OUTPUT FILES:
    1: Output bam file (<sample_name>_SM_bwa.bam)
    2: Output bam index (<sample_name>_SM_bwa.bam.bai)

OTHER OUTPUT FILES:
    1: Interim bam and sam files in output directory
    2: Trimmed reads files and corresponding trimming reports in folder 
        containing input sequencing files.
    3: Metadata files from indexing genome, placed in folder containing 
        reference genome fasta file.
