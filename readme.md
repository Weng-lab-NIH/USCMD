### Installing the Pipeline
You will need to create singularity images for each step. In each step's folder, there is a file named `step<number>.def`, which is a singularity definition file. You can use this file to create that step's singularity file by running the corresponding script named `build_step<number>.bash`. This script uses sylab's remote container building service, so you'll have to set up an account for that [here](https://cloud.sylabs.io/builder). Of course, you can build the container locally on your own machine if you have sudo permissions on it.

### Generating Synthetic Data
We have provided the metadata files and scripts needed to generate the synthetic reads data we used to test this pipeline in `syn_pipeline`. To generate this data yourself, you will need to have samtools and bwa installed. Then, compile `_2.0.0_generate_exome_reads_v2.cpp`. You should then be able to run `run_syn_data_pipeline_csv_input.bash` to get the exome reads files and possorted bam file needed to run the actual mutation calling pipeline. 

To generate your own synthetic donor, add three files to the `syn_pipeline/donor_metadata_csv` folder, named `<donor>_single_cell.csv`, `<donor>_SNPs.csv`, and `<donor>_single_cell.csv`. The donor_info file just contains age and sex info on the donor. The SNPs file specifies which chromosome and position SNPs should be itnroduced at. Note that the subsetted reference we used only has a single chromosome, so the chromosome ID is always 0. 

The single_cell file is the most complicated. For each donor, we generate multiple cells, each with their own cell id. Each cell has multiple areas_of_interest,  each with their own corresponding ID, chromosome number and position. At each area_of_interest we can have several UMIs (with their own ids). Each UMI can have multiple read_sets (again, each read set has its own ID). Each read_set has the same nucleotide specified for at the position we set for the area_of_interest. We can specify the number of reads generated for each read set.

**Please note the donor_id column must correspond with the names of the csv files.**

### Running the Pipeline
To run USCMD, you will need the following data:
- a reference human exome file (fasta format)
- two exome reads files (r1 and r2 in fastq.gz format)
- a possorted genome file output by cellranger (bam format)
- a list of barcodes identifying which cells in the above file to look for mutations in (.txt format, a barcode per line)

The following is an example call to the pipeline on synthetic donor SYN_M1:
```
run_pipeline.bash \
    --sample SYN_M1 \
    --barcode_list ./syn_data/sc_reads_v2/SYN_M1/used_barcodes.txt \
    --num_cores 2 \
    --ref_file ./syn_data/subsetted_reference/subset_genome.fa \
    --exome_dir ./syn_data/synthetic_exome_reads/SYN_M1 \
    --r1_filename SYN_M1_R1.fastq.gz \
    --r2_filename SYN_M1_R2.fastq.gz \
    --cellranger_bam ./syn_data/possorted_genomes_v2/SYN_M1/SYN_M1.bam \
    --pipeline_dir test_run_210928a/SYN_M1
```

| Argument Name | Description |
| -- | -- |
| sample | arbitrary sample name |
| barcode list | list of barcodes identifying which cells in the above file to look for mutations in (.txt format, a barcode per line) |
| num_cores | number of cpu cores to use |
| ref_file | reference human exome file (fasta format) |
| exome_dir | the directory containg the two exome reads files |
| r1_filename | the name of the R1 file within exome_dir |
| r2_filename | the name of the R2 file within exome_dir |
| cellranger_bam | possorted.bam output from cellranger |
| pipeline_dir | an empty directory to write output to |

Please note `run_pipeline.bash` assumes each step's container is contained in a folder named `pipeline_containers`.

We developed a dashboard to allow data visualization. To display this dashboard, you will need to input an HTTP port. By default, this is 80 on most systems. You will also need to input a CSV file with columns for sample_name, sex, age, and pipeline_dir. Here's an example of how to run the dashboard:

```
Rscript ./USCMD\ Dashboards.R \
  --donor_list_csv ./synthetic_donor_list.csv \
  --port 80 \
```

synthetic_donor_list.csv is in the dashboard folder of this repo.

### How it works (overview):
This method was implemented in ten steps:
(1) Align and Index Exome carries out the alignment of individual exome sequencing reads to the human reference genome to create a donor exome bam file. 

In parallel, (2) Filter Single Cell Data is used to split the combined aligned scRNAseq data output by the CellRanger pipeline (named possorted_genome_bam.bam in the cellranger output folder) into separate files for each of the individual cells.  

(3) Identify Exome SNPs is used to identify variants based on this individual exome.

(4) Create Donor-Specific Reference is used to filter these variants by checking if they are single nucleotide polymorphisms (SNPs), which were then used to create a donor-specific reference in which their SNPs are incorporated into the individual reference. 

(5) Call and Filter Cell Variants is used to call somatic mutations within the single cell bam files generated by the Filter Single Cell Data step based on the donor-specific reference.  

Variants can be classified as silent, missense, or nonsense mutations based on whether they lead to a change of amino acid or generate a stop codon. Mutations with known functional changes are then annotated in (6) Functional and Genetic Annotation of Variants.

This is followed by (7) Reformat Annotated Variants in which the read information associated with each variant call in the Call and Filter Cell Variants step is extracted. 

(8) UMI corrections is carried out on these cell variants and creates an output TSV (spreadsheet) file. 

(9) Score Mutations for single cells takes in these TSV files and scores the mutations based on multiple quality control metrics and outputs these scores as a CSV file. 

(10) Summarize Donor Exome and Mutations (10a) and Summarize Single-Cell Stats (10b) iterates through the donor inputs to the pipeline and outputs the exome depth, sequencing coverage, and sequencing reads for all the donors to a single CSV file, which can be used to adjust for differences in sequencing depth and exome coverage among the donors.
