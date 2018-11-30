# NGTAS_pipeline
Computational pipeline to analyse Next Generation Targeted Amplicon Sequence data

This is the computational pipeline for the analysis of Next Generation-Targeted Amplicon Sequencing (NG-TAS) data. A detailed description of the approach can be found in our manuscript: https://www.biorxiv.org/content/early/2018/07/15/366534. Required files are provided and raw data for our NA12878 benchmark dataset can be found here: https://figshare.com/articles/NGTAS_NA12878/7387370.

The pipeline is divided in two parts:
1.	Alignment and bam annotation. Set up to run for a group of demultiplexed fastq files as an array job. Each job will process a pair of files (read1 and read2 fastq) from a single replicate of a single sample. The pipeline will align the data and assign each pair to the amplicon they belong to, enabling the subsequent per-amplicon analyses. Two versions of the pipeline are available, using either Novoalign (1a) or BWA-MEM (1b) as aligners. 
In both cases, the first session of the script contains all the configuration commands that should be set by the user. In particular:

a)	Working directory, where all the needed subfolders will be created. It MUST contain a 'fastq' subfolder with the fastq files to be processed
b)	Path to all required software. This include: Novoalign/Novosort (v. 3) (or BWA-MEM v. 0.7.17), Samtools (v. 1.8), GATK (v. 3.6), Picard (v. 2.17), Java (v. 1.8), AnnotateBam and Amplistats java tools (provided).
c)	Path to required files. The files describing the 377 amplicons used in the paper are provided. They can be used as example to describe any set of amplicons. Reference genome is not provided.
d)	Other info and parameters like the ARRAY_ID and the number of nodes (N_NODES) typically acquired as environmental variable (e.g. SLURM_ARRAY_TASK_ID and SLURM_JOB_NUM_NODES in a Slurm HPC).

Please note that LIBRARY, BARCODE, FLOWCELL and LANE are extracted from the file name. If the filename structure is different, this part of the code has to be modified accordingly.
Output: 1) a filtered and annotated bam file; 2) two or three reports containing alignment and bam statistics; 3) mpileup-like output containing the coverage at each position for the four nucleotides. Stats are derived for each amplicon independently.


2.	Mutation calling. Set up to run as an array job for a set of sample groups (i.e. patients) as specified in the pairs_for_mutation_calling.txt spreadsheet (see below). Each job will process all samples belonging to the same group/patient and return an annotated .vcf containing details for all mutations called in at least one sample.
The first session of the script contains all the configuration commands that should be set by the user. In particular:

a)	Working directory, same as in part 1.
b)	Path to all required software. This include: Java (v. 1.8), Samtools (v. 1.8), GATK (v. 3.6), Picard (v. 2.17), Bedtools (v. 2.27), bcftools (v. 1.9), bgzip (included in HTSlib, v.1.8), Annovar (including annotation databases refGene, cytoBand, genomicSuperDups, esp6500siv2_all, snp138, ljb26_all), R (v. 3.5 with the package ‘VariantAnnotation’ installed).
c)	Path to required files. These include the amplicon description files similarly to part 1 (but excluding the primer regions) plus a set of R scripts (provided) and the pairs_for_mutation_calling.txt spreadsheet.
This spreadsheet must contain 8 columns with the following info:
•	Tumour/Plasma replicate file name (without extension)
•	Tumour/Plasma replicate ID
•	Normal replicate file name (without extension)
•	Normal replicate ID
•	Sample ID
•	Group/Patient ID
•	Study name (can be blank, not used by the pipeline)
•	DNA source (it is important that it is specified as ‘FFPE’ to enable the specific C>T/G>A filter in FFPE samples
An example spreadsheet to analyse the set of example raw data is provided.
d)	Other info and parameters like the ARRAY_ID and the number of nodes (N_NODES) typically acquired as environmental variable (e.g. SLURM_ARRAY_TASK_ID and SLURM_JOB_NUM_NODES in a Slurm HPC). The user can also specify the desired thresholds for the following mutation filtering parameters:
•	vaf_th_norm: max VAF for the alternative allele in the normal sample
•	vaf_th_tum: min VAF to call a somatic mutation
•	ratio: min Tumour (or plasma)/Normal VAF ratio to call a somatic mutation
•	cov_th_norm: min coverage (in a specific position) in the normal to call a somatic mutation
•	cov_th_tum: min coverage (in a specific position) in the tumour/plasma to call a somatic mutation
•	th_ffpe: min VAF to keep a C>T/G>A mutation in FFPE samples
Output: 1) txt file with the list of mutations called in each sample; 2) annotated .vcf file.
