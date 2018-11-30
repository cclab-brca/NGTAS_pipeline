#!/bin/bash 

# FASTQ TO  BAM - NOVOALIGN 

#################################### PIPELINE CONFIGURATION ####################################

# SPECIFY WORKING DIRECTORY
# NOTE: WD is the working directory where all the needed subfolders will be created. 
# 	It MUST contain a 'fastq' subfolder with the fastq files to be processed
WD=path/to/working/directory

# SPECIFY PATH TO SOFTWARE
BWA_PATH=path/to/bwa-0.7.17/
SAMTOOLS_PATH=path/to/samtools-1.8
GATK_PATH=path/to/GenomeAnalysisTK-3.6/
PICARD_PATH=path/to/picard/
JAVA18_PATH=path/to/jre1.8.0_171/bin/
NGTAS_JAVA_LIB=path/to/jar/

# SPECIFY PATH TO REQUIRED FILES	
TARGET=path/to/NGTAS_377amplicons_0based_b37.intervals
AMPLICONS=path/to/NGTAS_377amplicons_0based.bed
GENOME=path/to/human_g1k_v37_decoy.fasta
GENOME_INDEXED=path/to/indexed.bwa.human_g1k_v37_decoy.*

# OTHER INFO AND PARAMETERS
# NOTE: SLURM_ARRAY_TASK_ID and SLURM_JOB_NUM_NODES variables can have different names depending on the HPC system used
ARRAY_ID=$SLURM_ARRAY_TASK_ID
N_NODES=$SLURM_JOB_NUM_NODES

MEMORY="16g" # Memory available (gigabytes)
FASTQ1_SUFFIX=".r_1.fq.gz" # suffix of the fastq files
PLATFORM="Illumina" # Platform info to be passed to the bam file read group tag
MIN_INS_SIZE=60 # min fragment insert size required (pairs with smaller insent size will be filtered out)


################################ CREATE VARIABLES AND SUBFOLDERS ###############################

# Add software path to PATH
PATH=$BWA_PATH:$SAMTOOLS_PATH:$PATH

# Define input/output folders and create if not present
FASTQ_DIR=$WD/fastq
REPORTS_DIR=$WD/reports
REPORTS_ALIGNMENT_DIR=$REPORTS_DIR/alignment
REPORTS_ALIGNSTATS_DIR=$REPORTS_DIR/alignment.stats
REPORTS_BAMSTATS_DIR=$REPORTS_DIR/bam.stats 
BAM_DIR=$WD/bam 
OUTPUT_AMPLISTATS=$WD/amplistats/
DIRS=$$REPORTS_DIR:$REPORTS_ALIGNMENT_DIR:$REPORTS_BAMSTATS_DIR:$REPORTS_ALIGNSTATS_DIR:$BAM_DIR:$OUTPUT_AMPLISTATS
ARR=$(echo $DIRS | tr ":" "\n")
for x in $ARR
do
        if [ ! -d "$x" ]; then
        	mkdir -p $x
        fi
done

# Extract info for the fastq pair to be processed
FASTQ1_PATH=$(ls $FASTQ_DIR/*$FASTQ1_SUFFIX | sed -n "$ARRAY_ID"p)
FASTQ2_SUFFIX="${FASTQ1_SUFFIX/1/2}"
FASTQ2_PATH=${FASTQ1_PATH/$FASTQ1_SUFFIX/$FASTQ2_SUFFIX}
FASTQ1_FILENAME=${FASTQ1_PATH##*/}
SAMPLE_NAME=${FASTQ1_FILENAME%$FASTQ1_SUFFIX*}
# NOTE: Information extracted from SAMPLE_NAME need to be adapted to the contingent filename structure
LIBRARY=$(echo $SAMPLE_NAME | cut -d '.' -f 1)
BARCODE=$(echo $SAMPLE_NAME | cut -d '.' -f 2)
FLOWCELL=$(echo $SAMPLE_NAME | cut -d '.' -f 3)
LANE=$(echo $SAMPLE_NAME | cut -d '.' -f 4)


######################################## PIPELINE ############################################

echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Starting alignment pipeline"

# Create a local temporary folder where all input and output files will be saved
LOCAL="/tmp/Novo_$RANDOM"
mkdir "$LOCAL"
cp $GENOME_INDEXED $LOCAL
cp $FASTQ1_PATH $LOCAL
cp $FASTQ2_PATH $LOCAL

# Alignment using BWA-MEM
bwa mem -t $N_NODES \
	-r 1.5 \
	-R "@RG\tID:$SAMPLE_NAME\tPL:$PLATFORM\tPU:"$FLOWCELL.$LANE"\tLB:$LIBRARY\tSM:$BARCODE" \
	-M \
	$LOCAL/indexed.bwa.human_g1k_v37_decoy \
	$LOCAL/*$FASTQ1_SUFFIX \
	$LOCAL/*$FASTQ2_SUFFIX > $LOCAL/$SAMPLE_NAME.sam
samtools view -b -S $LOCAL/$SAMPLE_NAME.sam > $LOCAL/$SAMPLE_NAME.bam

# Sorting to compute stats
samtools sort $LOCAL/$SAMPLE_NAME.bam -o  $LOCAL/$SAMPLE_NAME.sorted.bam

# Compute alignment and bam statistics
echo "[exome:merge_bam_and_bamstats $(date +"%Y-%m-%d %T")] Starting CollectAlignmentMetrics"
"$JAVA18_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
	CollectAlignmentSummaryMetrics \
        R=$GENOME \
        I=$LOCAL/$SAMPLE_NAME.sorted.bam \
        O=$REPORTS_ALIGNSTATS_DIR/$SAMPLE_NAME.alignMetrics.txt

"$JAVA18_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
	CollectHsMetrics \
	BAIT_INTERVALS=$TARGET \
	TARGET_INTERVALS=$TARGET \
	I=$LOCAL/$SAMPLE_NAME.sorted.bam \
	O=$REPORTS_BAMSTATS_DIR/$SAMPLE_NAME.bamMetrics.txt \
	REFERENCE_SEQUENCE=$GENOME


# Filtering out reads with insert size smaller than 60 and sorting
samtools view -h -f 2 $LOCAL/$SAMPLE_NAME.bam | awk 'substr($0,1,1) == "@" || $9 > '$MIN_INS_SIZE' || $9 < -'$MIN_INS_SIZE' {print}' | samtools view -bS - > $LOCAL/$SAMPLE_NAME.filt.bam
samtools sort $LOCAL/$SAMPLE_NAME.filt.bam -o $LOCAL/$SAMPLE_NAME.filt.sorted.bam
samtools index $LOCAL/$SAMPLE_NAME.filt.sorted.bam


# Local Realignment
echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Running Local Realignment"

"$JAVA18_PATH"java -Xmx$MEMORY -jar "$GATK_PATH"GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R $GENOME \
   -I $LOCAL/$SAMPLE_NAME.filt.sorted.bam \
   -o $LOCAL/$SAMPLE_NAME.intervals \
   -nt $N_NODES
   
"$JAVA18_PATH"java -Xmx$MEMORY -jar "$GATK_PATH"GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R $GENOME \
   -I $LOCAL/$SAMPLE_NAME.filt.sorted.bam \
   -targetIntervals $LOCAL/$SAMPLE_NAME.intervals \
   -o $LOCAL/$SAMPLE_NAME.realn.bam


# Bam annotation - assign each read to the corresponding amplicon
"$JAVA18_PATH"java -Xmx$MEMORY -cp $NGTAS_JAVA_LIB/ngtas.jar:$NGTAS_JAVA_LIB/commons-cli-1.2.jar ngtas.AnnotateBam \
	-i $LOCAL/$SAMPLE_NAME.realn.bam \
	-o $BAM_DIR \
	-s $SAMTOOLS_PATH \
	-b $AMPLICONS \
	-t 3
# NOTE: -t option indicates how many bases the read starting point can differ from the amplicon starting point

# Compute per position coverage for each amplicon
"$JAVA18_PATH"java -Xmx$MEMORY -cp $NGTAS_JAVA_LIB/ngtas.jar:$NGTAS_JAVA_LIB/commons-cli-1.2.jar ngtas.AmpliStats \
	-i $BAM_DIR/$SAMPLE_NAME.realn.annotated.bam \
	-o $OUTPUT_AMPLISTATS \
	-s $SAMTOOLS_PATH \
	-r $GENOME \
	-b $AMPLICONS


# Remove temp files
rm -rf "$LOCAL"

echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Alignment complete"

# COPY SCRIPT
# scp /Users/callar01/OneDrive\ -\ Cancer\ Research\ UK\,\ Cambridge\ Institute/CRIdocs/ctDNA/pipeline_example/Scripts/1b-fastq_to_bam_BWA.sh callar01@clust1-headnode.cri.camres.org:/Users/callar01/ctDNA/pipeline_example

# RUN JOBS
# sbatch --array=1-9%9 --job-name ngtas --mem 16000 -n 8 --output /scratchb/cclab/callar01/log/ngtas.%A.%a.out.txt --error /scratchb/cclab/callar01/log/ngtas.%A.%a.err.txt /Users/callar01/ctDNA/pipeline_example/1b-fastq_to_bam_BWA.sh 

# srun --pty --mem 16000 -n 8 /usr/bin/bash


