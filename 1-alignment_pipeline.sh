#!/bin/bash 

# ALIGNMENT PIPELINE FOR NG-TAS DATA

#################################### INPUT INFO ####################################
# DIRECTORIES
WD=/path/to/working/directory
TEMP_DIR=$WD/temp
FASTQ_DIR=$WD/fastq
ANALYSIS_SUFFIX='.novo'
REPORTS_DIR=$WD/reports$ANALYSIS_SUFFIX 
REPORTS_ALIGNMENT_DIR=$REPORTS_DIR/alignment
REPORTS_ALIGNSTATS_DIR=$REPORTS_DIR/alignment.stats
REPORTS_BAMSTATS_DIR=$REPORTS_DIR/bam.stats 
BAM_DIR=$WD/bam$ANALYSIS_SUFFIX 
OUTPUT_MPILEUP=$WD/mpileup/
EXOME_DIRS=$TEMP_DIR:$REPORTS_DIR:$REPORTS_ALIGNMENT_DIR:$REPORTS_BAMSTATS_DIR:$REPORTS_ALIGNSTATS_DIR:$BAM_DIR:$OUTPUT_MPILEUP

# Ensure all paths exist. Create directory structure if not present
ARR=$(echo $EXOME_DIRS | tr ":" "\n")
for x in $ARR
do
        if [ ! -d "$x" ]; then
        	mkdir -p $x
        fi

done

# REQUIRED FILES	
TARGET=/path/to/NGTAS_377amplicons_0based_b37.intervals
AMPLICONS=/path/to/NGTAS_377amplicons_0based.bed
GENOME=/path/to/human_g1k_v37_decoy.fasta
GENOME_INDEXED=/path/to/indexed.human_g1k_v37_decoy.fasta

# SOFTWARE PATH
NOVO_PATH=/path/to/novocraft
SAMTOOLS_PATH=/path/to/samtools-1.1/bin
GATK_PATH=/path/to/GenomeAnalysisTK-3.5/
PICARD_PATH=/path/to/picard-tools-1.140/
JAVA17_PATH=/path/to/jre1.7.0_67/bin/
JAVA18_PATH=/path/to/jdk1.8.0_31/bin/
PATH=$NOVO_PATH:$SAMTOOLS_PATH:$PATH
CLASSPATH=.:/lustre/cclab/sammut01/java-code/tamseq-pipe/jar/commons-cli-1.2.jar:/lustre/cclab/sammut01/java-code/tamseq-pipe/jar/audrie.jar:/lustre/cclab/sammut01/java-code/tamseq-pipe/bin/

# OTHER INFO
N_CORES=8			# Number of cores available
MEMORY="16g"			# Available memory
FASTQ1_SUFFIX=".r_1.fq.gz"	# Suffix of FastQ file names
PLATFORM="Illumina"
####################################################################################

echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Starting alignment pipeline"
 
# fastq variables
# NOTE: LSB_JOBINDEX variable can have different names depending on the HPC system used
#	Information extracted from the filename need to be adapted to the contingent filename structure

FASTQ1_PATH=$(ls $FASTQ_DIR/*$FASTQ1_SUFFIX | sed -n "$LSB_JOBINDEX"p)		
FASTQ2_SUFFIX="${FASTQ1_SUFFIX/1/2}"
FASTQ2_PATH=${FASTQ1_PATH/$FASTQ1_SUFFIX/$FASTQ2_SUFFIX}
FASTQ1_FILENAME=${FASTQ1_PATH##*/}
SAMPLE_NAME=${FASTQ1_FILENAME%$FASTQ1_SUFFIX*}

LIBRARY=$(echo $SAMPLE_NAME | cut -d '.' -f 1)
BARCODE=$(echo $SAMPLE_NAME | cut -d '.' -f 2)
FLOWCELL=$(echo $SAMPLE_NAME | cut -d '.' -f 3)
LANE=$(echo $SAMPLE_NAME | cut -d '.' -f 4)


# Alignment
LOCAL_NOVO="/tmp/Novo_$RANDOM"
mkdir "$LOCAL_NOVO"
cp $GENOME_INDEXED $LOCAL_NOVO
cp $FASTQ1_PATH $LOCAL_NOVO
cp $FASTQ2_PATH $LOCAL_NOVO

echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Running Novoalign"
novoalign -d $LOCAL_NOVO/*.fasta \							# reference genome
	-f $LOCAL_NOVO/*$FASTQ1_SUFFIX $LOCAL_NOVO/*$FASTQ2_SUFFIX \			# input files
	-o SAM $'@RG\tID:readgroup\tSM:sample\tPU:platformunit\tLB:library' \		# output format and @RG annotation
	-l 30 \										# minimum number of matching bases
	-a ACACTGACGACATGGTTCTACA TACGGTAGCAGAGACTTGGTCT \				# adapter trimming
	-i PE 30-300 \									# paired end option
	-t 300 \									# alignment score threshold
	-c $N_CORES \									# multithreading
	-n 80 \										# hard trimming
	-k \										# base recalibration (on)
	-o FullNW \									# soft clipping
	2>$LOCAL_NOVO/$SAMPLE_NAME.novostats.txt | \					# novoalign report
samtools view -1 -bS - > $LOCAL_NOVO/$SAMPLE_NAME.bam


# SORTING
echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Running Novosort"
novosort  -m $MEMORY \
	-t $LOCAL_NOVO \
	-c $N_CORES \
	$LOCAL_NOVO/$SAMPLE_NAME.bam \
	-i \
	-o $LOCAL_NOVO/$SAMPLE_NAME.sorted.bam


# ALIGNMENT- AND BAM-STATS
echo "[exome:merge_bam_and_bamstats $(date +"%Y-%m-%d %T")] Starting CollectAlignmentMetrics"
"$JAVA17_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
	CollectAlignmentSummaryMetrics \
        R=$GENOME \
        I=$LOCAL_NOVO/$SAMPLE_NAME.sorted.bam \
        O=$REPORTS_ALIGNSTATS_DIR/$SAMPLE_NAME.alignMetrics.txt

"$JAVA17_PATH"java -Xmx2g -jar "$PICARD_PATH"picard.jar \
	CalculateHsMetrics \
	BAIT_INTERVALS=$TARGET \
	TARGET_INTERVALS=$TARGET \
	I=$LOCAL_NOVO/$SAMPLE_NAME.sorted.bam \
	O=$REPORTS_BAMSTATS_DIR/$SAMPLE_NAME.bamMetrics.txt \
	REFERENCE_SEQUENCE=$GENOME


# BAM FILTERING
# sam to bam and filtering
samtools view -h -f 2 $LOCAL_NOVO/$SAMPLE_NAME.bam | awk 'substr($0,1,1) == "@" || $9 > 60 || $9 < -60 {print}' | samtools view -bS - > $LOCAL_NOVO/$SAMPLE_NAME.filt.bam
samtools sort $LOCAL_NOVO/$SAMPLE_NAME.filt.bam $LOCAL_NOVO/$SAMPLE_NAME.sort
samtools index $LOCAL_NOVO/$SAMPLE_NAME.sort.bam


# Local Realignment (OPTIONAL)
echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Running Local Realignment"

"$JAVA17_PATH"java -Xmx$MEMORY -jar "$GATK_PATH"GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R $GENOME \
   -I $LOCAL_NOVO/$SAMPLE_NAME.sort.bam \
   -o $LOCAL_NOVO/$SAMPLE_NAME.intervals \
   -nt $N_CORES
   
"$JAVA17_PATH"java -Xmx$MEMORY -jar "$GATK_PATH"GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R $GENOME \
   -I $LOCAL_NOVO/$SAMPLE_NAME.sort.bam \
   -targetIntervals $LOCAL_NOVO/$SAMPLE_NAME.intervals \
   -o $LOCAL_NOVO/$SAMPLE_NAME.realn.bam



# BAM ANNOTATION
"$JAVA18_PATH"java -Xmx$MEMORY -cp $CLASSPATH AnnotateBam -i $LOCAL_NOVO/$SAMPLE_NAME.realn.bam -o $BAM_DIR -s $SAMTOOLS_PATH -b $AMPLICONS -t 3

# MPILEUP
"$JAVA18_PATH"java -Xmx$MEMORY -cp $CLASSPATH AmpliStats -i $BAM_DIR/$SAMPLE_NAME.realn.annotated.bam -o $OUTPUT_MPILEUP -s $SAMTOOLS_PATH -r $GENOME -b $AMPLICONS

# Move output and remove temp files
mv $LOCAL_NOVO/$SAMPLE_NAME.novostats.txt $REPORTS_ALIGNMENT_DIR
rm -rf "$LOCAL_NOVO"

echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Alignment complete"



