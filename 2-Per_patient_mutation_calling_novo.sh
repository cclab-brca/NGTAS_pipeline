#!/bin/bash 

# MUTATION CALLING PIPELINE FOR NGTAS DATA USING MUTECT2

#################################### INPUT INFO ####################################

# DIRECTORIES
WD=/path/to/working/directory
TEMP_DIR=$WD/temp
ANALYSIS_SUFFIX='.novo'
CALLER='mutect2' 
BAM_DIR=$WD/bam$ANALYSIS_SUFFIX 
CALLS_DIR=$WD/reports$ANALYSIS_SUFFIX/jointcalls
EXOME_DIRS=$TEMP_DIR:$CALLS_DIR:$CALLER_BAM_DIR

# Ensure all paths exist. Create directory structure if not present
ARR=$(echo $EXOME_DIRS | tr ":" "\n")
for x in $ARR
do
        if [ ! -d "$x" ]; then
        	mkdir -p $x
        fi

done

# REQUIRED FILES	
GENOME=/path/to/human_g1k_v37_decoy.fasta
GENOME_DIC=/path/to/human_g1k_v37_decoy.dict

TARGET_BED=/path/to/NGTAS_377amplicons_no_primers_0based.bed
TARGET_INT=/path/to/NGTAS_377amplicons_no_primers_0based_b37.intervals
MATCHES=/path/to/pairs_for_mutation_calling.txt
R_SCRIPTS=/path/to/Rscripts/folder

# SOFTWARE PATH
JAVA17_PATH=/path/to/jre1.7.0_67/bin/
JAVA18_PATH=/path/to/jdk1.8.0_31/bin/
GATK_PATH=/path/to/GenomeAnalysisTK-3.6/
BEDTOOLS_PATH=/path/to/bedtools2/bin
PICARD_PATH=/path/to/picard-tools-1.140/
BCF_PATH=/path/to/bcftools-1.2/
SAMTOOLS_PATH=/path/to/samtools-1.2/
BGZIP_PATH=/path/to/htslib-1.3/
ANNOVAR_PATH=/path/to/annovar/
PATH=$BEDTOOLS_PATH:$BCF_PATH:$SAMTOOLS_PATH:$BGZIP_PATH:$PATH

# FILTERING PARAMETERS AND OTHER INFO
vaf_th_norm=0.01	# min VAF in tumour/plasma
vaf_th_tum=0.01		# min VAF in buffy coat/normal sample
ratio=5			# min ratio between tumour VAF/normal VF
cov_th_norm=100		# min coverage in the normal
cov_th_tum=100		# min coverage in the tumour/plasma
th_ffpe=0.2		# min VAF to keep a C>T/A>G mutation in FFPE samples

N_CORES=8		# Number of cores available
MEMORY="8g"		# Available memory
####################################################################################

echo "[exome:mutation_calling $(date +"%Y-%m-%d %T")] Starting mutation calling"

# DEFINE VARIABLES 
GROUP_LIST=$(less $MATCHES | cut -f 6 | sort -u | xargs)
GROUP=$(echo $GROUP_LIST | cut -d ' ' -f "$LSB_JOBINDEX")	# NOTE: LSB_JOBINDEX variable can have different names depending on the HPC system used
SAMPLES=$(awk '$6 == "'"$GROUP"'" { print $5 }' $MATCHES | uniq)
NSAMPLES=$(echo $SAMPLES | wc -w)
NAMPLI=$(wc -l $COMMON_PATH/$TARGET_BED | cut -d ' ' -f 1)
ALL_T_FILENAMES=$(awk '$6 == "'"$GROUP"'" { print $1 }' $MATCHES)
ALL_N_FILENAMES=$(echo $(awk '$6 == "'"$GROUP"'" { print $3 }' $MATCHES) | xargs -n1 | sort -u)
ALL_FILENAMES1=$ALL_T_FILENAMES' '$ALL_N_FILENAMES
FILE_ARRAY1=($ALL_FILENAMES1)
from="%"
to=".realn.annotated.ba*"
ALL_FILENAMES2=${FILE_ARRAY1[@]/$from/$to}
FILE_ARRAY2=($ALL_FILENAMES2)
from="#"
to="$BAM_DIR/"
ALL_FILENAMES3=${FILE_ARRAY2[@]/$from/$to}

# move files
LOCAL="/tmp/mut_$RANDOM/"
mkdir "$LOCAL"
cp $GENOME  $LOCAL
cp $COMMON_PATH/$TARGET_BED  $LOCAL
cp $COMMON_PATH/$TARGET_INT  $LOCAL
cp $ALL_FILENAMES3 $LOCAL

# identify mutations in each amplicon
>$LOCAL$GROUP.called.txt
for i in $(seq 1 $NAMPLI)
do
	iAMPLI=$(sed -n "$i"p "$LOCAL"*.bed)
	AMPLINAME=$(echo $iAMPLI | cut -d ' ' -f 4)
	INTERVAL=$(echo $iAMPLI | cut -d ' ' -f 1):$(echo $iAMPLI | cut -d ' ' -f 2)-$(echo $iAMPLI | cut -d ' ' -f 3)
	# process all the samples from the same patient	
	for j in $(seq 1 $NSAMPLES)
	do
		SAMPLE=$(echo $SAMPLES | cut -d ' ' -f $j)
		TUMORFILES=$(awk '$5 == "'"$SAMPLE"'" && $6 == "'"$GROUP"'" { print $1 }' $MATCHES)
		TUMORNAMES=$(awk '$5 == "'"$SAMPLE"'" && $6 == "'"$GROUP"'" { print $2 }' $MATCHES)
		NORMALFILES=$(awk '$5 == "'"$SAMPLE"'" && $6 == "'"$GROUP"'" { print $3 }' $MATCHES)
		NORMALNAMES=$(awk '$5 == "'"$SAMPLE"'" && $6 == "'"$GROUP"'" { print $4 }' $MATCHES)
		SOURCE=$(awk '$5 == "'"$SAMPLE"'" && $6 == "'"$GROUP"'" { print $8 }' $MATCHES | uniq)
		# process all the replicates for the same sample
		for h in $(seq 1 3)
		do
			TUMORFILE=$(echo $TUMORFILES | cut -d ' ' -f $h)
			TUMORNAME=$(echo $TUMORNAMES | cut -d ' ' -f $h)
			NORMALFILE=$(echo $NORMALFILES | cut -d ' ' -f $h)
			NORMALNAME=$(echo $NORMALNAMES | cut -d ' ' -f $h)

			# generate single amplicon bam file
			samtools view -b -h -r $AMPLINAME "$LOCAL$TUMORFILE".realn.annotated.bam -o "$LOCAL$TUMORFILE".$AMPLINAME.bam
			samtools view -b -h -r $AMPLINAME "$LOCAL$NORMALFILE".realn.annotated.bam -o "$LOCAL$NORMALFILE".$AMPLINAME.bam

			# MODIFY BAM HEADER
			"$JAVA17_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
				AddOrReplaceReadGroups \
				VALIDATION_STRINGENCY="SILENT" \
				INPUT="$LOCAL$TUMORFILE".$AMPLINAME.bam \
				OUTPUT="$LOCAL$TUMORFILE".$AMPLINAME.rg.bam \
				RGID=$AMPLINAME \
				RGLB="library" \
				RGPL="platform" \
				RGPU="unit" \
				RGDS=$INTERVAL \
				RGSM=$TUMORNAME
			"$JAVA17_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
				AddOrReplaceReadGroups \
				VALIDATION_STRINGENCY="SILENT" \
				INPUT="$LOCAL$NORMALFILE".$AMPLINAME.bam \
				OUTPUT="$LOCAL$NORMALFILE".$AMPLINAME.rg.bam \
				RGID=$AMPLINAME \
				RGLB="library" \
				RGPL="platform" \
				RGPU="unit" \
				RGDS=$INTERVAL \
				RGSM=$NORMALNAME
			samtools index "$LOCAL$TUMORFILE".$AMPLINAME.rg.bam
			samtools index "$LOCAL$NORMALFILE".$AMPLINAME.rg.bam

			# MUTECT 2
			"$JAVA18_PATH"java -Xmx$MEMORY -jar "$GATK_PATH"GenomeAnalysisTK.jar \
				-T MuTect2 \
				-R $GENOME \
				-I:tumor $LOCAL$TUMORFILE.$AMPLINAME.rg.bam \
				-I:normal $LOCAL$NORMALFILE.$AMPLINAME.rg.bam \
				-L $LOCAL$TARGET_INT \
				-stand_emit_conf 0.0 \
				--minPruning 5 \
				-o $LOCAL$TUMORNAME.$NORMALNAME.$AMPLINAME.mut2.vcf \
				-nct 1
			# filtering
			R < $R_SCRIPTS/filter_mutect2.R $LOCAL $TUMORNAME.$NORMALNAME.$AMPLINAME $vaf_th_norm $ratio $cov_th_norm $cov_th_tum --no-save
		done
		# select mutations called in 2 out of 3 replicates
		R1=$LOCAL$(echo $TUMORNAMES | cut -d ' ' -f 1).$(echo $NORMALNAMES | cut -d ' ' -f 1).$AMPLINAME.mut2.filt.vcf
		R2=$LOCAL$(echo $TUMORNAMES | cut -d ' ' -f 2).$(echo $NORMALNAMES | cut -d ' ' -f 2).$AMPLINAME.mut2.filt.vcf
		R3=$LOCAL$(echo $TUMORNAMES | cut -d ' ' -f 3).$(echo $NORMALNAMES | cut -d ' ' -f 3).$AMPLINAME.mut2.filt.vcf
		if [ -s $R1 ] && [ -s $R2 ]
		then
			bedtools intersect -u -a $R1 -b $R2 -header | bgzip -cf > $LOCAL$SAMPLE.$AMPLINAME.R1R2.vcf.gz
			bcftools index -tf $LOCAL$SAMPLE.$AMPLINAME.R1R2.vcf.gz
		fi

		if [ -s $R1 ] && [ -s $R3 ]
		then		
			bedtools intersect -u -a $R1 -b $R3 -header | bgzip -cf > $LOCAL$SAMPLE.$AMPLINAME.R1R3.vcf.gz
			bcftools index -tf $LOCAL$SAMPLE.$AMPLINAME.R1R3.vcf.gz
		fi

		if [ -s $R2 ] && [ -s $R3 ]
		then		
			bedtools intersect -u -a $R2 -b $R3 -header | bgzip -cf > $LOCAL$SAMPLE.$AMPLINAME.R2R3.vcf.gz
			bcftools index -tf $LOCAL$SAMPLE.$AMPLINAME.R2R3.vcf.gz
		fi
		if [ -s $LOCAL$SAMPLE.$AMPLINAME.R1R2.vcf.gz ] || [ -s $LOCAL$SAMPLE.$AMPLINAME.R1R3.vcf.gz ] || [ -s $LOCAL$SAMPLE.$AMPLINAME.R2R3.vcf.gz ]
		then
			bcftools concat -aD $LOCAL$SAMPLE.$AMPLINAME.R*R*.vcf.gz -o $LOCAL$SAMPLE.$AMPLINAME.vcf
		fi

		if [ -s $LOCAL$SAMPLE.$AMPLINAME.vcf ]
		then
			# FILTERING PER SAMPLE
			# remove FILTER information
			bcftools annotate -x "FILTER" $LOCAL$SAMPLE.$AMPLINAME.vcf -o $LOCAL$SAMPLE.$AMPLINAME.nofilt.vcf

			# sort VCF (contigs)
			"$JAVA17_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
			SortVcf \
			I=$LOCAL$SAMPLE.$AMPLINAME.nofilt.vcf \
			O=$LOCAL$SAMPLE.$AMPLINAME.sort.vcf \
			SEQUENCE_DICTIONARY=$GENOME_DIC

			# bam list
			SAMPLEFILES=($(echo $TUMORFILES) $(echo $NORMALFILES))
			from="%"
			to=.$AMPLINAME.rg.bam
			BAM_ARRAY=(${SAMPLEFILES[@]/$from/$to})
			from="#"
			to="-I $LOCAL"
			BAM_LIST=${BAM_ARRAY[@]/$from/$to}

			# Haplotypecaller
			"$JAVA18_PATH"java -Xmx$MEMORY -jar "$GATK_PATH"GenomeAnalysisTK.jar \
				-T HaplotypeCaller \
				-R $GENOME \
				$BAM_LIST \
				-gt_mode GENOTYPE_GIVEN_ALLELES \
				-alleles $LOCAL$SAMPLE.$AMPLINAME.sort.vcf \
				-L $LOCAL$SAMPLE.$AMPLINAME.sort.vcf \
				-U ALLOW_SEQ_DICT_INCOMPATIBILITY \
				--downsampling_type NONE \
				--minPruning 2 \
				--min_base_quality_score 20 \
				-stand_emit_conf 0.0 \
				-stand_call_conf 0.0 \
				-nct $N_CORES \
				-out_mode EMIT_ALL_SITES \
				-o $LOCAL$SAMPLE.$AMPLINAME.hc.vcf

			# filtering in R
			R < $R_SCRIPTS/filter_persample.R $LOCAL $SAMPLE.$AMPLINAME $TUMORNAMES $NORMALNAMES $SOURCE $vaf_th_norm $vaf_th_tum $ratio  $th_ffpe --no-save
			if [ -s $LOCAL$SAMPLE.$AMPLINAME.called.txt ]
			then 
				cat $LOCAL$SAMPLE.$AMPLINAME.called.txt >> $LOCAL$GROUP.called.txt
			fi
		fi
		if [ -s $LOCAL$SAMPLE.$AMPLINAME.sort.filt.vcf ]
		then
			bgzip -cf $LOCAL$SAMPLE.$AMPLINAME.sort.filt.vcf > $LOCAL$SAMPLE.$AMPLINAME.filt.vcf.gz
			bcftools index -tf $LOCAL$SAMPLE.$AMPLINAME.filt.vcf.gz
		fi
	done

	nfiles=$(ls $LOCAL*$AMPLINAME.filt.vcf.gz | wc -l)
	if [ $nfiles -gt 0 ]
	then
		# combine calls from all samples
		bcftools concat -aD $LOCAL*$AMPLINAME.filt.vcf.gz -o $LOCAL$GROUP.$AMPLINAME.vcf

		# remove FILTER information
		bcftools annotate -x "FILTER" $LOCAL$GROUP.$AMPLINAME.vcf -o $LOCAL$GROUP.$AMPLINAME.nofilt.vcf

		# sort VCF (contigs)
		"$JAVA17_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
		SortVcf \
		I=$LOCAL$GROUP.$AMPLINAME.nofilt.vcf \
		O=$LOCAL$GROUP.$AMPLINAME.sort.vcf \
		SEQUENCE_DICTIONARY=$GENOME_DIC
		# compress for downstream concat
		bgzip -cf $LOCAL$GROUP.$AMPLINAME.sort.vcf > $LOCAL$GROUP.$AMPLINAME.sort.vcf.gz
		bcftools index -tf $LOCAL$GROUP.$AMPLINAME.sort.vcf.gz

		# bam list
		BAM_ARRAY=($LOCAL*.$AMPLINAME.rg.bam)
		from="#"
		to="-I "
		BAM_LIST=${BAM_ARRAY[@]/$from/$to}

		# Haplotypecaller
		"$JAVA18_PATH"java -Xmx$MEMORY -jar "$GATK_PATH"GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-R $GENOME \
			$BAM_LIST \
			-gt_mode GENOTYPE_GIVEN_ALLELES \
			-alleles $LOCAL$GROUP.$AMPLINAME.sort.vcf \
			-L $LOCAL$GROUP.$AMPLINAME.sort.vcf \
			-U ALLOW_SEQ_DICT_INCOMPATIBILITY \
			--downsampling_type NONE \
			--minPruning 2 \
			--min_base_quality_score 20 \
			-stand_emit_conf 0.0 \
			-stand_call_conf 0.0 \
			-nct $N_CORES \
			-out_mode EMIT_ALL_SITES \
			-o $LOCAL$GROUP.$AMPLINAME.hc.vcf

		# Annovar
		LOCAL_TEMP="$LOCAL"temp
		mkdir $LOCAL_TEMP
		perl "$ANNOVAR_PATH"table_annovar.pl $LOCAL$GROUP.$AMPLINAME.hc.vcf \
		"$ANNOVAR_PATH"humandb/ \
		-buildver hg19 \
		-out $LOCAL$GROUP.$AMPLINAME \
		-remove \
		-protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,snp138,ljb26_all \
		-operation g,r,r,f,f,f,f \
		-nastring . \
		--onetranscript \
		--tempdir $LOCAL_TEMP \
		-vcfinput

		# add amplicon info and bgzip
		R < $R_SCRIPTS/Add_amplicon_info_perpatient_calls.R $LOCAL$GROUP.$AMPLINAME.hg19_multianno.vcf --no-save
		bgzip -cf $LOCAL$GROUP.$AMPLINAME.hg19_multianno.vcf.ampli > $LOCAL$GROUP.$AMPLINAME.hg19_multianno.vcf.gz
		bcftools index -tf $LOCAL$GROUP.$AMPLINAME.hg19_multianno.vcf.gz
	fi

	rm $LOCAL*mut2*
	rm $LOCAL*.hc.vcf*
	rm $LOCAL*.sort.vcf
	rm $LOCAL*.sort.filt.vcf
	rm $LOCAL*.sort.vcf.idx
	rm $LOCAL*$AMPLINAME.vcf
	rm $LOCAL*$AMPLINAME.nofilt.vcf
	rm $LOCAL*$AMPLINAME.filt.vcf.gz*
	rm $LOCAL*$AMPLINAME.R*R*.vcf.gz*
	rm $LOCAL*.hg19_multianno.vcf
	rm $LOCAL*.hg19_multianno.vcf.ampli
	rm $LOCAL*$AMPLINAME.bam
	rm $LOCAL*$AMPLINAME.rg.bam*
	rm $LOCAL*$AMPLINAME.called.txt
done

# merge calls from all amplicons
TBI_ARRAY=($LOCAL*.hg19_multianno.vcf.gz.tbi)
from=".tbi"
to=""
FILE_LIST=${TBI_ARRAY[@]/$from/$to}
# HC output
bcftools concat -a $FILE_LIST -o $CALLS_DIR/$GROUP.hc.annot.vcf
# mutation list vcf
bcftools concat -a $LOCAL$GROUP*.sort.vcf.gz -o $CALLS_DIR/$GROUP.mutlist.vcf
# mutation list txt
cp $LOCAL$GROUP.called.txt $CALLS_DIR

rm -r $LOCAL

