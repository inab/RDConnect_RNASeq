set -e

####### SAMPLE VARIANTS #####
ID='testID'
SO='coordinate'
RGID=$ID
RGLB='library'
RGPL='platform'
RGPU='machine'
RGSM=$ID

####### TOOLS VALUES ########
GATK='/home/mbosio/software/GATK/GATK_3.6.0/'
PICARD='/home/mbosio/software/Picard/v2.17.10/picard.jar'
DBSNP='/home/mbosio/software/GATK/GATK_3.6.0/resources/dbsnp_138.hg19.snps.vcf'
DBINDEL='/home/mbosio/software/GATK/GATK_3.6.0/resources/dbsnp_138.hg19.indels.vcf'
REF='hs37d5.fa'
CODE='/data/projects/RDconnect/code/'

###### EXECUTION VAL. #######
Execution_dir="tmp"
Current_dir=$(pwd)
RAM='20g'
CPU='1'
BAM_IN='input.bam'
TMPDIR='/data/projects/RDconnect/analysis/javatmpdir/'
SEX='m'

export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR $JAVA_TOOL_OPTIONS"


###### COMMANDLINE ARGS  ######

  ## Set some flags
  #  -b : bam file
  #  -r : reference file
  #  -s : sample id
  #  -n : CPUS
  #  -g : gigabyte memory



  while getopts 'b:r:s:n:g:t:f' flag; do
    case "${flag}" in
      b) BAM_IN="${OPTARG}" ;;
      s) ID="${OPTARG}" ;;
      n) CPU="${OPTARG}" ;;
      g) RAM="${OPTARG}g" ;;
      t) Execution_dir="${OPTARG}" ;;
      *) error "Unexpected option ${flag}" ;;
    esac
  done

  RGID=$ID
  RGSM=$ID
  RGSM=$(echo "${ID%_*}")

  

mkdir -p $Execution_dir
cd $Execution_dir


#Step 1 Add replace read group
  echo "-------------------------------------------------"
  echo "Step 1 AddReplace read group"
  echo "-------------------------------------------------"

if  [  -s Step1_done.txt  ] ; then
      echo 'Step1_out.bam already exists'
else
  time java -Xmx$RAM -jar $PICARD  AddOrReplaceReadGroups \
     I=$BAM_IN \
     O=Step1_out.bam \
     SO=coordinate \
     RGID=$RGID \
     RGLB=$RGLB \
     RGPL=$RGPL \
     RGPU=$RGPU \
     RGSM=$RGSM

 echo 'done' > Step1_done.txt
fi



#Step 2  Mark duplicates
  echo "-------------------------------------------------"
  echo "Step 2 Duplicate marking"
  echo "-------------------------------------------------"

if  [  -s Step2_done.txt  ] ; then
      echo 'Step2_out.bai already exists'
else  
  time java -Xmx$RAM -jar $PICARD MarkDuplicates \
     I=Step1_out.bam \
     O=Step2_out.bam  \
     CREATE_INDEX=true \
     VALIDATION_STRINGENCY=SILENT \
     M=output.metrics \
     use_jdk_deflater=true \
     use_jdk_inflater=true
 
 echo 'done' > Step2_done.txt

fi 

#Step3  Split and trim 
  echo "-------------------------------------------------"
  echo "Step 3 Split n trim"
  echo "-------------------------------------------------"


if  [  -s Step3_done.txt  ] ; then         
      echo 'Step3_out.bam already exists'
else
  time java  -Xmx$RAM -jar $GATK/GenomeAnalysisTK.jar \
     -T SplitNCigarReads \
     -R $REF \
     -I Step2_out.bam \
     -o Step3_out.bam \
     -rf ReassignOneMappingQuality \
     -RMQF 255 \
     -RMQT 60 \
     -U ALLOW_N_CIGAR_READS 

 echo 'done' > Step3_done.txt


#     -U ALLOW_N_CIGAR_READS
fi

#Step 4  Indel realingment [Marginal benefit]
  echo "-------------------------------------------------"
  echo "Step 4 Indel realingment"
  echo "-------------------------------------------------"

if  [  -s Step4_done.txt  ] ; then         
      echo 'Step4_out.bam already exists'
else

 time java -Xmx$RAM -jar $GATK/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R $REF \
   -I Step3_out.bam \
   -o Step4.intervals \
   -known $DBINDEL \
   --minReadsAtLocus 6 \
   --maxIntervalSize 200 \
   --downsampling_type NONE  

#   -U ALLOW_SEQ_DICT_INCOMPATIBILITY


 time java -Xmx$RAM -jar $GATK/GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R $REF \
   -I Step3_out.bam \
   -known $DBINDEL \
   -targetIntervals Step4.intervals \
   --maxReadsForRealignment 10000 \
   --consensusDeterminationModel USE_SW \
   --downsampling_type NONE \
   -o Step4_out.bam 

 echo 'done' > Step4_done.txt

fi


#Step5 BQSR   [Marginal benefit]
  echo "-------------------------------------------------"
  echo "Step 5 BQSR "
  echo "-------------------------------------------------"

if  [  -s Step5_done.txt  ] ; then         
      echo 'Step5_out.bam already exists'
else
  time java -Xmx$RAM -jar $GATK/GenomeAnalysisTK.jar \
     -T BaseRecalibrator \
     -nct 10 \
     --default_platform illumina \
     -cov ReadGroupCovariate \
     -cov QualityScoreCovariate \
     -cov CycleCovariate \
     -knownSites $DBSNP \
     -cov ContextCovariate -R $REF \
     -I Step4_out.bam \
     --downsampling_type NONE \
     -o recalibration_report.grp 

   
  

  time java -Xmx$RAM -jar $GATK/GenomeAnalysisTK.jar \
   -T PrintReads \
   -R $REF \
   -I Step4_out.bam\
   -BQSR recalibration_report.grp \
   -o Step5_out.bam   

echo 'done' > Step5_done.txt

fi 


  echo "-------------------------------------------------"
  echo "Step 6  Variant Calling"
  echo "-------------------------------------------------"

if  [  -s Step6_done.txt  ] ; then
      echo 'Step6   already done'
else

  time java -Xmx$RAM  -jar $GATK/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
        -R $REF \
        -I Step5_out.bam \
         -dontUseSoftClippedBases \
         -o Step6_out.g.vcf \
         --pair_hmm_implementation LOGLESS_CACHING \
         -nct 10 \
         -dt NONE \
         -rf BadCigar \
         --never_trim_vcf_format_field \
         -A AlleleBalance \
         -A BaseCounts \
         -A BaseQualityRankSumTest \
         -A ChromosomeCounts \
         -A ClippingRankSumTest \
         -A Coverage \
         -A DepthPerAlleleBySample \
         -A DepthPerSampleHC \
         -A FisherStrand \
         -A GCContent \
         -A HaplotypeScore \
         -A HardyWeinberg \
         -A HomopolymerRun \
         -A ClippingRankSumTest \
         -A LikelihoodRankSumTest \
         -A LowMQ \
         -A MappingQualityRankSumTest \
         -A MappingQualityZero \
         -A MappingQualityZeroBySample \
         -A NBaseCount \
         -A QualByDepth \
         -A RMSMappingQuality \
         -A ReadPosRankSumTest \
         -A StrandBiasBySample \
         -A StrandOddsRatio \
         -A VariantType \
         -ploidy 2 \
         --min_base_quality_score 10 \
         -ERC GVCF \
         -variant_index_type LINEAR \
         -variant_index_parameter 128000 \
         --GVCFGQBands 20 \
         --GVCFGQBands 25 \
         --GVCFGQBands 30 \
         --GVCFGQBands 35 \
         --GVCFGQBands 40 \
         --GVCFGQBands 45 \
         --GVCFGQBands 50 \
         --GVCFGQBands 70 \
         --GVCFGQBands 90 \
         --GVCFGQBands 99 \
         --standard_min_confidence_threshold_for_calling 30 \
         --standard_min_confidence_threshold_for_emitting 10 

       echo 'Done' > Step6_done.txt

fi


  echo "-------------------------------------------------"
  echo "Step 7 Compression"
  echo "-------------------------------------------------"
    if  [  -s $ID"_out.g.vcf.gz"  ] ; then         
      echo 'Step 7 already done'
  else
      bgzip Step6_out.g.vcf 
      mv Step6_out.g.vcf.gz $ID"_out.g.vcf.gz"
      tabix -p vcf $ID"_out.g.vcf.gz"
   fi 



  echo "-------------------------------------------------"
  echo "Step 8 Variant filtering"
  echo "-------------------------------------------------"

if  [  -s Step8_done.txt   ] ; then         
      echo 'Step8_out.filtered.vcf  already exists'
else
 echo 'doing step 8'
 gzip -d  $ID"_out.g.vcf.gz"
 bgzip $ID"_out.g.vcf"
 tabix -p vcf $ID"_out.g.vcf.gz"
 
 time java -Xmx$RAM -jar $GATK/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $REF \
  -V $ID"_out.g.vcf.gz" \
  -window 35 -cluster 3 \
  -filterName FS -filter "FS > 30.0" \
  -filterName QD -filter "QD < 2.0" \
  -o Step8_out.filtered.g.vcf
 
  echo 'Done' > Step8_done.txt
fi

  echo "-------------------------------------------------"
  echo "Step 9 ASE Read counter"
  echo "-------------------------------------------------"
if  [  -s Step9_done.txt   ] ; then
      echo 'Step9   already completed.'
else

 java -Xmx$RAM -jar $GATK/GenomeAnalysisTK.jar \
   -T GenotypeGVCFs \
   -R $REF \
   --variant Step8_out.filtered.g.vcf \
   -o Step9_input.vcf


  time java -Xmx$RAM -jar $GATK/GenomeAnalysisTK.jar \
    -T ASEReadCounter \
    -R $REF \
    -I Step5_out.bam \
    -sites Step9_input.vcf \
    -U ALLOW_N_CIGAR_READS \
    -minDepth 10 \
    --minBaseQuality 20 \
    --minMappingQuality 30 \
    --includeDeletions \
    --outputFormat CSV \
    -o $ID"_ASEReadCounter.table"

  # now calc pvalue per variant if needed 
  # python $CODE/ASE_pval.py $ID"_ASEReadCounter.table" > $ID"_ASEReadCounter_pval.csv" 
  echo 'Done' > Step9_done.txt
fi


  echo "-------------------------------------------------"
  echo "Step 10 Cleanup"
  echo "-------------------------------------------------"
  
  #remove internediate products [bam files mostly]