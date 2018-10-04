#Ejemplo de linia de comando para STAR:
set -e
# STAR 2.5.3a
# Rsem 1.3.0
# gcc 4.8.5

########################################################################################

READ1='/data/projects/cms/RNASeq/JohnDawson.NEW0961.R1.fastq.gz'
READ2='/data/projects/cms/RNASeq/JohnDawson.NEW0961.R2.fastq.gz'
OUTFOLDER='output'
PREFIX='aligned'
ID='SAMPLEID'
SM='SMID'
STARINDEX='/data/projects/RDconnect/analysis/alignmnent_quant/indexes'
RSEM_REFERENCE_NAME='/data/projects/RDconnect/analysis/alignmnent_quant/indexes/indexes'
TMPDIR='tmp'
# from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
REF='/data/resources/fasta/hg19/h37d5/hsapiens.hs37d5.fasta'
########################################################################################


  while getopts 'A:B:n:m:j:i:s:r:t:o:' flag; do
    case "${flag}" in
      A) READ1="${OPTARG}" ;;
      B) READ2="${OPTARG}" ;;
      n) CPU="${OPTARG}" ;;
      m) RAM="${OPTARG}g" ;;
      j) RSEM="${OPTARG}" ;;
      i) SAMPLEID="${OPTARG}" ;;
      s) STARINDEX="${OPTARG}" ;;
      r) REF="${OPTARG}" ;;
      t) TMPDIR="${OPTARG}";;
      o) OUTFOLDER="${OPTARG}";;
      *) error "Unexpected option ${flag}" ;;
    esac
  done





##########

mkdir -p $OUTFOLDER
mkdir -p $TMPDIR
ID=$SAMPLEID
SM=$SAMPLEID

#module load STAR/2.5.3a
#### my star is 2.1.0

STAR --runThreadN $CPU \
     --outSAMunmapped Within \
     --genomeDir $STARINDEX     \
     --readFilesIn $READ1 $READ2 \
     --outFileNamePrefix $OUTFOLDER/$SAMPLEID \
     --readFilesCommand zcat     \
     --quantMode TranscriptomeSAM \
     --outFilterType BySJout     \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8     \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --alignIntronMin 20     \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000     \
     --outTmpDir $TMPDIR/STAR_ \
     --outSAMtype BAM SortedByCoordinate     \
     --outSAMattrIHstart 0 \
     --outSAMattributes NH HI NM MD AS nM     \
     --outSAMattrRGline ID:$ID SM:$SM



rsem-calculate-expression -p 8 \
     --no-bam-output \
     --bam       \
     --temporary-folder $TMPDIR/rsem \
     --paired-end "$OUTFOLDER/$SAMPLEID"Aligned.toTranscriptome.out.bam  \
       $RSEM_REFERENCE_NAME \
       $OUTFOLDER/$PREFIX

