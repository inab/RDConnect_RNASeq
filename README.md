# RDConnect_RNASeq

## Step 1 : Alignment
From fastq files for paired end files will run 
-  STAR 2.5.3a to generate alignment to genome and transcriptome
-  RSEM xxxxx to quantify gene and isoform usage 

It requires to have previously generated STAR indexes apt for the task.
Same goes for RSEM indexes. Both folders are to be added as input for the call
We require also that the reference fasta file is processed as by GATK requirements (with samtools faidx and Picard)

sh rd_connect_align.sh \
         -A f1.fq.gz \
         -B f2.fq.gz \
         -n CPUS \
         -m gigabyte memory \
         -j RSEM_indexes_folder \
         -i sampleID \
         -s STAR_indexes_folder \
         -r reference_fasta_file \
         -t tmp_folder \
         -o results_folder

## Step 2: Variant Calling
- This step starts from the aligned bam file to the genome and produces
- One processed bam file (with adjusted scores, duplicate removed, etc.)
- One gVCF file with all variants and non variant sites
- One VCF file with all filtered variants [these  are the variants for RD-Connect site]
- It also produces an Allele Specific Report for each variant, 

sh rd_connect_rna_call.sh \
            -b : my_bam.bam \
            -r : reference file \
            -s : sample id \
            -n : CPUS \
            -g : gigabyte memory \


## Requirements
- GATK 3.6.0
- STAR 2.5.3a
- RSEM xxxx
- Picard xxxx
- Tabix

## Nextflow execution / Docker
Both steps come with nextflow implementation (in Nextflow folder) which allow for full reproducibility, counting with a Docker implementation using publicly available images.

