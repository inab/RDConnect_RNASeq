/*
 *
 *   This file takes care to call a VCF call + BAM files and ASE values for RNA seq data already aligned.
 */
 
/* 
 * 
 * Authors:
 * - Mattia Bosio <mattia-.bosio85@gmail.com>
 */ 

 
/*
 * Default pipeline parameters. They can be overriden on the command line eg. 
 * given `params.foo` specify on the run command line `--foo some_value`.  
*/



OUTPREFIX = params.OUTPREFIX
log.info """
         
                   P I P E L I N E    
         ===================================
         ## Input data and prefix
         Input BAM file:  ${params.bam}
         Input BAI file:  ${params.bai}
         Output folder :  ${params.outdir}
         Sample ID     :  ${params.SAMPLEID}
         ...
         
         ## Reference resuources
         Genome fasta :  ${params.ref}
         Fasta index  :  ${params.ref_fai}
         Faste Dict   :  ${params.ref_dict}
         
         ## Paths to tools
         GATK             : ${params.GATK}
         RNASeqCaller:  ${params.code_path}  
         
         """
         .stripIndent()

code_path = params.code_path

bam_file_ch  = Channel.fromPath(params.bam)
bai_file_ch  = Channel.fromPath(params.bai)


ref_file_ch  = Channel.fromPath(params.ref)
ref_file_ch.into { ref_file_1; ref_file_2;ref_file_3; ref_file_4 ; ref_file_5; ref_file_6; ref_file_7}

ref_fai_ch  = Channel.fromPath(params.ref_fai)
ref_fai_ch.into {ref_fai_1; ref_fai_2 ; ref_fai_3; ref_fai_4; ref_fai_5; ref_fai_6; ref_fai_7}

ref_dict_ch  = Channel.fromPath(params.ref_dict)
ref_dict_ch.into {ref_dict_1; ref_dict_2 ; ref_dict_3 ; ref_dict_4; ref_dict_5; ref_dict_6; ref_dict_7}

dbindel_file = Channel.fromPath(params.DBINDEL)
dbsnp_file = Channel.fromPath(params.DBSNP)


/*
*  Step 1 : AddReplace read group"
*  Step 2 : Duplicate marking
*  Step 3 : Split n trim
*  Step 4 : Indel realingment
*  Step 5 : BQSR
*  Step 6 : Variant Calling
*  Step 7 : Compression
*  Step 8 : Variant filtering
*  Step 9 : ASE Read counter
*/
 
process normalize {
    tag "AddReplaceRG"
    label 'picard'
    input:
    file bam_in from bam_file_ch
    file bai_in from bai_file_ch
    file fasta_ref from ref_file_1
    file fai_ref from ref_fai_1
    file dict_ref from ref_dict_1
     
    output:
    file("Step1_out.bam") into step1_bam
 
    script:
    """
    java  -jar $params.PICARD  AddOrReplaceReadGroups \
     I=$bam_in \
     O=Step1_out.bam \
     SO=coordinate \
     RGID=$params.SAMPLEID \
     RGLB=$params.RGLB \
     RGPL=$params.RGPL \
     RGPU=$params.RGPU \
     RGSM=$params.SAMPLEID \
     USE_JDK_DEFLATER=true \
     USE_JDK_INFLATER=true
    """ 
}




process Duplicate_marking {
    tag "Duplicate marking"
				label 'picard'
    input:
    file bam1_in from step1_bam
				
    output:
    file( "Step2_out.bam") into step2_bam
				file( "Step2_out.bai") into step2_bai

    script:
    gbmem = "${task.memory.toGiga()}g"
    """
    java -Xmx${gbmem} -jar $params.PICARD MarkDuplicates \
     I=$bam1_in \
     O=Step2_out.bam  \
     CREATE_INDEX=true \
     VALIDATION_STRINGENCY=SILENT \
     M=output.metrics    \
     USE_JDK_DEFLATER=true \
     USE_JDK_INFLATER=true
    """  
}  
  
process SplitNTrim {
    tag "Split n trim"
				label 'GATK'
    input:
    file bam2_in from step2_bam
    file bai2_in from step2_bai
    file fasta_ref2 from ref_file_2
    file fai_ref2 from ref_fai_2
    file dict_ref2 from ref_dict_2
				
    output:
    file( "Step3_out.bam") into step3_bam
				file( "Step3_out.bai") into step3_bai


    script:
    """
       java  -jar $params.GATK \
     -T SplitNCigarReads \
     -R $fasta_ref2 \
     -I $bam2_in \
     -o Step3_out.bam \
     -rf ReassignOneMappingQuality \
     -RMQF 255 \
     -RMQT 60 \
     -U ALLOW_N_CIGAR_READS 
    """
    
}
  
////////////////////////////////////////////////////////
process IndelRealign {
    tag "Indel realingment"
				label 'GATK'
    input:
    file bam3_in from step3_bam
    file bai3_in from step3_bai
    file fasta_ref3 from ref_file_3
    file fai_ref3 from ref_fai_3
    file dict_ref3 from ref_dict_3
    file dbindel from dbindel_file
				
    output:
				file( "Step4_out.bam") into step4_bam


    script:
    """
       java -jar $params.GATK \
      -T RealignerTargetCreator \
      -R $fasta_ref3 \
      -I $bam3_in \
      -o Step4.intervals \
      -known $dbindel \
      --minReadsAtLocus 6 \
      --maxIntervalSize 200 \
      --downsampling_type NONE  
   
   #   -U ALLOW_SEQ_DICT_INCOMPATIBILITY
      
    java -jar $params.GATK \
      -T IndelRealigner \
      -R $fasta_ref3\
      -I $bam3_in \
      -known $dbindel \
      -targetIntervals Step4.intervals \
      --maxReadsForRealignment 10000 \
      --consensusDeterminationModel USE_SW \
      --downsampling_type NONE \
      -o Step4_out.bam
    """
  }
////////////////////////////////////////////////////////
  process BQSR {
    tag "BQSR"
				label 'BQSR'
    publishDir params.outdir, mode:'copy'

    input:
    file bam4_in from step4_bam
    file fasta_ref4 from ref_file_4
    file fai_ref4 from ref_fai_4
    file dict_ref4 from ref_dict_4
    file dbsnp from dbsnp_file
				
    output:
				file( "Step5_out.bam") into step5_bam

    script:
    """
     java -jar $params.GATK \
      -T BaseRecalibrator \
      -nct $task.cpus \
      --default_platform illumina \
      -cov ReadGroupCovariate \
      -cov QualityScoreCovariate \
      -cov CycleCovariate \
      -knownSites $dbsnp \
      -cov ContextCovariate  \
      -R $fasta_ref4 \
      -I $bam4_in \
      --downsampling_type NONE \
      -o recalibration_report.grp
     
     java -jar $params.GATK \
     -T PrintReads \
     -R $fasta_ref4 \
     -I $bam4_in \
     -BQSR recalibration_report.grp \
     -o Step5_out.bam   

    """

     

  }
 
step5_bam.into{step5_bam_1;step5_bam_2;step5_bam_3 }

 
////////////////////////////////////////////////////////
 process VarCalling {
    tag "VarCalling"
				label 'VarCalling'
    
    input:
    file bam5_in from step5_bam_1
    file fasta_ref5 from ref_file_5
    file fai_ref5 from ref_fai_5
    file dict_ref5 from ref_dict_5
				
    output:
				file( "Step6_out.g.vcf") into step6_vcf

    script:    
    """
    gbmem = "${task.memory.toGiga()}g"
    java -Xmx${gbmem} -jar $params.GATK \
     -T HaplotypeCaller \
        -R $fasta_ref5 \
        -I $bam5_in \
         -dontUseSoftClippedBases \
         -o Step6_out.g.vcf \
         --pair_hmm_implementation LOGLESS_CACHING \
         -nct $task.cpus \
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

    """
    
 }
 
 ////////////////////////////////////////////////////////
 process Compression {
    tag "Compression"
				
    publishDir params.outdir, mode:'copy'
    input:
    file vcf_in from step6_vcf
    				
    output:
				file( "Step6_out.g.vcf.gz") into step7_vcf_gz
    file( "Step6_out.g.vcf.gz.tbi") into step7_vcf_gz_tbi
    
    script:
    """
      bgzip $vcf_in
      tabix $vcf_in".gz"
    """
 }
 
 ////////////////////////////////////////////////////////
 process Filter {
  tag "Filter"
  label 'GATK'
  publishDir params.outdir, mode:'copy'
    input:
    file gzvcf_in from step7_vcf_gz
    file gzvcf_tbi_in from step7_vcf_gz_tbi
    file fasta_ref6 from ref_file_6
    file fai_ref6 from ref_fai_6
    file dict_ref6 from ref_dict_6
    
    output:
    file "Step8_out.filtered.g.vcf" into step8_vcf
  
   script:
   """
      java -jar $params.GATK \
       -T VariantFiltration \
       -R $fasta_ref6 \
       -V $gzvcf_in \
       -window 35 -cluster 3 \
       -filterName FS -filter "FS > 30.0" \
       -filterName QD -filter "QD < 2.0" \
       -o Step8_out.filtered.g.vcf
  
   """
 
 }
 
 ////////////////////////////////////////////////////////
 process ASE {
  tag " ASE Read counter"
  label 'GATK'
  
  publishDir params.outdir, mode:'copy'
    input:
       file vcf_filtered from step8_vcf
       file bam5_in_2  from step5_bam_2
       file fasta_ref7 from ref_file_7
       file fai_ref7   from ref_fai_7
       file dict_ref7  from ref_dict_7
    output:
      
      
    script:
    """
    java -jar $params.GATK \
      -T GenotypeGVCFs \
      -R $fasta_ref7 \
      --variant $vcf_filtered \
      -o Step9_input.vcf

    java -jar $params.GATK \
     -T ASEReadCounter \
     -R $fasta_ref7 \
     -I $bam5_in_2 \
     -sites Step9_input.vcf \
     -U ALLOW_N_CIGAR_READS \
     -minDepth 10 \
     --minBaseQuality 20 \
     --minMappingQuality 30 \
     --includeDeletions \
     --outputFormat CSV \
     -o $params.SAMPLEID"_ASEReadCounter.table"
    """
 }
 
 /* num = Channel.from(params.chromosomes.splitCsv())*/

workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report\n" : "Oops .. something went wrong" )

} 
