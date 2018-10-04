/*
 *
 *   This file takes care to align paired end FQ files.
 */
 
/* 
 * 
 * Authors:
 * - Mattia Bosio <mattia-.bosio@bsc.es>
 */ 

 
/*
 * Default pipeline parameters. They can be overriden on the command line eg. 
 * given `params.foo` specify on the run command line `--foo some_value`.  
*/


/ my inputs - outputs/

OUTPREFIX = params.OUTPREFIX
log.info """
         
                   P I P E L I N E    
         ===================================
         ## Input data and prefix
         Input READ1 file:  ${params.fq1}
         Input READ1 file:  ${params.fq2}
         Output folder :  ${params.outdir}
         Sample ID     :  ${params.SAMPLEID}
         
         ## Reference resuources
         Genome fasta :  ${params.ref}
         Fasta index  :  ${params.ref_fai}
         Fasta Dict   :  ${params.ref_dict}
         STAR_index   :  ${params.STARINDEX}
         RSEM_index   :  ${params.RSEMINDEX}
       
        
         """
         .stripIndent()


fq1_ch  = Channel.fromPath(params.fq1)
fq2_ch  = Channel.fromPath(params.fq2)


ref_file_ch  = Channel.fromPath(params.ref)
//ref_file_ch.into { ref_file_1; ref_file_2;ref_file_3; ref_file_4 ; ref_file_5; ref_file_6; ref_file_7}

ref_fai_ch  = Channel.fromPath(params.ref_fai)
//ref_fai_ch.into {ref_fai_1; ref_fai_2 ; ref_fai_3; ref_fai_4; ref_fai_5; ref_fai_6; ref_fai_7}

ref_dict_ch  = Channel.fromPath(params.ref_dict)
// ref_dict_ch.into {ref_dict_1; ref_dict_2 ; ref_dict_3 ; ref_dict_4; ref_dict_5; ref_dict_6; ref_dict_7}

idxes_star = Channel.fromPath(params.STARINDEX +"/*")
idxes_rsem = Channel.fromPath(params.RSEMINDEX +"/*")


 
process Align {
    tag "Align"
    label 'STAR'
    
    publishDir params.outdir, mode:'copy'
    
    input:
    file fq1 from fq1_ch
    file fq2 from fq2_ch
    file fasta_ref from ref_file_ch
    file fai_ref from ref_fai_ch
    file dict_ref from ref_dict_ch
    file all_indexes from idxes_star.collect() 
     
    output:
    file(params.SAMPLEID + "Aligned.toTranscriptome.out.bam") into transcriptome_bam
    file(params.SAMPLEID + "Aligned.sortedByCoord.out.bam") into aligned_bam
    script:
    """
    STAR --runThreadN $task.cpus \
     --outSAMunmapped Within \
     --genomeDir ./    \
     --readFilesIn $fq1 $fq2 \
     --outFileNamePrefix $params.SAMPLEID \
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
     --outTmpDir ./STAR_ \
     --outSAMtype BAM SortedByCoordinate     \
     --outSAMattrIHstart 0 \
     --outSAMattributes NH HI NM MD AS nM     \
     --outSAMattrRGline ID:$params.SAMPLEID SM:$params.SAMPLEID

    """ 
}




process RSEM {
    tag "RSEM"
				label 'RSEM'
    
    publishDir params.outdir, mode:'copy'
    
    input:
    file bam1_in from aligned_bam
				file all_indexes from idxes_rsem.collect()
    
    output:
    file( params.SAMPLEID + ".isoforms.results") into isoforms_file
				file( params.SAMPLEID +".genes.results") into genes_file


    script:
    """
				 rsem-calculate-expression \
       -p $task.cpus \
       --no-bam-output \
       --bam       \
       --temporary-folder ./rsem \
       --paired-end $bam1_in \
         ./ \
         $params.$SAMPLEID

    """  
}  
 

workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report\n" : "Oops .. something went wrong" )

} 
