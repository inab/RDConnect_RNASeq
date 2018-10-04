# First specify all references and tools path in the shell script
# for GATK - PICARD 
# dbSNP data from GATK [ see how to generate the files from GATK guidelines]



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


sh rd_connect_rna_call.sh \
            -b : my_bam.bam \
            -r : reference file \
            -s : sample id \
            -n : CPUS \
            -g : gigabyte memory \


