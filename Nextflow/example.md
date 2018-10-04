
## Run Alignment and quantification pipeline
Set all parameters in the params file to ensure they are correct for your file system or override them from command line as per Nextflow specs. If you prefer to run it locally, the -profile directive can be changed
To edit resources allocation to each step, edit the config file as appropriate.

    nextflow run RD_connect_RNA_align.nf\
        -c align.nextflow.config \
        -params-file params_align.yaml \
        -profile docker
    

## Run Variant Calling pipeline
Set all parameters in the params file to ensure they are correct for your file system or override them from command line as per Nextflow specs. If you prefer to run it locally, the -profile directive can be changed
To edit resources allocation to each step, edit the config file as appropriate.

    nextflow run RD_connect_RNA_call.nf \
      -params-file params_variant_call.yaml  \
      -profile docker \
      -c variant_call.nextflow.config

## Docker images
- broadinstitute/picard:2.18.2
- broadinstitute/gatk3:3.6-0
- quay.io/biocontainers/tabix:0.2.5--1
- quay.io/biocontainers/star:2.5.3a--0
- quay.io/biocontainers/rsem:1.3.0--boost1.64_2

