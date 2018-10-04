## Run Alignment and quantification pipeline

nextflow run RD_connect_RNA_align.nf\
    -c align.nextflow.config \
    -params-file params_align.yaml \
    -profile docke
   

## Run Variant Calling pipeline

nextflow run RD_connect_RNA_call.nf \
  -params-file params_variant_call.yaml  \
  -profile docker \
  -c variant_call.nextflow.config

