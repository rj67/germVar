#!/bin/bash

bam_list=$ARGUMENTS['bam_list']
bed_file=$ARGUMENTS['bed_file']
out_dir=$ARGUMENTS['out_dir']
project=$ARGUMENTS['project']

#out_dir='/cbio/cslab/home/rj67/working/haplo_redux'

$JAVA_PATH/java -Xmx8g 
  -jar $GATK_PATH/GenomeAnalysisTK.jar \
  -R $REF_PATH/human_g1k_v37.fasta \
  -T HaplotypeCaller \
  --dbsnp $REF_PATH/dbsnp_138.b37.vcf \
  -stand_call_conf 30.0 \
  -stand_emit_conf 30.0 \
  -I $bam_list
  -minPruning 2  \
  --maxNumHaplotypesInPopulation 256 \
  --max_alternate_alleles 6 \
  --pcr_indel_model CONSERVATIVE \
  -A VariantType \
  --mergeVariantsViaLD \
  -bamout $out_dir/$project.joint.bam \
  -o $out_dir/$project.haplo.vcf \
  -L $bed_file
  -nct 4
  


#-o $out_dir/$project.$bednum.raw.vcf \

