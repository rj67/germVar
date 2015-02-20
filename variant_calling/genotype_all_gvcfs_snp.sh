#!/bin/bash

#chrs=(`cut -f 1 $bed_prefix.bed | sort | uniq `)
#chr=${chrs[$1-1]}
#echo $chr
chr=$1

# I/O
project='norm_trunc'
out_dir='/home/rj67/group/results/Aug_5'
in_dir='/home/rj67/group/results/merged_gvcf'
bed_prefix='/home/rj67/group/results/trunc_bed/norm_joint_high'

#out_file=$project.$chr
out_file='nonsyn.tp53.2'

# reference file
Map_file='/home/rj67/seq/Mapability/wgEncodeCrgMapabilityAlign75mer_exome.bedGraph' 

$JAVA_PATH/java -Xmx12g \
  -jar $GATK_PATH/GenomeAnalysisTK.jar \
  -R $REF_PATH/human_g1k_v37.fasta \
  -T GenotypeGVCFs \
  --dbsnp $REF_PATH/dbsnp_138.b37.vcf \
  -V $in_dir/ACC.merged.gvcf  \
  -V $in_dir/BLCA.merged.gvcf\
  -V $in_dir/BRCA1.merged.gvcf\
  -V $in_dir/BRCA2.merged.gvcf\
  -V $in_dir/BRCA3.merged.gvcf\
  -V $in_dir/CESC.merged.gvcf\
  -V $in_dir/COAD.merged.gvcf\
  -V $in_dir/ESCA.merged.gvcf\
  -V $in_dir/GBM.merged.gvcf\
  -V $in_dir/HNSC.merged.gvcf\
  -V $in_dir/KICH.merged.gvcf\
  -V $in_dir/KIRC.merged.gvcf\
  -V $in_dir/KIRP.merged.gvcf\
  -V $in_dir/LGG.merged.gvcf\
  -V $in_dir/LIHC.merged.gvcf\
  -V $in_dir/LUAD.merged.gvcf\
  -V $in_dir/OV.merged.gvcf\
  -V $in_dir/PAAD.merged.gvcf\
  -V $in_dir/PCPG.merged.gvcf\
  -V $in_dir/PRAD.merged.gvcf\
  -V $in_dir/READ.merged.gvcf\
  -V $in_dir/SARC.merged.gvcf\
  -V $in_dir/SKCM.merged.gvcf\
  -V $in_dir/STAD.merged.gvcf\
  -V $in_dir/THCA.merged.gvcf\
  -V $in_dir/UCEC1.merged.gvcf\
  -V $in_dir/UCEC2.merged.gvcf\
  -V $in_dir/UCS.merged.gvcf\
  -o $out_dir/$out_file.vcf \
  -nt 1 \
  -L /home/rj67/tp53.bed
  #-L $bed_prefix.$chr.bed

# split multiallelic variants, put into GT vcf file.
cat $out_dir/$out_file.vcf | vcfAnnoAlt.py | sed '1,/^##/d' | \
  vcfbreakmulti | $BCFTOOLS_PATH/bcftools norm -f $REF_PATH/human_g1k_v37.fasta -O v - | \
  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar varType - | vcfTrimMNP.py  | sed '1,/^##/d' | \
  vcfkeepgeno - "GT" "AD" "DP" > $out_dir/$out_file.GT.vcf


## remove the genotype information, annotate variant
$BCFTOOLS_PATH/bcftools view -G $out_dir/$out_file.GT.vcf | vcflength | \
 $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
  -noStats -t -no-downstream -no-upstream -no-intergenic -no-utr -no-intron -onlyTr $REF_PATH/list_goi_GRCh37.75_long_transcript.txt -v GRCh37.75 - | \
 $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | grep -v 'EFF=sequence_feature'  | vcfConvSnpEff.py | sed '1,/^##/d' |  \
 $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar annotate -id $REF_PATH/dbsnp_138.b37.vcf - | \
 vcfAnnoSTR.py -r $REF_PATH/human_g1k_v37.fasta  | sed '1,/^##/d' | \
 vcfannotate -b $Map_file -k Mappability > $out_dir/$out_file.tmp.vcf

# for snp
#perl $VEP_PATH/variant_effect_predictor.pl --cache --offline --everything --humdiv --no_stats --coding_only --plugin LoF,human_ancestor_fa:/home/rj67/.vep/Plugins/human_ancestor.fa.gz  \
#  -i  $out_dir/$project.$chr.tmp.vcf -vcf -o STDOUT | vcfConvVEP.py | sed '1,/^##/d' | vcfFixAllele.py | sed '1,/^##/d' > $out_dir/$project.$chr.VAR.vcf

#perl $VEP_PATH/variant_effect_predictor.pl --cache --offline --everything --humdiv --no_stats --plugin LoF,human_ancestor_fa:/home/rj67/.vep/Plugins/human_ancestor.fa.gz  \
#  -i  $out_dir/$project.$chr.tmp.vcf -vcf -o STDOUT | vcfConvVEP.py | sed '1,/^##/d' | vcfFixAllele.py | sed '1,/^##/d' > $out_dir/$project.$chr.VAR.vcf

#sed -i '1s/^/##fileformat=VCFv4.1\n/' $out_dir/$project.$chr.VAR.vcf
#sed -i '1s/^/##fileformat=VCFv4.1\n/' $out_dir/$project.$chr.GT.vcf
