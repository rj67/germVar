#!/bin/bash

studies=(`ls ./gvcf_lists/*.list | cut -f 3 -d "/" | cut -f 1 -d "_" `)
study=${studies[$1-1]}

ref_dir='/cbio/cslab/home/rj67/resources/Ref_v37'
out_dir='/cbio/cslab/home/rj67/working/haplo_redux/joint_vcf'
java_dir='/cbio/cslab/home/leew1/local/jre1.7.0_45/bin'
in_dir='/cbio/cslab/home/rj67/working/haplo_redux/merged_gvcfs'
gatk_dir='/cbio/cslab/home/rj67/local/GenomeAnalysisTK-3.0-0'
snpEff_dir='/cbio/cslab/home/rj67/local/snpEff'
bcftools_dir='/cbio/cslab/home/rj67/local/bin'
X1kG_dir='/cbio/cslab/home/rj67/resources/1000G'

bed_file='/cbio/cslab/home/rj67/working/haplo_redux/candidate_gene_hg19R_exons.bed'

$java_dir/java -Xmx8g -Djava.io.tmpdir=/scratch \
  -jar $gatk_dir/GenomeAnalysisTK.jar \
  -R $ref_dir/human_g1k_v37.fasta \
  -T GenotypeGVCFs \
  --dbsnp $ref_dir/dbsnp_138.b37.vcf \
  -V $in_dir/$study.merged.gvcf \
  -o $out_dir/$study.joint.vcf \
  -nt 4 \
  -L $bed_file  

cat $out_dir/$study.joint.vcf | vcfAnnoAlt.py | vcfbreakmulti | $bcftools_dir/bcftools norm -f $ref_dir/human_g1k_v37.fasta -O v - | \
 $java_dir/java -jar $snpEff_dir/SnpSift.jar varType - | vcfTrimMNP.py  | \
 $java_dir/java -Xmx4g -jar $snpEff_dir/snpEff.jar  -c $snpEff_dir/snpEff.config  \
     -noStats -t -canon -no-downstream -no-upstream -no-intergenic -v GRCh37.74 -  | \
  $snpEff_dir/scripts/vcfEffOnePerLine.pl | vcfConvSnpEff.py | vcfSelExon.py | \
 $java_dir/java -jar $snpEff_dir/SnpSift.jar annotate -noId -info 1kG_AF $X1kG_dir/ALL.wgs.phase1_release_v3.20101123.snps_indels.sites.vcf - | \
 vcfAnnoPfam.py  > $out_dir/$study.anno.vcf

vcffilter -f "Impact = HIGH" -f "QD > 1."  $out_dir/$study.anno.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2}}' >  $out_dir/$study.high.bed
