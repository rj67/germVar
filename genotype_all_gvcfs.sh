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
# program
java_dir='/home/rj67/opt/jdk1.7.0_25/bin'
gatk_dir='/home/rj67/opt/GenomeAnalysisTK-3.2-2'
snpEff_dir='/home/rj67/opt/snpEff'
bcftools_dir='/home/rj67/opt/bcftools-master'
annovar_dir='/home/rj67/opt/annovar'
# reference
ref_dir='/home/rj67/seq/Ref_v37'
X1kG_file='/home/rj67/seq/1000G/coding/ALL.wgs.integrated_phase1_release_v3_coding_annotation.20101123.snps_indels.sites.rename.vcf'
ESP_file='/home/rj67/seq/ESP6500/ESP6500SI-V2-SSA137.updatedProteinHgvs.snps_indels.split.vcf'
Pfam_file='/home/rj67/seq/Pfam/Pfam_query_results.tsv'
HGNC_file='/home/rj67/seq/Pfam/HGNC_UniProt_ESTranscript.tsv'
Map_file='/home/rj67/seq/Mapability/wgEncodeCrgMapabilityAlign75mer_exome.bedGraph' 

#$java_dir/java -Xmx12g -Djava.io.tmpdir=/scratch \
$java_dir/java -Xmx12g \
  -jar $gatk_dir/GenomeAnalysisTK.jar \
  -R $ref_dir/human_g1k_v37.fasta \
  -T GenotypeGVCFs \
  --dbsnp $ref_dir/dbsnp_138.b37.vcf \
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
  -o $out_dir/$project.$chr.vcf \
  -nt 1 \
  -L $bed_prefix.$chr.bed

cat $out_dir/$project.$chr.vcf | vcfAnnoAlt.py | sed '1,/^##/d' | vcfbreakmulti | $bcftools_dir/bcftools norm -f $ref_dir/human_g1k_v37.fasta -O v - | \
 $java_dir/java -jar $snpEff_dir/SnpSift.jar rmRefGen | vcfkeepgeno - "GT" "AD" "DP" "GQ" | \
 $java_dir/java -jar $snpEff_dir/SnpSift.jar varType - | vcfTrimMNP.py  | sed '1,/^##/d' | \
 $java_dir/java -Xmx4g -jar $snpEff_dir/snpEff.jar  -c $snpEff_dir/snpEff.config  \
     -noStats -t -canon -no-downstream -no-upstream -no-intergenic -v GRCh37.74 -  | \
  $snpEff_dir/scripts/vcfEffOnePerLine.pl | vcfConvSnpEff.py | sed '1,/^##/d' |  \
 $java_dir/java -jar $snpEff_dir/SnpSift.jar filter "( ( Impact = 'HIGH' ) | ( Impact = 'MODERATE' )) " | \
 $java_dir/java -jar $snpEff_dir/SnpSift.jar annotate -noId -info X1kG_AF,X1kG_LDAF,X1kG_ERATE $X1kG_file - | \
 $java_dir/java -jar $snpEff_dir/SnpSift.jar annotate -noId -info ESP_AC,ESP_AN,ESP_AF $ESP_file - | \
 $java_dir/java -jar $snpEff_dir/SnpSift.jar annotate -id $ref_dir/dbsnp_138.b37.vcf - | \
 vcfAnnoPfam.py -p $Pfam_file -t $HGNC_file | sed '1,/^##/d' | vcfAnnoSTR.py -r $ref_dir/human_g1k_v37.fasta  | sed '1,/^##/d' | \
 vcfannotate -b $Map_file -k Mappability > $out_dir/$project.$chr.anno.vcf

# use annovar to add alternative annotation
table_annovar.pl  $out_dir/$project.$chr.anno.vcf $annovar_dir/humandb/ -out $out_dir/$project.$chr -remove -protocol ensGene -operation g  -buildver hg19 -vcfinput
sed -i '1s/^/##INFO=<ID=AAChange.ensGene,Number=1,Type=String,Description="annovar AAChange annotation">\n/' $out_dir/$project.$chr.hg19_multianno.vcf
sed -i '1s/^/##INFO=<ID=GeneDetail.ensGene,Number=1,Type=String,Description="annovar Splicing annotation">\n/' $out_dir/$project.$chr.hg19_multianno.vcf
#annotate the vcf file
$java_dir/java -jar $snpEff_dir/SnpSift.jar annotate -noId -info AAChange.ensGene,GeneDetail.ensGene $out_dir/$project.$chr.hg19_multianno.vcf $out_dir/$project.$chr.anno.vcf > $out_dir/$project.$chr.annovar.vcf 
rm $out_dir/$project.$chr.hg19_multianno.vcf
rm $out_dir/$project.$chr.avinput

sed -i '1s/^/##fileformat=VCFv4.1\n/' $out_dir/$project.$chr.annovar.vcf
