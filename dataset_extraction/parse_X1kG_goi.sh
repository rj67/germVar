X1kG_file_prefix='ALL.autosomes.phase3_shapeit2_mvncall_integrated_v4.20130502.sites'
#X1kG_file_prefix='ALL.wgs.phase1_release_v3.20101123.snps_indels.sites'

goi_bed='candidate_gene_hg19R_exons.bed'
Pfam_file='/home/rj67/seq/Pfam/Pfam_query_results.tsv'
HGNC_file='/home/rj67/seq/Pfam/HGNC_UniProt_ESTranscript.tsv'
Map_file='/home/rj67/seq/Mapability/wgEncodeCrgMapabilityAlign75mer_exome.bedGraph' 
 
# intersect with goi list, remove structural variants, annotate multiallele variants, split, left align
cat $X1kG_file_prefix.vcf | $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar intervals $goi_bed |  \
  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar filter "!( exists SVTYPE )"  | \
  vcfAnnoAlt.py |  sed '1,/^##/d' | vcfbreakmulti | \
  $BCFTOOLS_PATH/bcftools norm -f $REF_PATH/human_g1k_v37.fasta -O v - | \
  vcflength | java -jar $SNPEFF_PATH/SnpSift.jar varType - | \
  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
     -noStats -t -canon -no-downstream -no-upstream -no-intergenic -v GRCh37.74 -  | \
  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | vcfConvSnpEff.py | sed '1,/^##/d' |  \
 vcfAnnoPfam.py -p $Pfam_file -t $HGNC_file | sed '1,/^##/d' | vcfAnnoSTR.py -r $REF_PATH/human_g1k_v37.fasta  | sed '1,/^##/d' | \
 vcfannotate -b $Map_file -k Mappability > $X1kG_file_prefix.goi.vcf

