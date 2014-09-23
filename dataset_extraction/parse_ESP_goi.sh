# filter out long variants, annotate variant type

ESP_file_prefix='ESP6500SI-V2-SSA137.updatedProteinHgvs.snps_indels'

goi_bed=$REF_PATH/'candidate_gene_hg19_exons.bed'
Pfam_file='/home/rj67/seq/Pfam/Pfam_query_results.tsv'
HGNC_file='/home/rj67/seq/Pfam/HGNC_UniProt_ESTranscript.tsv'
Map_file='/home/rj67/seq/Mapability/wgEncodeCrgMapabilityAlign75mer_exome.bedGraph' 
 
# intersect with goi list, remove structural variants, annotate multiallele variants, split, left align
cat $ESP_file_prefix.vcf | $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar intervals $goi_bed | \
  splitESPMultiAllele.py |  sed '1,/^##/d'| vcfAnnoAlt.py |  sed '1,/^##/d' | vcfbreakmulti | \
  $BCFTOOLS_PATH/bcftools norm -f $REF_PATH/human_g1k_v37.fasta -O v - | \
  vcflength | java -jar $SNPEFF_PATH/SnpSift.jar varType - | \
  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
     -noStats -t -canon -no-downstream -no-upstream -no-intergenic -v GRCh37.74 -  | \
  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | vcfConvSnpEff.py | sed '1,/^##/d' |  \
 vcfAnnoPfam.py -p $Pfam_file -t $HGNC_file | sed '1,/^##/d' | vcfAnnoSTR.py -r $REF_PATH/human_g1k_v37.fasta  | sed '1,/^##/d' | \
 vcfannotate -b $Map_file -k Mappability > $ESP_file_prefix.goi.vcf

