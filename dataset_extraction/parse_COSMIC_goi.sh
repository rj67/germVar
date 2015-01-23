# filter out long variants, annotate variant type

COSMIC_prefix='/home/rj67/seq/COSMIC/CosmicCodingMuts.v71' 
COSMIC_out='/home/rj67/seq/COSMIC/CosmicCodingMuts.v71.goi' 

goi_bed=$REF_PATH/'candidate_gene_hg19_exons.bed'

chr=$1

tabix -h $COSMIC_prefix.vcf.gz $chr  | $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar intervals $goi_bed |  \
  vcflength |  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar filter " ( length <= 100 ) & (length > -100)" | \
  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar varType -  | vcfTrimMNP.py | sed '1,/^##/d' |  \
  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
   -noStats -t -canon -no-downstream -no-upstream -no-intergenic -no-utr -no-intron -v GRCh37.75 -  | \
  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | grep -v 'EFF=sequence_feature'| vcfConvSnpEff.py | sed '1,/^##/d' > $COSMIC_out.$chr.vcf

