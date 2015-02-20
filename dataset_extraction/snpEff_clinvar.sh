
clinvar_prefix='/home/rj67/seq/ClinVar/clinvar_txt.20150109'

$JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config   \
  -noStats -t -no-downstream -no-upstream -no-intergenic -no-utr -no-intron -v GRCh37.75 -onlyTr $REF_PATH/list_goi_GRCh37.75_long_transcript.txt $clinvar_prefix.vcf  | \
  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | grep -v 'EFF=sequence_feature'| vcfConvSnpEff.py | sed '1,/^##/d' | \
  vcfFixAllele.py | sed '1,/^##/d'  > $clinvar_prefix.snpEff.vcf

tabx_index.sh $clinvar_prefix.snpEff.vcf

