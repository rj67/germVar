if [ $# -lt 1 ] ; then
  echo "usage: $0 list_of_vcf_to_evaluate"
  exit 1
fi

for filename in $@
do 
(
  prefix=${filename/.vcf/''}
  
  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
        -noStats -t -no-downstream -no-upstream -no-intergenic -no-utr -no-intron -onlyTr $REF_PATH/list_goi_CCDS_transcript.txt -v GRCh37.75 $prefix.vcf | \
        $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | grep -v 'EFF=sequence_feature' | vcfConvSnpEff.py | sed '1,/^##/d'  > $prefix.snpEff.vcf
  
) 
done

