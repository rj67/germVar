# filter out long variants, annotate variant type

ESP_file_prefix='/home/rj67/seq/ESP6500/ESP6500SI-V2-SSA137.updatedProteinHgvs.snps_indels'

goi_bed=$REF_PATH/'candidate_gene_hg19_exons.bed'
Map_file='/home/rj67/seq/Mapability/wgEncodeCrgMapabilityAlign75mer_exome.bedGraph'

# intersect with goi list, remove structural variants, annotate multiallele variants, split, left align
cat $ESP_file_prefix.vcf | $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar intervals $goi_bed | \
  splitESPMultiAllele.py |  sed '1,/^##/d'| vcfAnnoAlt.py |  sed '1,/^##/d' | vcfbreakmulti | \
  $BCFTOOLS_PATH/bcftools norm -f $REF_PATH/human_g1k_v37.fasta -O v - | \
  vcflength | java -jar $SNPEFF_PATH/SnpSift.jar varType - |  vcfTrimMNP.py  | sed '1,/^##/d' |\
  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
  -noStats -t -no-downstream -no-upstream -no-intergenic -no-utr -no-intron -onlyTr $REF_PATH/list_goi_CCDS_transcript.txt -v GRCh37.75 - | \
  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | grep -v 'EFF=sequence_feature'  | vcfConvSnpEff.py | sed '1,/^##/d' |  \
  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar filter "( Impact = 'HIGH' )"  | \
  vcfAnnoSTR.py -r $REF_PATH/human_g1k_v37.fasta  | sed '1,/^##/d' | \
  vcfannotate -b $Map_file -k Mappability > $ESP_file_prefix.goi.tmp.vcf

perl $VEP_PATH/variant_effect_predictor.pl --cache --offline --sift b -polyphen b --humdiv --ccds --uniprot  --hgvs --symbol --numbers --canonical --protein --biotype  \
    --no_stats --plugin LoF,human_ancestor_fa:/home/rj67/.vep/Plugins/human_ancestor.fa  -i $ESP_file_prefix.goi.tmp.vcf -vcf -o STDOUT | \
    vcfConvVEP.py | sed '1,/^##/d' | vcfFixAllele.py | sed '1,/^##/d' |\
    vcfAnnoANC.py --anc_ref ~/.vep/Plugins/human_ancestor.fa --ref $REF_PATH/human_g1k_v37.fasta | sed '1,/^##/d' >  $ESP_file_prefix.goi.lof.vcf
~                                        
