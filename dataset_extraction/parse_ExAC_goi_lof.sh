ExAC_file_prefix='/cbio/cslab/home/rj67/resources/ExAC/ExAC.r0.3.sites.vep'

goi_bed=$REF_PATH/'candidate_gene_hg19_exons.bed'
Map_file='/cbio/cslab/home/rj67/resources/Mapability/wgEncodeCrgMapabilityAlign75mer_exome.bedGraph' 

chr=$1

tabix -h $ExAC_file_prefix.vcf.gz $chr |\
  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar filter "(( FILTER= 'PASS')) " |\
  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar rmInfo - GQ_HIST DP_HIST |\
  grep -v '^##GATK\|^##GVCFBlock\|^##reference=\|^##source\|^##SnpSift\|^##VEP\|^##LoF' |\
  $JAVA_PATH/java -Xmx2g -jar $SNPEFF_PATH/SnpSift.jar intervals $goi_bed | \
  vcfAnnoAlt.py |  vcfbreakmulti | \
  $BCFTOOLS_PATH/bcftools norm -f $REF_PATH/human_g1k_v37.fasta -O v - | \
  vcflength | $JAVA_PATH/java -Xmx2g  -jar $SNPEFF_PATH/SnpSift.jar varType - |  vcfTrimMNP.py  | \
  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
  -noStats -t -no-downstream -no-upstream -no-intergenic -no-utr -no-intron -no REGULATION -onlyTr $REF_PATH/list_goi_GRCh37.75_CCDS_transcript.txt -v GRCh37.75 - | \
  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | grep -v 'EFF=sequence_feature'  | vcfConvSnpEff.py |  \
  $JAVA_PATH/java -Xmx2g -jar $SNPEFF_PATH/SnpSift.jar filter "( Impact = 'HIGH' )"  | \
  vcfAnnoSTR.py -r $REF_PATH/human_g1k_v37.fasta  | sed '1,/^##/d' | \
  vcfannotate -b $Map_file -k Mappability > $ExAC_file_prefix.$chr.tmp.vcf

perl $VEP_PATH/variant_effect_predictor.pl --cache --offline --ccds --uniprot  --hgvs --symbol --numbers --canonical --protein --biotype  \
    --no_stats --plugin LoF,human_ancestor_fa:/cbio/cslab/home/rj67/.vep/Plugins/human_ancestor.fa  -i $ExAC_file_prefix.$chr.tmp.vcf -vcf -o STDOUT | \
    vcfConvVEP.py | vcfFixAllele.py > $ExAC_file_prefix.$chr.lof.vcf 


