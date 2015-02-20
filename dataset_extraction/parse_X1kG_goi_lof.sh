#X1kG_file_prefix='/home/rj67/seq/1000G/ALL.wgs.phase1_release_v3.20101123.snps_indels.sites'
X1kG_file_prefix='/home/rj67/seq/1000G/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites'

goi_bed=$REF_PATH/'candidate_gene_hg19_exons.bed'
#Pfam_file='/home/rj67/seq/Pfam/Pfam_query_results.tsv'
#HGNC_file='/home/rj67/seq/Pfam/HGNC_UniProt_ESTranscript.tsv'
Map_file='/home/rj67/seq/Mapability/wgEncodeCrgMapabilityAlign75mer_exome.bedGraph' 
 
# intersect with goi list, remove structural variants, annotate multiallele variants, split, left align
#gunzip -c $X1kG_file_prefix.vcf.gz | $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar intervals $goi_bed | \
#  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar filter "!( exists SVTYPE )"  | \
#  vcfAnnoAlt.py |  sed '1,/^##/d' | vcfbreakmulti | \
#  $BCFTOOLS_PATH/bcftools norm -f $REF_PATH/human_g1k_v37.fasta -O v - | \
#  vcflength | java -jar $SNPEFF_PATH/SnpSift.jar varType - |  vcfTrimMNP.py  | sed '1,/^##/d' |\
#  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
#  -noStats -t -no-downstream -no-upstream -no-intergenic -no-utr -no-intron -onlyTr $REF_PATH/list_goi_CCDS_transcript.txt -v GRCh37.75 - | \
#  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | grep -v 'EFF=sequence_feature'  | vcfConvSnpEff.py | sed '1,/^##/d' |  \
#  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar filter "( Impact = 'HIGH' )"  | \
#  vcfAnnoSTR.py -r $REF_PATH/human_g1k_v37.fasta  | sed '1,/^##/d' | \
#  vcfannotate -b $Map_file -k Mappability > $X1kG_file_prefix.goi.tmp.vcf

perl $VEP_PATH/variant_effect_predictor.pl --cache --offline --sift b -polyphen b --humdiv --ccds --uniprot  --hgvs --symbol --numbers --canonical --protein --biotype  \
    --no_stats --coding_only --plugin LoF,human_ancestor_fa:/home/rj67/.vep/Plugins/human_ancestor.fa.gz  -i $X1kG_file_prefix.goi.tmp.vcf -vcf -o STDOUT | \
    vcfConvVEP.py | sed '1,/^##/d' | vcfFixAllele.py | sed '1,/^##/d' > $X1kG_file_prefix.goi.lof.vcf 

#  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config  \
#     -noStats -t -canon -no-downstream -no-upstream -no-intergenic -v GRCh37.74 -  | \
#  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | vcfConvSnpEff.py | sed '1,/^##/d' |  \
# vcfAnnoPfam.py -p $Pfam_file -t $HGNC_file | sed '1,/^##/d' | vcfAnnoSTR.py -r $REF_PATH/human_g1k_v37.fasta  | sed '1,/^##/d' | \
# vcfannotate -b $Map_file -k Mappability > $X1kG_file_prefix.goi.vcf

