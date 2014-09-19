ClinVar_file_prefix='/home/rj67/seq/ClinVar/clinvar_20140902'

goi_bed=$REF_PATH'/candidate_gene_hg19_exons.bed'
Pfam_file='/home/rj67/seq/Pfam/Pfam_query_results.tsv'
Transcript_file='/home/rj67/seq/Pfam/HGNC_UniProt_ESTranscript.tsv'
 
# intersect with goi list, remove structural variants, annotate multiallele variants, split, left align
cat $ClinVar_file_prefix.vcf | java -jar $SNPEFF_PATH/SnpSift.jar intervals $goi_bed |  \
   vcfbreakmulti | $BCFTOOLS_PATH/bcftools norm -f $REF_PATH/human_g1k_v37.fasta -O v - | \
  vcflength | java -jar $SNPEFF_PATH/SnpSift.jar varType - | \
  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config   -noStats -t -canon -no-downstream -no-upstream -no-intergenic -v GRCh37.74 -  |
  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | vcfConvSnpEff.py | sed '1,/^##/d' | \
  vcfAnnoPfam.py -p $Pfam_file -t $Transcript_file | sed '1,/^##/d'  > $ClinVar_file_prefix.goi.vcf
