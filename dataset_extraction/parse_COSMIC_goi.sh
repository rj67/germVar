# filter out long variants, annotate variant type

COSMIC_file='/home/rj67/seq/COSMIC/CosmicCodingMuts_v70.vcf' 
COSMIC_out_file='/home/rj67/seq/COSMIC/CosmicCodingMuts_v70.anno.vcf' 

goi_bed=$REF_PATH/'candidate_gene_hg19_exons.bed'
Pfam_file='/home/rj67/seq/Pfam/Pfam_query_results.tsv'
Transcript_file='/home/rj67/seq/Pfam/HGNC_UniProt_ESTranscript.tsv'

vcf-sort -c $COSMIC_file  | $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar intervals $goi_bed |  \
  vcflength |  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar filter " ( length <= 100 ) & (length > -100)" | \
  $JAVA_PATH/java -jar $SNPEFF_PATH/SnpSift.jar varType -  | vcfTrimMNP.py | sed '1,/^##/d' |  \
  $JAVA_PATH/java -Xmx4g -jar $SNPEFF_PATH/snpEff.jar  -c $SNPEFF_PATH/snpEff.config   -noStats -t -canon -no-downstream -no-upstream -no-intergenic -v GRCh37.74 -  | \
  $SNPEFF_PATH/scripts/vcfEffOnePerLine.pl | vcfConvSnpEff.py | sed '1,/^##/d' | \
  vcfAnnoPfam.py -p $Pfam_file -t $Transcript_file | sed '1,/^##/d' > $COSMIC_out_file

Rscript ~/opt/germVar-scripts/variant_manipulation/convVCFDataFrame.R $COSMIC_out_file
