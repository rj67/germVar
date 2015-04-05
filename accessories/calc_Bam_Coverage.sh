java -Xmx2g -jar /home/rj67/opt/GATK/GenomeAnalysisTK.jar \
   -R /home/rj67/seq/Ref_v37/human_g1k_v37.fasta \
   -T DepthOfCoverage \
   -o orig_bam_coverage \
   -I orig_bam.list \
   -L /home/rj67/seq/Refseq_exon/GRCh37_refseq_exons_pm100_uniq.bed 

