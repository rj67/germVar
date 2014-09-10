grep -v 'hap\|Un_gl' candidate_gene_hg19R_exons.bed | sort -k1,1V -k2,2n | mergeBed -i - > tmp
mv tmp candidate_gene_hg19R_exons.bed
