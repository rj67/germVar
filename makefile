######################################
#     query cgHub for all 25 cancers
######################################
for study in ACC BLCA LGG BRCA CESC COAD ESCA GBM HNSC KICH KIRC KIRP LIHC \
  LUAD OV PAAD PCPG PRAD READ SARC SKCM STAD THCA UCS UCEC
do 
   bsub -q short -W 40 -o $study.out -e $study.err -J $study ./query_parse.sh $study ./Sep_10_CGHub_query
   sleep 10
done


######################################
#     extract RNASeq data from TCGA data portal download
#####################################
cd /home/rj67/group/RNASeqV2
# split all disease id into study1/study2
cat study1/2 | xargs -i Rscript /home/rj67/opt/germVar-scripts/dataset_extraction/extractStudyRNASeqV2.R {}


######################################
#     lift over CCDS exon bed file
#####################################
# run first half of summary_CCDS.R to write out the bed file
#cd ~/opt/LiftOver
#liftOver CCDS_goi.exons.hg38.bed hg38ToHg19.over.chain CCDS_goi.exons.hg38Tohg19.bed CCDS_goi.exons.hg38Tohg19.unmapped.bed


######################################
#    extract relevant variants from COSMIC     
#####################################
for i in `seq 1 22`; do bsub -q short -W 40 -o $i.out -e $i.err -J $i ./dataset_extraction/parse_COSMIC_goi.sh $i; done 
# also do the above for X
vcf-concat CosmicCodingMuts.v71.goi.*.vcf | vcf-sort -c > CosmicCodingMuts.v71.goi.vcf
tabix_index.sh CosmicCodingMuts.v71.goi.vcf

######################################
#    annotate Clinvar_vcf written from Clinvar_txt
#####################################
./dataset_extraction/snpEff_clinvar.sh 

######################################
#   split the likely deleterious SNP list by chromosome
#####################################
for chr in `cut -f 1 nsSNP.MHS_Clinvar.bed | sort | uniq`; do grep -w $chr nsSNP.MHS_Clinvar.bed > nsSNP.MHS_Clinvar.$chr.bed ; done

######################################
#    Submit variant calling on nsSNP
#####################################
for i in `seq 1 22`; do bsub -q short -W 2:00 -R "rusage[mem=14000]" -o $i.out -e $i.err -J $i ./variant_calling/genotype_all_gvcfs_snp.sh $i; done 
# also do the above for X



