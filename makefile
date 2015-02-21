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
#   Single sample Haplotypecaller call in gvcf mode
#####################################
#!/bin/sh
#PBS -N single_gvcf
#PBS -l mem=13gb,walltime=4:00:00,nodes=1:ppn=1
#PBS -q batch
#PBS -j oe
##PBS -k oe
cd $PBS_O_WORKDIR
bam_file=`head -$PBS_ARRAYID /cbio/cslab/home/rj67/bam/Feb_2_later_bam_uid.list | tail -1 | cut -f 2`
anal_id=`head -$PBS_ARRAYID /cbio/cslab/home/rj67/bam/Feb_2_laterbam_uid.list | tail -1 | cut -f 1`
out_dir='/cbio/cslab/home/rj67/vcf/single_gvcfs'
bed_file=$REF_PATH/'candidate_gene_hg19_exons.bed'

haplotypeCaller_gvcf.sh -l=$bam_file -b=$bed_file -o=$out_dir -n=$anal_id

######################################
#   Merge single gvcf into small batches
# #####################################
#PBS -N Combine
#PBS -l mem=17gb,walltime=8:00:00,nodes=1:ppn=1
#PBS -q batch
#PBS -j oe
##PBS -k oe
cd $PBS_O_WORKDIR
studies=(`ls ./gvcf_lists/*.list | cut -f 3 -d "/" | cut -f 1,2 -d "_" `)
study=${studies[$PBS_ARRAYID]}

/cbio/cslab/home/rj67/local/germVar/variant_calling/combine_gvcfs.sh $study

######################################
#   Genotype the batch merged gvcfs, extract the lof position
#   # #####################################
#PBS -N Genotype_batch
#PBS -l mem=16gb,walltime=8:00:00,nodes=1:ppn=4
#PBS -q batch
#PBS -j oe
##PBS -k oe
cd $PBS_O_WORKDIR
studies=(`ls ./gvcf_lists/*.list | cut -f 3 -d "/" | cut -f 1,2 -d "_" `)
study=${studies[$PBS_ARRAYID]}

/cbio/cslab/home/rj67/local/germVar/variant_calling/genotype_batch_gvcfs.sh $study



######################################
#    Submit variant calling on nsSNP
#####################################
for i in `seq 1 22`; do bsub -q short -W 2:00 -R "rusage[mem=14000]" -o $i.out -e $i.err -J $i ./variant_calling/genotype_all_gvcfs_snp.sh $i; done 
# also do the above for X



