######################################################
# tabix index a list of vcf
######################################################

if [ $# -lt 1 ] ; then
  echo "usage: $0 list_of_vcf_to_index"
  exit 1
fi

chrs=(`seq 1 22`)
chrs+=("X")

for vcf in $@
do 
  for chr in "${chrs[@]}"
  do  
    vcfpre=${vcf/.vcf.gz/''}
    echo $chr
    tabix -h $vcfpre.vcf.gz $chr | grep -v '^##FILTER' | bgzip | tee > $vcfpre.$chr.vcf.gz
    tabix -f -p vcf $vcfpre.$chr.vcf.gz
  done
done
