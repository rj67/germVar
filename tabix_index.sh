######################################################
# tabix index a list of vcf
######################################################

if [ $# -lt 1 ] ; then
  echo "usage: $0 list_of_vcf_to_index"
  exit 1
fi

for vcf in $@
do 
  cat $vcf | grep -v '^##FILTER' | bgzip | tee > $vcf.gz
  tabix -f -p vcf $vcf.gz
done
