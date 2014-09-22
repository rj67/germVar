filenames=$(</dev/stdin)

for filename in $filenames
do

  RGHead=$(samtools view -H $filename | grep -m 1 '^@RG')
   
  RGSM=$(echo $RGHead | tr " " "\n" | grep '^SM:' | cut -f 2 -d ':')

  echo $RGSM


done

