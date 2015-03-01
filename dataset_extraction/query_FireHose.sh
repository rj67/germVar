for study in ACC BLCA LGG BRCA CESC COAD ESCA GBM HNSC KICH KIRC KIRP LIHC \
  LUAD LUSC OV PAAD PCPG PRAD READ SARC SKCM STAD THCA UCS UCEC
do
  # Focal deletion
  wget http://gdac.broadinstitute.org/runs/analyses__2014_10_17/reports/cancer/$study-TP/CopyNumber_Gistic2/del_genes.conf_99.txt
  mv del_genes.conf_99.txt $study.del_genes.conf_99
done
