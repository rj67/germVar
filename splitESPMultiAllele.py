#!/usr/bin/env python

# parse the ESP vcf file, split the multiallelic, pipe to bcftools to left align indels 
# ./splitESPMultiAllele.py | /home/rj67/opt/bcftools-master/bcftools norm -f /home/rj67/seq/Ref_v37/human_g1k_v37.fasta -O v - > ESP6500SI-V2-SSA137.updatedProteinHgvs.snps_indels.split.vcf

import vcf
from math import log10
from math import floor
import sys
from vcf.model import _Record as VcfRecord
from vcf.parser import _Info as VcfInfo

ESP_file='./ESP6500SI-V2-SSA137.updatedProteinHgvs.snps_indels.vcf'

vcf_reader = vcf.Reader(open(ESP_file, 'r'))

vcf_reader.infos['ALT_idx'] = VcfInfo( 'ALT_idx', '1', 'String', 'index for the alternative alleles')
vcf_reader.infos['ALT_pos'] = VcfInfo( 'ALT_pos', 1, 'String', 'original postition for the multiallele complex')
vcf_reader.infos['ALT_num'] = VcfInfo( 'ALT_num', 1, 'Integer', 'number of alternative allele for the multiallele complex')

vcf_reader.infos['ESP_AC'] = VcfInfo( 'ESP_AC', 1, 'Integer', 'minor allele count')
vcf_reader.infos['ESP_AN'] = VcfInfo( 'ESP_AN', 1, 'Integer', 'total allele count')
vcf_reader.infos['ESP_AF'] = VcfInfo( 'ESP_AF', 1, 'Float', 'minor allele frequency')

vcf_reader.infos['ESP_EA_AF'] = VcfInfo( 'ESP_EA_AF', 1, 'Float', 'European minor allele frequency')
vcf_reader.infos['ESP_AA_AF'] = VcfInfo( 'ESP_AA_AF', 1, 'Float', 'African minor allele frequency')

writer = vcf.Writer(sys.stdout, vcf_reader, lineterminator='\n')

# info field that might be multi-allelic
for Record in vcf_reader:
  
  # when encountering multiallele, modify record and print line
  TACs = map(int, Record.INFO['TAC'])
  EA_ACs = map(int, Record.INFO['EA_AC'])
  AA_ACs = map(int, Record.INFO['AA_AC'])
  
  if len(Record.ALT) > 1 :
    ALTs = [ alt.sequence for alt in Record.ALT ]
    NALT = len(Record.ALT)

    for i in range(NALT):
      
      Record.add_info('ESP_AC', TACs[i])
      Record.add_info('ESP_AN', sum(TACs))
      ESP_AF = float(TACs[i]) / float(sum(TACs))
      ESP_EA_AF = float(EA_ACs[i]) / float(sum(EA_ACs))
      ESP_AA_AF = float(AA_ACs[i]) / float(sum(AA_ACs))
      if not ESP_AF == 0.0:
         ESP_AF = round(ESP_AF, 1-int(floor(log10(ESP_AF))))
      Record.add_info('ESP_AF', ESP_AF)
      if not ESP_AA_AF == 0.0:
         ESP_AA_AF = round(ESP_AA_AF, 1-int(floor(log10(ESP_AA_AF))))
      Record.add_info('ESP_AA_AF', ESP_AA_AF)
      if not ESP_EA_AF == 0.0:
         ESP_EA_AF = round(ESP_EA_AF, 1-int(floor(log10(ESP_EA_AF))))
      Record.add_info('ESP_EA_AF', ESP_EA_AF)
      #  Record.add_info('ESP_EA_AF', round(ESP_EA_AF, 1-int(floor(log10(ESP_EA_AF)))))
      #  Record.add_info('ESP_AA_AF', round(ESP_AA_AF, 1-int(floor(log10(ESP_AA_AF)))))
      
      Record.add_info('ALT_idx', str(i+1))
      Record.add_info('ALT_pos', str(Record.CHROM) + ':' + str(Record.POS))
      Record.add_info('ALT_num', NALT)

      del Record.ALT[1:] 
      Record.ALT[0].sequence = ALTs[i]

      writer.write_record(Record)

  else:

    Record.add_info('ESP_AC', TACs[0])
    Record.add_info('ESP_AN', sum(TACs))
    ESP_AF = float(TACs[0]) / float(sum(TACs))
    ESP_EA_AF = float(EA_ACs[0]) / float(sum(EA_ACs))
    ESP_AA_AF = float(AA_ACs[0]) / float(sum(AA_ACs))
    if not ESP_AF == 0.0:
       ESP_AF = round(ESP_AF, 1-int(floor(log10(ESP_AF))))
    Record.add_info('ESP_AF', ESP_AF)
    if not ESP_AA_AF == 0.0:
       ESP_AA_AF = round(ESP_AA_AF, 1-int(floor(log10(ESP_AA_AF))))
    Record.add_info('ESP_AA_AF', ESP_AA_AF)
    if not ESP_EA_AF == 0.0:
       ESP_EA_AF = round(ESP_EA_AF, 1-int(floor(log10(ESP_EA_AF))))
    Record.add_info('ESP_EA_AF', ESP_EA_AF)
    
    writer.write_record(Record)

