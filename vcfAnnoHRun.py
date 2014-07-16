#!/usr/bin/env python

import vcf
import sys
from vcf.parser import _Info as VcfInfo

vcf_reader = vcf.Reader(sys.stdin)

vcf_reader.infos['HRUN'] = VcfInfo( 'HRUN', 0, 'Flag', 'whether there is HRUN')
vcf_reader.infos['HRUN_length'] = VcfInfo( 'HRUN_legnth', 1, 'Integer', 'length of HRUN')

writer = vcf.Writer(sys.stdout, vcf_reader, lineterminator='\n')

for Record in vcf_reader:

  if abs(Record.INFO['length'][0]) == 1 : 
    if Record.INFO['length'][0] == 1 :
      RU = Record.ALT[0].sequence[1]
    elif Record.INFO['length'][0] == -1 :
      RU = Record.REF[1]
    
    primer_strip = Record.INFO['Primer'].lstrip(RU)
    if len(Record.INFO['Primer']) - len(primer_strip) >= 5 :
      Record.add_info('HRUN', True)
      Record.add_info('HRUN_legnth',len(Record.INFO['Primer']) - len(primer_strip))

  writer.write_record(Record)
#    for i in range(len(Record.ALT)):
#      INFO_multi =[]
#      for key in multi_keys:
#        if Record.INFO.get(key, None):
#          INFO_multi.append(key + "=" + str(Record.INFO[key][i]))
#      #INFO_multi = [key + "=" + str(Record.INFO[key][i]) for key in multi_keys]
#
#      # check whether RPA annotation is present, if yes, add it, the RPA includes reference so offset by 1
#      if Record.INFO.get('RPA', None ):
#        INFO_multi.append('RPA='+ str(Record.INFO["RPA"][i+1]))
#      
#      INFO_out = ';'.join(INFO_other + INFO_multi + ['ALT_id='+str(i+1)+"-"+str(Record.CHROM)+':'+str(Record.POS)])
#
#      if i > 0 :
#        FILTER = "MultiAlt"
#      print "\t".join([str(Record.CHROM), str(Record.POS), ID, str(Record.REF), str(Record.ALT[i]), str(Record.QUAL), FILTER, INFO_out, Record.FORMAT])
#        
#  # otherwise, use Pyvcf to print line
#  else:
#     writer.write_record(Record)
