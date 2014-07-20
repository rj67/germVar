#!/usr/bin/env python

import vcf
import sys
from vcf.parser import _Info as VcfInfo

vcf_reader = vcf.Reader(sys.stdin)

vcf_reader.infos['ALT_idx'] = VcfInfo( 'ALT_idx', 'A', 'String', 'index for the alternative alleles')
vcf_reader.infos['ALT_pos'] = VcfInfo( 'ALT_pos', 1, 'String', 'original postition for the multiallele complex')
vcf_reader.infos['ALT_num'] = VcfInfo( 'ALT_num', 1, 'Integer', 'number of alternative allele for the multiallele complex')

writer = vcf.Writer(sys.stdout, vcf_reader, lineterminator='\n')

# info field that might be multi-allelic
for Record in vcf_reader:

  # when encountering multiallele, annotate alt sequence
  if len(Record.ALT) > 1 :

    ALT_idx = ','.join([ str(x) for x in range(1, len(Record.ALT)+ 1)])
    ALT_pos = str(Record.CHROM) + ':' + str(Record.POS)
    ALT_num = str(len(Record.ALT))

    Record.add_info('ALT_idx',ALT_idx)
    Record.add_info('ALT_pos',ALT_pos)
    Record.add_info('ALT_num',ALT_num)
    writer.write_record(Record)

  # otherwise, use Pyvcf to print line
  else:
    writer.write_record(Record)

