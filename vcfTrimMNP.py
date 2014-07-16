#!/usr/bin/env python

###################################################################

#vcfTrimMNP: trim trailing common nucleotides in MNP

###################################################################

import vcf
import sys
from vcf.model import _Record as VcfRecord
from vcf.parser import _Info as VcfInfo
from vcf.model import _Substitution as VcfSubstitution

if __name__ == '__main__':

  vcf_reader = vcf.Reader(sys.stdin)
  
  writer = vcf.Writer(sys.stdout, vcf_reader, lineterminator='\n')
  
  for Record in vcf_reader:
  
    # rely on SnpSift VARTYPE annotation
    # check for MNP flag is set and VARTYPE is SNP
    if Record.INFO.get('MNP', False) :
      if "SNP" in Record.INFO['VARTYPE'] :
  
        Record.REF = Record.REF[0]
        Record.ALT[0].sequence = str(Record.ALT[0])[0]
  
    writer.write_record(Record)
