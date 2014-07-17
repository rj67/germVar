#!/usr/bin/python
import xml.etree.ElementTree as ET
from pandas import *
import urllib

# read in all the HGNC proteins, query the Pfam, record the family info

uniprots = open('all_HGNC_UniProt.txt','r').read().split()
#print uniprots[0:100]

lists = []
for i in range(len(uniprots)):
  uniprot = uniprots[i]
  print i
  url = 'http://pfam.sanger.ac.uk/protein/'+uniprot+'?output=xml'
  
  tree = ET.parse(urllib.urlopen(url))
     
  root = tree.getroot()
  
  #print root.tag, root.attrib
  
  for child in root:
    #print child.tag, child.attrib
    matches = child.find('{http://pfam.sanger.ac.uk/}matches')
    if not matches is None:
      
      for match in matches.findall('{http://pfam.sanger.ac.uk/}match'):
        
        for loc in match:
          
          entry = dict(match.attrib.items() + loc.attrib.items())
          entry['uniprot'] = uniprot
          lists.append(entry)
        #$print match.find('location').tag, match.find('location').attrib
    else:
      print uniprot, 'no match'

df = DataFrame(lists, columns=['uniprot', 'type', 'id', 'start', 'end', 'evalue', 'accession','ali_start', 'ali_end', 'bitscore', 'hmm_start', 'hmm_end']) 

#print df.to_string()
df.to_csv('pfam_hgnc.txt', sep = '\t', index = False)

#for record in root.findall('entry'):
#  print record.find('description').text
