#!/usr/bin/python

import xml.etree.ElementTree as ET
import argparse
import sys
import os
import numpy as np
from pandas import *
from os import listdir

parser = argparse.ArgumentParser()
parser.add_argument('study',help="The TCGA study abbreviation")
parser.add_argument('output_dir',help="directory for the xml and output file")
 
opts = parser.parse_args()

xml_files = [ opts.output_dir + '/' + opts.study + '.xml' ]
print xml_files
output_csv_file = opts.output_dir + '/' + opts.study+'.csv'

pandas.set_option('max_colwidth',200)

lists = []
for file in xml_files:
  tree = ET.parse(file)
  root = tree.getroot()
  
  for record in root.getiterator('Result'):
    uuid = record.find('analysis_id').text
    participant = record.find('participant_id').text
    disease = record.find('disease_abbr').text
    legacy_sample_id = record.find('legacy_sample_id').text
    sample_id = legacy_sample_id.rsplit('-',2)[0]
    upload_date = record.find('upload_date').text
    sample_type = record.find('sample_type').text
    platform = record.find('platform')
    if platform != None :
      platform = platform.text
      if platform == None : platform = '-'
    else:
      platform = '-'
    analyte = record.find('analyte_code').text
    refassem = record.find('refassem_short_name').text
    refassem_abbr = "19" if ("19" in refassem or "37" in refassem) else "18"
    center = record.find('center_name').text
    ###
    file = record.find('files').find('file')
    filename = file.find('filename').text
    filesize = file.find('filesize').text
    ###
    analysis = record.find('analysis_xml').find('ANALYSIS_SET').find('ANALYSIS')
    anl_title = analysis.find('TITLE').text
    processing = '|'.join(set([elem.text for elem in analysis.find('ANALYSIS_TYPE').find('REFERENCE_ALIGNMENT').find('PROCESSING').getiterator('PROGRAM')]))
    #processing = '|'.join(set([elem.text for elem in record.getinterator('PROGRAM')]))
    
    ###
    experiment = record.find('experiment_xml').find('EXPERIMENT_SET').find('EXPERIMENT')
    exp_title = experiment.find('TITLE')
    if exp_title != None:
      exp_title = exp_title.text
      if exp_title == None : exp_title = '-'
    else:
      exp_title = '-'
    exp_library = experiment.find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_NAME').text
    ####
    exp_length = '-'
    for elem in experiment.getiterator('SPOT_LENGTH'):
        exp_length = elem.text
  
    #convert file size
    filesize = round(float(filesize)/pow(1024,3),1)
 
    #if platform not found check in filename and exp_title
    if platform == '-':
      if 'Illumina' in filename: platform = 'ILLUMINA'
    if platform == '-':
      if 'Illumina' in exp_title: platform = 'ILLUMINA'
    if platform == '-':
      if 'SOLiD' in filename: platform = 'ABI_SOLID'
   
    lists.append((uuid,participant,sample_id,upload_date,disease,sample_type,platform,analyte,refassem_abbr,center,filename,filesize,exp_library,exp_length,processing,legacy_sample_id,anl_title,exp_title,refassem))


df=DataFrame(lists,columns=['analysis','participant','sample_id','upload_date','disease','sample_type','platform','analyte','refassem_abbr','center','filename','filesize','exp_library','exp_length','processing','legacy_sample_id', 'anl_title','exp_title','refassem'])

df=df.sort_index(by=['disease','center','participant','sample_id'])
#print df.to_string(index=False,index_names=False)

df.to_csv(output_csv_file)

