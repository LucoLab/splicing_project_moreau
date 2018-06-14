
'''

:date: June 14, 2018
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Launch Creation of the json configuration files

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import os 
import re
import pathlib

###########################################################################################################
########################################   Functions   ####################################################
###########################################################################################################

def parse_listing(file) : 
    
    dictSamples = {}
    count=0
    with open(file) as lines:
            for line in lines: 
                if(count==0):
                    count+=1
                    continue
                elements = line.strip().split("\t")
                #print(elements)
                #dictSamples[elements[2]]= ["epurated_"+elements[1]+"_1_val_1.fq.gz" ,"epurated_"+elements[1]+"_2_val_2.fq.gz"]
                dictSamples[elements[0].replace("\"","")]= [elements[0].replace("\"","")+"_R1.fastq.gz" ,elements[0].replace("\"","")+"_R2.fastq.gz"]

    lines.close()

    return dictSamples
###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will create conf for all samples.  
    /home/jean-philippe.villemin/anaconda3/bin/python3 /home/jean-philippe.villemin/splicing_project_moreau/src/createJsonConfForPsi.py  -l /home/jean-philippe.villemin/samples.txt -i /eqmoreaux/data/global_data/MMC/RNAseq/FASTQ/ -s /home/jean-philippe.villemin/splicing_project_moreau/bash/whippet_filter_eventType_wrapped_for_psiOnly.sh -o /home/jean-philippe.villemin/moreau_splicing/
    Example : 
    /home/jean-philippe.villemin/anaconda3/bin/python3 ./src/createJsonConfForPsi.py \
    -l ../samples.csv \
    -i /eqmoreaux/data/global_data/MMC/RNAseq/FASTQ/ \
    -s /home/jean-philippe.villemin/splicing_project_moreau/bash/whippet_filter_eventType_wrapped_for_psiOnly.sh
    -o /home/jean-philippe.villemin/moreau_splicing/

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-l","--listing",action="store",help="Path to a listing",required=True,type=str,dest='listing')
    parser.add_argument("-s","--scriptBash",action="store",help="Path to another bash script",required=True,type=str,dest='scriptBash')
    parser.add_argument("-o","--outputDir",action="store",help="Path to directory output",required=True,type=str,dest='outputDir')
    parser.add_argument("-i","--inputDir",action="store",help="Path to a listing",required=True,type=str,dest='inputDir')

    parameters = parser.parse_args()
    
    config_for_all_json = parameters.outputDir+"/config/"
    output = parameters.outputDir+"/output/"

    pathlib.Path(config_for_all_json).mkdir(parents=True, exist_ok=True) 
    pathlib.Path(output).mkdir(parents=True, exist_ok=True) 

    mywholeList = open(parameters.outputDir+"listing.json.per.sample.txt","w")

    dictSamples = parse_listing(parameters.listing)
    
    count=0
    
    for sample in dictSamples :
        count+=1
        lineToPrint ='{"'+sample+'" : {  "R1": "'+dictSamples[sample][0]+'","R2": "'+dictSamples[sample][1]+'"}},'
        pathToConfJson=config_for_all_json+sample+".json"
        
        mywholeList.write(pathToConfJson+"\n")
        
        sampleFile=open(pathToConfJson,"w")
        sampleFile.write("{\n")
        sampleFile.write('"organism":"human",\n')    
        sampleFile.write('"path_to_output" : "'+output+'",\n')
        sampleFile.write('"path_to_input" : "'+parameters.inputDir+'",\n')
        sampleFile.write('"path_to_cleaner" : "'+parameters.scriptBash+'",\n')
        sampleFile.write('"files" :\n')
        sampleFile.write("{\n")
        sampleFile.write('"SAMPLE1" : [\n')
        sampleFile.write(lineToPrint[:-1]+"\n")
        sampleFile.write("]\n")
        sampleFile.write("},\n")
        sampleFile.write('"analysis":{\n')
        sampleFile.write('"PROJECT": { "SAMPLE": "SAMPLE1"}\n')
        sampleFile.write("}\n")
        sampleFile.write("}\n")
        
    sampleFile.close()    
    mywholeList.close()


    


    