
'''

:date: September 9, 2017
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Launch Splicing with Whippet

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import os 
import re

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
                elements = line.strip().split(";")
                #print(elements)
                #dictSamples[elements[2]]= ["epurated_"+elements[1]+"_1_val_1.fq.gz" ,"epurated_"+elements[1]+"_2_val_2.fq.gz"]
                dictSamples[elements[2]]= ["epurated_"+elements[1]+"_1.fastq" ,"epurated_"+elements[1]+"_2.fastq"]

    lines.close()

    return dictSamples
###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will create conf for all samples.  
    Example : 
    python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/createJsonConfForPsi.py -l /home/jean-philippe.villemin/lakitu/PROJECT/BEAUTY/BEAUTY_JP.csv 

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-l","--listing",action="store",help="Path to a listing",required=True,type=str,dest='listing')

    parameters = parser.parse_args()
    
    dictSamples = parse_listing(parameters.listing)
    count=0
    mywholeList = open("/home/jean-philippe.villemin/lakitu/PROJECT/BEAUTY/listing.conf.per.sample.txt","w")
    for sample in dictSamples :
        count+=1
        lineToPrint ='{"'+sample+'" : {  "R1": "'+dictSamples[sample][0]+'","R2": "'+dictSamples[sample][1]+'"}},'
        pathToConf="/home/jean-philippe.villemin/lakitu/PROJECT/BEAUTY/conf_per_sample/"+sample+".json"
        pathToConfShenron="/home/luco/PROJECT/BEAUTY/conf_per_sample/"+sample+".json"
        mywholeList.write(pathToConfShenron+"\n")
        
        sampleFile=open(pathToConf,"w")

        sampleFile.write("{\n")
        sampleFile.write('"organism":"human",\n')    
        sampleFile.write('"path_to_output" : "/home/luco/PROJECT/BEAUTY/output/",\n')
        sampleFile.write('"path_to_input" : "/home/luco/commun/BEAUTY/dataClean/",\n')
        sampleFile.write('"path_to_cleaner" : "/home/luco/code/bash/whippet_filter_eventType_wrapped_for_psiOnly.sh",\n')
        sampleFile.write('"files" :\n')
        sampleFile.write("{\n")
        sampleFile.write('"SAMPLE1" : [\n')
        sampleFile.write(lineToPrint[:-1]+"\n")
        sampleFile.write("]\n")
        sampleFile.write("},\n")
        sampleFile.write('"analysis":{\n')
        sampleFile.write('"BREAST_CANCER_BEAUTY": { "SAMPLE": "SAMPLE1"}\n')
        sampleFile.write("}\n")
        sampleFile.write("}\n")
        
    sampleFile.close()    
    mywholeList.close()


    


    