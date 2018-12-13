
'''

:date: July 10, 2018
:platform: Debian 

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Create Matrice psi from several samples

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import os 
import re
from os import listdir
from os.path import isfile, join
import codecs
import numpy as np
import scipy.stats as stats

###########################################################################################################
########################################   Functions   ####################################################
###########################################################################################################

def parse_listing(file) : 
    """
    Parse list of patient and group. Two column files.
  
    Args:
        file (string): Path to file.

    Returns:
        dictSamples (dict): Dict samples to group. 'EX167687': {'group': 'Triple Negative'}, 'EX167686': {'group': 'Lum A'},

    """
    
    dictSamples = {}
    count       =  0
    list_groups = []
    
    with open(file) as lines:
        for line in lines:
             
            if(count==0):
                count+=1
               
                list_groups = line.strip().split("\t")
                #print(list_groups)
                continue
            
            elements = line.strip().split("\t")
            #bcr_patient_barcode    Call
            #TCGA-A8-A08F    LumB
            #TCGA-A8-A09K    LumA
            
            #Subject_ID    Clinical.Molecular.Subtype
            #EX181420    Triple Negative
            #EX181336    Lum B
            #EX181261    Lum Unk
            
            if elements[0] not in dictSamples :
                dictSamples[elements[0]] = { }
            
            for name_group in list_groups :
                #print(name_group)
                #print(list_groups.index(name_group))
                #print(elements)
                dictSamples[elements[0]][name_group] = elements[list_groups.index(name_group)]
          
    lines.close()

    return dictSamples,list_groups

def parse_exons(file) : 
    """
    Parse list of exons used by whippet to determine exons skipping CE.
    CE.whippet.bed is always ordered the same way but can be different depending annotation gtf file used.

    Args:
        file (string): Path to file CE.whippet.bed.

    Returns:
        dictExons (dict): Dict exons to list of features.

        '115345': {'whippet-coord': ['chr12:52171657-52171713'], 'gene-symbol': ['KRT80']}, 
        '115346': {'whippet-coord': ['chr12:52678541-52678756'], 'gene-symbol': ['KRT1']},
    """ 
    dictExons = {}
    with open(file) as lines:
        for line in lines: 
        
            elements = line.strip().split("\t")
            #chr17    20426790    20427079    139766    0    + DDX5
        
            indice_output_whippet = str(int(elements[3]))
        
            if indice_output_whippet not in dictExons : 
                dictExons[indice_output_whippet]= { }
        
            if "whippet-coord"  not in dictExons[indice_output_whippet] :
                dictExons[indice_output_whippet]["whippet-coord"] = [ ]
            if "gene-symbol"  not in dictExons[indice_output_whippet] :
                dictExons[indice_output_whippet]["gene-symbol"] = [ ]
        
            dictExons[indice_output_whippet]["gene-symbol"].append(elements[6])
            dictExons[indice_output_whippet]["whippet-coord"].append(elements[0]+":"+elements[1]+"-"+elements[2])

        lines.close()      
  
    return dictExons

def readAllPsis(dictExons,dictPatients,psiFile,idPatient) :
    """
    Grep the psi values from psi files.
  
    Args:
        dictPatients (dict): dictionnary, patients to group.
        psiFile (string): path to psi file to parse.
        idPatient (string): idPatient , also name of the file to parse.

    """

    #ensembl_gene_id,hgnc_symbol,gene_biotype,coordinates,strand,event,psi,totalReads,complexity,entropy
    #ENSG00000000003,TSPAN6,protein_coding,chrX:100635558-100635746,-,CE,1,122,K0,0
    
    indiceExonWhippet = 0 
    with open(psiFile) as lines:

        for line in lines: 
            if (indiceExonWhippet == 0):
                indiceExonWhippet+=1
                continue

            elements = line.strip().split(",")
            
            psi     = elements[6]
            dictPatients[idPatient][str(indiceExonWhippet)]= psi
            

            indiceExonWhippet+=1        
                
    lines.close()
    
    return 

###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

#ensembl_gene_id,hgnc_symbol,gene_biotype,coordinates,strand,event,psi,totalReads,complexity,entropy
#ENSG00000000003,TSPAN6,protein_coding,chrX:100635558-100635746,-,CE,NA,NA,NA,NA

#GENE    TpM    Read_Counts
#ENSG00000162851.7    30.2    446.0    


if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will create a clean matrice with all psi.  
    Example : 
    /home/jean-philippe.villemin/anaconda3/bin/python3 /home/jean-philippe.villemin/splicing_project_moreau/src/prepareDataForHeatmap.py  -l /home/jean-philippe.villemin/moreau_splicing_test/ID_2_GROUP.tsv -d /home/jean-philippe.villemin/moreau_splicing_test/output/ -e /home/jean-philippe.villemin/CE.whippet.tsv -t CE
        
    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-l","--listing",action="store",help="Path to all annotation",required=True,type=str,dest='listing')
    parser.add_argument("-d","--dir",action="store",help="Path to dir",required=True,type=str,dest='dir')
    parser.add_argument("-e","--exons",action="store",help="Path to exons",required=True,type=str,dest='exonfile')
    parser.add_argument("-t","--type",action="store",help="TYPE_EVENT",required=True,type=str,dest='type')

    parameters = parser.parse_args()
    
    dictPatients,list_groups = parse_listing(parameters.listing)

    #/home/luco/PROJECT/TCGA/TCGA/CE.whippet.tsv
    dictExons    = parse_exons(parameters.exonfile)

    
    allSampleFiles = [f for f in listdir(parameters.dir) if (isfile(join(parameters.dir, f)) and f.endswith(parameters.type+'.psiannoted.csv') and not f.startswith('.')) ]

    columCells = []
    
    result = open(os.path.dirname(parameters.listing)+"/"+parameters.type+".output.tsv","w")
    

    for psiFile in allSampleFiles : 
        
        elements = psiFile.strip().split(".")
        
        print(elements[0])

        columCells.append(elements[0])
        
        readAllPsis(dictExons,dictPatients,parameters.dir+""+psiFile,elements[0])

    
    
    list_several_group = [[] for x in range(0,len(list_groups))]

    count_i=0
    for name_group in list_groups :
            list_several_group[count_i].append("")
            list_several_group[count_i].append("")
            count_i+=1


    for personID in columCells  : 
     
        count=0
        for name_group in list_groups :
            list_several_group[count].append(dictPatients[personID][name_group])
            count+=1

    for group_listed in list_several_group :   
        #print(group_listed)    
        result.write("\t".join(group_listed)+"\n")

    for indiceWhippet in sorted(dictExons.keys()) :
        
        listingToPrint = []
        for personID in columCells  : 
            #listingToPrint.append(dictExons[indiceExonWhippet][personID])
            listingToPrint.append(dictPatients[personID][indiceWhippet])

        result.write(str(dictExons[indiceWhippet]["gene-symbol"][0])+"\t"+str(dictExons[indiceWhippet]["whippet-coord"][0])+ "\t" + ("\t".join(listingToPrint))+"\n" )

    result.close() 
   