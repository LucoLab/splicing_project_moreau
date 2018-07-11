
'''

:date: July 10, 2018
:platform: Debian 

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Filter Matrice of psi 

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
from sklearn.cluster import KMeans
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn import metrics

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
    
    with open(file) as lines:
        for line in lines:
             
            if(count==0):
                count+=1
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
                
            dictSamples[elements[0]]["group"] = elements[1]
          
    lines.close()

    return dictSamples

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

def read_genes(path2file) :
    
    genes =  set ()
       
    with open(path2file) as lines:
        for line in lines: 
       
            genes.add(line.strip())
            
            
    lines.close()
            
    return list(genes)


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
    This script will clean the matrice with all psi.  
    Example : 
    python3 /home/luco/code/python/filterHeatmap.py  -m /home/luco/PROJECT/BEAUTY/output.tsv -g /home/luco/PROJECT/BEAUTY/genes.txt

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-m","--matrice",action="store",help="Matrice with all psi.",required=True,type=str,dest='matrice')
    parser.add_argument("-g","--genes",action="store",help="Genes Symbol you want to look at.",required=False,type=str,dest='genes')

    parameters = parser.parse_args()
    #['ZNF836', 'chr19:52168057-52168152', '1', '1', '1', '1', '1', '1', 'NA', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', 'NA', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', 'NA', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0.96', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', 'NA', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', 'NA', '1', '1', '1', '1', '1', '1', 'NA', '1', '1', '1', '1', '1', '1']

    genes = []
    if (parameters.genes) :
        genes = read_genes(parameters.genes)
    
    result = open(os.path.dirname(parameters.matrice)+"/"+"output_filtered.tsv","w")

    countLine = 0
    with open(parameters.matrice) as lines:
        for line in lines: 
            
            elements = line.strip().split("\t")

            if(countLine in [0,1]):
                result.write(line)
                print(line)
                countLine+=1
                continue
            
            gene = elements[0]
            if (gene in genes ) :
                result.write(line)
            
            '''
            Clustering Part
            '''

            array=np.asarray(elements[2:])
            array = pd.to_numeric(array, errors='coerce')
            array = array[~np.isnan(array)]
            #array= array.reshape(1, -1) #if it contains a single sample
            array= array.reshape(-1, 1) # single feature
            
            if len(array) == 0:
                countLine+=1
                continue

            db = DBSCAN(eps=0.3, min_samples=10).fit(array)
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_
            # Number of clusters in labels, ignoring noise if present.
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
           
            if(n_clusters_>= 2) :
                print(array)
                print(elements[0])
                print(elements[1])
                print('Estimated number of clusters: %d' % n_clusters_)
                result.write(line)

            '''
            End
            '''
  
            countLine+=1
            
        lines.close()      
  
    result.close() 