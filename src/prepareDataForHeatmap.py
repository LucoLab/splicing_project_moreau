
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
from os import listdir
from os.path import isfile, join
import codecs
import numpy as np
import scipy.stats as stats
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
                #1413527;SRR3069909;EX177770;Mixed lobular ductal;no;ER+/HER2+;T1;N0;100;NA;75+%;75+%;3+;Not done/indeterminate;30 - 39
                if elements[2] not in dictSamples :
                    dictSamples[elements[2]] = { }
                    
                dictSamples[elements[2]]["histo_type"]          = elements[3]
                dictSamples[elements[2]]["response2Treatment"]  = elements[4]
                dictSamples[elements[2]]["molecular_type"]      = elements[5]
                dictSamples[elements[2]]["T.stage"]             = elements[6]
                dictSamples[elements[2]]["N.stage"]             = elements[7]
                
    lines.close()

    return dictSamples

def parse_exons(file) : 
    
    dictSamples = {}
    with open(file) as lines:
            for line in lines: 

                elements = line.strip().split("\t")
                #print(elements)
                #chr1    15529087    15529272    33    0    +    T1    chr1    15529087    15529272    -0.385_DNAJC16_ENSG00000116138_0.86_2    -0.385    +    185
                #chr1    77947452    77947552    584    0    -    T5    chr1    77947452    77947552    0.225_FUBP1_ENSG00000162613_0.49_22    0.225    -    100
                #chr1    77947452    77947552    584    0    -    T1    chr1    77947452    77947547    -0.309_FUBP1_ENSG00000162613_0.49_22    -0.309    -    95
                #chr1    156938417    156938513    684    0    -    T5    chr1    156938417    156938513    0.795_ARHGEF11_ENSG00000132694_0.15_37    0.795    -    96
                #chr1    156938417    156938513    684    0    -    T1    chr1    156938417    156938513    0.481_ARHGEF11_ENSG00000132694_0.15_37    0.481    -    96

                #"ATTENTION CORRECTIF EXON AVANT ds CE.EXON."
                indice_output_whippet = str(int(elements[3]))
                timepointLine         = elements[6]
                nameLine              = elements[10]
                #m            = re.search('^.*BAIAP2.*$', nameLine) # Ok you miss somes of them ENSEMBL_PARY that you count in the sum but that nothing compared to the whole thing
                #if(m) :
                #    print(nameLine)
                #    print(nameLine)
                if indice_output_whippet not in dictSamples : 
                    dictSamples[indice_output_whippet]= { }
                    
                if "timePoint"  not in dictSamples[indice_output_whippet] :
                    dictSamples[indice_output_whippet]["timePoint"] = [ ]
                if "name"  not in dictSamples[indice_output_whippet] :
                    dictSamples[indice_output_whippet]["name"] = [ ]
                
                if "whippet-coord"  not in dictSamples[indice_output_whippet] :
                    dictSamples[indice_output_whippet]["whippet-coord"] = [ ]
                
                if "emt-coord"  not in dictSamples[indice_output_whippet] :
                    dictSamples[indice_output_whippet]["emt-coord"] = [ ]   
                       
                dictSamples[indice_output_whippet]["timePoint"].append(timepointLine)
                dictSamples[indice_output_whippet]["name"].append(nameLine)
                
                dictSamples[indice_output_whippet]["whippet-coord"].append(elements[0]+":"+elements[1]+"-"+elements[2])
                dictSamples[indice_output_whippet]["emt-coord"].append(elements[7]+":"+elements[8]+"-"+elements[9])

                
    lines.close()      
    
    print ("Exons whippet Overlap : "+str(len(dictSamples.keys())))
    
    for indice in sorted(dictSamples.keys()):
        print(indice)
        print(dictSamples[indice]["name"])
        
        dictSamples[indice]["timePoint"] =  test_timepoint(dictSamples[indice])
        
        #test_name_overlap(dictSamples[indice])

        #if (len(dictSamples[indice]["name"])==2) :
        dictSamples[indice]["name"] = reformat(dictSamples[indice])
            
        print(dictSamples[indice]["whippet-coord"])
        print(dictSamples[indice]["emt-coord"]) 
        print(dictSamples[indice]["name"])
        print(dictSamples[indice]["timePoint"])
        
    return dictSamples

def test_name_overlap(dictSampleForOneIndice):
    
    restOfTheNameList =[]
    for name in dictSampleForOneIndice["name"] :
    
        elements_name = name.strip().split("_")
        restOfTheName = "_".join(elements_name[1:len(elements_name)])
        restOfTheNameList.append(restOfTheName)
    
    if(len(list(set(restOfTheNameList)))==1) :
        pass
    else :
        pass
        #return restOfTheNameList[0] 
        #raise AssertionError("Overlap T1 and T5 are not the same, but between T1 and T5 that can be the case between exon can be the same but not flanking", dictSampleForOneIndice)


    
def test_timepoint(dictSampleForOneIndice):
    
    new_timepoint = "timepoint"
    if (len(dictSampleForOneIndice["timePoint"])==1) :
        return "TimePoint: "+str(dictSampleForOneIndice["timePoint"][0])
    elif (len(dictSampleForOneIndice["timePoint"])==2) :
        
        elementsForTimePoint =[]
        if ("T1" in dictSampleForOneIndice["timePoint"] and "T5" in dictSampleForOneIndice["timePoint"]):
            return "TimePoint: Both"
        else : 
            raise AssertionError("Two overlaps are T1 or T5 ", dictSampleForOneIndice)
    else : 
        print(dictSampleForOneIndice)
        raise AssertionError("There is more than two overlaps ", dictSampleForOneIndice)
    
        
        
def reformat (dictSampleForOneIndice) :
    
    new_name = "kikou"
    #['0.256_CS_ENSG00000062485_0.33_3', '0.219_CS_ENSG00000062485_0.
    elementsForNewName =[]
    elementsForNewName2 =[]

    for name in dictSampleForOneIndice["name"] :
        print (name)
        elements_name = name.strip().split("_")
        elementsForNewName.append(elements_name[0])
        elementsForNewName2.append(elements_name[3])

        restOfTheName = "_".join(elements_name[1:len(elements_name)-2])+"_"+elements_name[-1]
        
    new_name = "|".join(elementsForNewName)+"_"+("|".join(list(set(elementsForNewName2))))+"_"+restOfTheName   
    
    return  new_name

def readAllPsis(dictPatients,dictExons,psiFile,idPatient) :
    

    #ensembl_gene_id,hgnc_symbol,gene_biotype,coordinates,strand,event,psi,totalReads,complexity,entropy
    #ENSG00000000003,TSPAN6,protein_coding,chrX:100635558-100635746,-,CE,1,122,K0,0
    count = 0 
    with open(psiFile) as lines:
    #print(psiFile)
    #with codecs.open(psiFile, "r",encoding='utf-8', errors='ignore') as lines:

            for line in lines: 
                if (count == 0):
                    count+=1
                    continue

                elements = line.strip().split(",")
                
                psi     = elements[6]
                bioType = elements[2]
              

                if str(count) in dictExons.keys() :
                    #if "matrice"  not in dictExons[str(count)] :
                        #dictExons[str(count)]["matrice"] = [ ]
                    if idPatient  not in dictExons[str(count)] :
                        dictExons[str(count)][idPatient] = "OK"
                        
                    dictExons[str(count)][idPatient] = psi

                
                    #dictExons[str(count)]["matrice"].append(psi)
                count+=1        
                
    lines.close()
    
    return 
###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################
#ensembl_gene_id,hgnc_symbol,gene_biotype,coordinates,strand,event,psi,totalReads,complexity,entropy
#ENSG00000000003,TSPAN6,protein_coding,chrX:100635558-100635746,-,CE,NA,NA,NA,NA
#ENSG00000000003,TSPAN6,protein_coding,chrX:100635178-100635252,-,CE,NA,NA,NA,NA
#GENE    TpM    Read_Counts
#ENSG00000162851.7    30.2    446.0    
#ENSG00000020633.18    118.39999999999999    3751.9999999999995    
#ENSG00000273004.1    0.0    0.0

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will create conf for all samples.  
    Example : 
    python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/prepareDataForHeatmap.py  -l /home/jean-philippe.villemin/lakitu/PROJECT/BEAUTY/BEAUTY_JP.csv  -d /home/jean-philippe.villemin/lakitu/PROJECT/BEAUTY/output/ -e /home/jean-philippe.villemin/data/data/intersect.80perc.exons.T1.T5.EMT.tsv 

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-l","--listing",action="store",help="Path to all annotation",required=True,type=str,dest='listing')
    parser.add_argument("-d","--dir",action="store",help="Path to dir",required=True,type=str,dest='dir')
    parser.add_argument("-e","--exons",action="store",help="Path to exons",required=True,type=str,dest='exonfile')

    parameters = parser.parse_args()
    
    dictExons = parse_exons(parameters.exonfile)
    
    dictPatients = parse_listing(parameters.listing)
    allSampleFiles = [f for f in listdir(parameters.dir) if (isfile(join(parameters.dir, f)) and f.endswith('CE.psiannoted.csv') and not f.startswith('.')) ]

    columCells = []
    
    result = open(os.path.dirname(parameters.listing)+"/"+"outputBeauty.tsv","w")
    badCovered = open(os.path.dirname(parameters.listing)+"/"+"badCovered.tsv","w")

    print(os.path.dirname(parameters.listing)+"/"+"outputBeauty.tsv")
    for psiFile in allSampleFiles : 
        
        #EX175364.CE.psiannoted.csv
        elements = psiFile.strip().split(".")
        print(elements[0])
        columCells.append(elements[0])
        
        readAllPsis(dictPatients,dictExons,parameters.dir+""+psiFile,elements[0])
    
    # Tcheck Sample too Low
    blackListedPerson = []
    for personID in columCells  : 
        
        events =[]
        for id in dictExons :
      
            events.append(dictExons[id][personID])
                
        total_NA=events.count("NA")

        if (total_NA/len(events) >= 0.5) :
            blackListedPerson.append(personID)
            badCovered.write(personID+"\t"+str(total_NA)+"\n")
    #if (listingToPrint.count("NA"))
    badCovered.close()
    
    patientLine = ["",""]
    histo_typeLine = ["",""]
    mol_typeLine = ["",""]
    t_typeLine = ["",""]
    n_typeLine = ["",""]
    response_typeLine = ["",""]

    for personID in columCells  : 
          #print(personID)
          
          if personID in blackListedPerson : continue
          
          patientLine.append("Patient: "+personID)
          histo_typeLine.append("Histo: "+dictPatients[personID]["histo_type"])
          mol_typeLine.append("MolType: "+dictPatients[personID]["molecular_type"])
          response_typeLine.append("Treatment: "+dictPatients[personID]["response2Treatment"])
          n_typeLine.append("N.stage: "+dictPatients[personID]["N.stage"])
          t_typeLine.append("T.stage: "+dictPatients[personID]["T.stage"])

          #print(dictPatients[personID]["histo_type"])
          #print(dictPatients[personID]["response2Treatment"])
          #print(dictPatients[personID]["molecular_type"])
          #print(dictPatients[personID]["T.stage"])
          #print(dictPatients[personID]["N.stage"])
    result.write("\t".join(patientLine)+"\n")
    result.write("\t".join(histo_typeLine)+"\n")
    result.write("\t".join(response_typeLine)+"\n")
    result.write("\t".join(t_typeLine)+"\n")
    result.write("\t".join(n_typeLine)+"\n")
    result.write("\t".join(mol_typeLine)+"\n")

    for id in dictExons :
        
        listingToPrint =[]
        for personID in columCells  : 
            
            if personID in blackListedPerson : continue

            listingToPrint.append(dictExons[id][personID])
            

        result.write(str(dictExons[id]["name"])+"\t"+str(dictExons[id]["timePoint"])+"\t"+ ("\t".join(listingToPrint))+"\n" )
        
        #print(dictExons[id]["name"])
        #print (dictExons[id]["timePoint"] )
        #print(dictExons[id]["matrice"]) 
            
    #chrY    9810070    9810109    167211    0    +       
    result.close() 
   