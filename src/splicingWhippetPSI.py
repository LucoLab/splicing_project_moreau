
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

def write_subprocess_log(completedProcess,logger):
    """
    Write in log the stdout or stderr of subprocess.
    Tcheck if everything was ok.
  
    Args:
        completedProcess (obj): Instance of CompletedProcess send by subprocess.run().
        logger (obj): Instance of logging().
  
    """
    try :
        completedProcess.check_returncode()
        logger.info(completedProcess.stdout)
    except subprocess.CalledProcessError as exc:
                logger.error("===> Exception Caught : ")
                logger.error(exc)  
                logger.error("====> Standard Error : ")
                logger.error(completedProcess.stderr) 
                

def create_logger(config):
    """
    Define a logger instance to write in a file and stream.
  
    Args:
        config (obj): Configuration instance.

    Returns:
        logger (obj): Logger instance to log messages in file and output mainstream.
    
        
    """
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(config.parameters['path_to_output']+"/"+'activity.log', 'a', 1000000, 1)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    logger.addHandler(stream_handler)

    return logger

def create_geneList(bed_file,output) : 
    
    pattern1 ="^0_(.*)_(.*)_(.*)_0$"
    setGenes = set()
    with open(bed_file) as lines:
            for line in lines: 
             
                elements = line.strip().split("\t")
                #chr16    30888713    30888958    0_ENSG00000260083_MIR762HG_0_0    0    +

                #['"ENSG00000000419","DPM1","protein_coding","chr20:50955186-50955285","-","CE",0.2,58,"K0",0']
                found         = re.search(pattern1,  elements[3])
                symbol=str(found.group(2))
                #print("_"+symbol+"_")
                setGenes.add("_"+symbol+"_")
                #print("\t".join([chromosome,str(int(start)-1),end,name,signal,strand])+"\n")
    lines.close()
    f = open(output,'w')
    for gene in setGenes  :
        f.write(str(gene)+"\n")
    f.close()
###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will launch rmats on fastq.  
    Example : 
    python3 /home/luco/code/python/splicingWhippetPSI.py -c /home/luco/PROJECT/BEAUTY/BREAST_CANCER.json

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')

    #parser.add_argument("-a","--analysis",action="store",help="Name of analysis",required=True,type=str,dest='analysis')

    parameters = parser.parse_args()

    config = custom_parser.Configuration(parameters.file_config,"json")
   
    logger = create_logger(config)

    print('=========> START:')
    
    print("file_config : "+parameters.file_config)
    print("path_to_output : "+config.parameters['path_to_output'])

    index_path="None"
    if (config.parameters['organism']=="human") :
        index_path="/home/jean-philippe.villemin/data/data/index_whippet/human/human_index_whippet"
    if (config.parameters['organism']=="mouse") :
        index_path="/home/jean-philippe.villemin/data/data/index_whippet/mouse/mouse_index_whippet"
 
    logger.info("organism : "+config.parameters['organism'])    
    logger.info(index_path)   
    print("    ")
    pattern ="^(.*):(\d+)-(\d+)$"
    '''
    ######################################################################
    DataFrame Object Construction : First Part
    ######################################################################
    '''
    for analyse_key,analyse_values in config.parameters["analysis"].items(): 
       
        #if(analyse_key==parameters.analysis) :
        logger.info(analyse_key)
        for element_key,element_value in analyse_values.items() : 
            logger.info("-> "+element_key+" :  "+element_value)
            for sample_group in config.parameters["files"][element_value] : 
                logger.info("====>")
                logger.info(element_value)
        
                for hash_key,hash_values in sample_group.items() : 
                        
                        # you don't to redo something already done.
                        if(os.path.isfile(config.parameters['path_to_output']+hash_key+".psi.gz") or os.path.isfile(config.parameters['path_to_output']+hash_key+".psi")) : 
                                logger.info("No need to process again. File "+hash_key+".psi.gz already exist")
                                continue
                                
                        logger.info(hash_key)
                        logger.info(hash_values.get("R1"))
                        logger.info(hash_values.get("R2"))
                        if(hash_values.get("R2")=="None" ) : 
                            execLine = "julia /home/jean-philippe.villemin/bin/Whippet.jl/bin/whippet-quant.jl -x "+index_path+" "+config.parameters['path_to_input']+hash_values.get("R1")+" -o "+config.parameters['path_to_output']+hash_key
                        else : 
                            execLine = "julia /home/jean-philippe.villemin/bin/Whippet.jl/bin/whippet-quant.jl -x "+index_path+"  "+config.parameters['path_to_input']+hash_values.get("R1")+" "+config.parameters['path_to_input']+hash_values.get("R2")+" -o "+config.parameters['path_to_output']+hash_key
        
                        logger.info(execLine)
        
                        whippetquant = subprocess.run((execLine),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)                
                        write_subprocess_log(whippetquant,logger)
                            


    logger.info("START TO PARSE ANALYSIS FOR PSI ====>")
    
    hashFiles = {}   

    for analyse_key,analyse_values in config.parameters["analysis"].items(): 
            #if(analyse_key==parameters.analysis) :
        logger.info(analyse_key) 
        logger.info("TEST :    "+analyse_values.get("SAMPLE"))
        
        samples = [] 
        
        for sample_groupTest in config.parameters["files"][analyse_values.get("SAMPLE")] : 
            for hash_keyTest,hash_valuesTest in sample_groupTest.items() : 
                samples.append(hash_keyTest)
                
                if(os.path.isfile(config.parameters['path_to_output']+hash_keyTest+".psi")) : 
                    logger.info("No need to process again. File "+hash_key+".psi.gz already exist")
                else :
                    
                    logger.info("GUNZIP ====>")
                    gunzip_command = "gunzip  -k "+config.parameters['path_to_output']+hash_keyTest+".psi.gz"
                    logger.info(gunzip_command)
                    gunzip = subprocess.run((gunzip_command),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(gunzip,logger)

                logger.info(config.parameters['path_to_output']+hash_keyTest+".psi") 

                for event in ["CE"] : #,"AA","AD","RI"

                    filter_command = config.parameters["path_to_cleaner"]+" "+ event +" "+hash_keyTest+" "+config.parameters['path_to_output']
                    logger.info("EVENT : "+event)
                    logger.info(filter_command)
                    filter = subprocess.run((filter_command),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(filter,logger)
                    
                    Rcommand ="Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/annotSyymbol.R --organism="+config.parameters['organism']+" --file="+config.parameters['path_to_output']+hash_keyTest+".20."+event+".psi"
                    logger.info(Rcommand)
                    R = subprocess.run((Rcommand),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(R,logger)
                    
                    Rcommand2 ="Rscript /home/jean-philippe.villemin/code/RNA-SEQ/Rscript/annotSyymbol.R --organism="+config.parameters['organism']+" --file="+config.parameters['path_to_output']+hash_keyTest+".80."+event+".psi"
                    logger.info(Rcommand2)
                    R2 = subprocess.run((Rcommand2),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(R2,logger)
                    
                    
                hashFiles[hash_keyTest] = {}   
                for event in ["CE"] :#,"AA","AD","RI"
                    
                    hashFiles[hash_keyTest][event] = {}
                    
                    for quality in ["20","80"] :
                        
                        hashFiles[hash_keyTest][event][quality] = {}
                        hashFiles[hash_keyTest][event][quality] = config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".sorted.bed"
                        
                        count = 0
                        logger.info(config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".psiannoted.csv")

                        f = open(config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".bed",'w')
                        countLines = 0
                        
                        with open(config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".psiannoted.csv") as lines:
                            for line in lines: 
                                if (count ==0 ):
                                    count=1 
                                    continue
                                elements = line.strip().split(",")

                                #['"ENSG00000000419","DPM1","protein_coding","chr20:50955186-50955285","-","CE",0.2,58,"K0",0']
                                lineCoordinates         = re.search(pattern,  elements[3])
                                chromosome      =   lineCoordinates.group(1)
                                start           =   lineCoordinates.group(2)
                                end             =   lineCoordinates.group(3)
                                strand          =   elements[4]
                                name            = "0_"+elements[0]+"_"+elements[1]+"_"+elements[6]+"_0"
                                signal          =   elements[6]
                                
                                #print("\t".join([chromosome,str(int(start)-1),end,name,signal,strand])+"\n")
                                f.write("\t".join([chromosome,str(int(start)-1),end,name,signal,strand])+"\n")
                                countLines+=1
                            lines.close()
                        f.close()
                        
                        print("number of lines  : "+str(countLines))

                        proc_sort = subprocess.run(("sort -k1,1V -k2,2n  "+config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".bed"+ " > "+config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".sorted.bed"),shell=True, check=True)

                        subprocess.run(("rm "+config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".psiannoted.csv"),shell=True)
                        #subprocess.run(("rm "+config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".bed"),shell=True)
                        subprocess.run(("rm "+config.parameters['path_to_output']+hash_keyTest+"."+quality+"."+event+".psi"),shell=True)
            print("")

        path_to_irimia_always_inclus="/home/jean-philippe.villemin/data/data/raw.irimia.108samples.inclus.hg38.bed"
        for event in ["CE"] :#,"AA","AD","RI"
                  
            for quality in ["20","80"] :

                base=[]
                for index_sample in range(1,len(samples))  : 
                    base.append(hashFiles[samples[index_sample]][event][quality]) 

                proc_intersect = subprocess.run(("bedtools intersect -u -f 1 -r  -a "+hashFiles[samples[0]][event][quality]+" -b "+",".join(base)+ " > "+config.parameters['path_to_output']+event+".intersect."+analyse_key+"."+quality+".bed"),shell=True, check=True) 
                write_subprocess_log(proc_intersect,logger)
                
                proc_intersect = subprocess.run(("bedtools intersect -u -a "+config.parameters['path_to_output']+event+".intersect."+analyse_key+"."+quality+".bed -b "+path_to_irimia_always_inclus+ " > "+config.parameters['path_to_output']+event+".intersect.Constitutif."+analyse_key+"."+quality+".bed"),shell=True, check=True) 
                write_subprocess_log(proc_intersect,logger)
            
                proc_intersect = subprocess.run(("bedtools intersect -v  -a "+config.parameters['path_to_output']+event+".intersect."+analyse_key+"."+quality+".bed -b "+path_to_irimia_always_inclus+ " > "+config.parameters['path_to_output']+event+".intersect.notConstitutif."+analyse_key+"."+quality+".bed"),shell=True, check=True) 
                write_subprocess_log(proc_intersect,logger)
                
                ##################################################################################################

                
                proc_count0 = subprocess.getoutput(("wc -l  "+config.parameters['path_to_output']+event+".intersect."+analyse_key+"."+quality+".bed"))
                logger.info("Lines All "+event+" "+quality+" : "+str(proc_count0.split(" ")[0]))
                
                proc_count2 = subprocess.getoutput(("wc -l  "+config.parameters['path_to_output']+event+".intersect.notConstitutif."+analyse_key+"."+quality+".bed"))
                logger.info("Lines notConstitutif "+event+" "+quality+" :  "+str(proc_count2.split(" ")[0]))
                
                proc_count = subprocess.getoutput(("wc -l  "+config.parameters['path_to_output']+event+".intersect.Constitutif."+analyse_key+"."+quality+".bed"))
                logger.info("Lines Constitutif "+event+" "+quality+" : "+str(proc_count.split(" ")[0]))
               
                ##################################################################################################

                pathToGeneListnotConstitutif=config.parameters['path_to_output']+event+".intersect.notConstitutif."+analyse_key+"."+quality+".GENES.txt"
                create_geneList(config.parameters['path_to_output']+event+".intersect.notConstitutif."+analyse_key+"."+quality+".bed",pathToGeneListnotConstitutif) 
                
                # I want constitutif exons from genes found with Constitutif exons (irimia)
                proc_grep = subprocess.getoutput(("grep -F -f "+pathToGeneListnotConstitutif+" "+config.parameters['path_to_output']+event+".intersect.Constitutif."+analyse_key+"."+quality+".bed")) 
                
                file = open(config.parameters['path_to_output']+event+".intersect.Constitutif.butsameGenesAsNotConstitutif."+analyse_key+"."+quality+".bed",'w')
                file.write(proc_grep)
                file.close()
                
                proc_count3 = subprocess.getoutput(("wc -l  "+config.parameters['path_to_output']+event+".intersect.Constitutif.butsameGenesAsNotConstitutif."+analyse_key+"."+quality+".bed"))
                logger.info("Lines intersect.Constitutif.butsameGenesAsNotConstitutif "+event+" "+quality+" : "+str(proc_count3.split(" ")[0]))
            
                #proc_merge = subprocess.run(("bedtools merge -i "+parameters.path+parameters.out+".final.sorted.bed "+ "  -c 4,5,6 -o collapse,mean,collapse > "+parameters.path+parameters.out+"_TMP_final.sorted.merged.bed"),shell=True, check=True) 
                #write_subprocess_log(proc_merge,logger)
                
    