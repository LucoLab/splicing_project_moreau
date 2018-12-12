
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
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(config.parameters['path_to_output']+"/"+'activity.log', 'a', 1000000, 1)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    return logger


###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''

    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will launch rmats on fastq.
    Example :
    python3 /home/luco/code/python/splicingWhippetPSI.py --config $1 -i /home/jean-philippe.villemin/index_whippet_human -r /home/jean-philippe.villemin/splicing_project_moreau/Rscript/annotSymbol.R -j /home/jean-philippe.villemin/julia-d55cadc350/bin/julia 
    '''),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')
    parser.add_argument("-i","--index",action="store",help="Path to index directory.",required=True,type=str,dest='index')
    parser.add_argument("-r","--Rpath",action="store",help="Path to index directory.",required=True,type=str,dest='Rpath')
    parser.add_argument("-j","--juliaPath",action="store",help="Path to index directory.",required=True,type=str,dest='juliaPath')

    parameters = parser.parse_args()

    config = custom_parser.Configuration(parameters.file_config,"json")

    logger = create_logger(config)

    logger.info('')
    logger.info('=========> START:')

    logger.info("file_config : "+parameters.file_config)
    logger.info("path_to_output : "+config.parameters['path_to_output'])

    index_path=parameters.index
   

    logger.info("organism : "+config.parameters['organism'])
    logger.info("index_path : "+index_path)
    '''
    ######################################################################
    DataFrame Object Construction : First Part
    ######################################################################
    '''
    for analyse_key,analyse_values in config.parameters["analysis"].items():

        for element_key,element_value in analyse_values.items() :
            logger.debug("-> "+element_key+" :  "+element_value)
            for sample_group in config.parameters["files"][element_value] :
                logger.info("====> Run whippetquant")
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
                            execLine = parameters.juliaPath+" ~/.julia/v0.6/Whippet/bin/whippet-quant.jl -x "+index_path+" "+config.parameters['path_to_input']+hash_values.get("R1")+" -o "+config.parameters['path_to_output']+hash_key
                        else :
                            execLine = parameters.juliaPath+" ~/.julia/v0.6/Whippet/bin/whippet-quant.jl -x "+index_path+"  "+config.parameters['path_to_input']+hash_values.get("R1")+" "+config.parameters['path_to_input']+hash_values.get("R2")+" -o "+config.parameters['path_to_output']+hash_key

                        logger.info(execLine)

                        whippetquant = subprocess.run((execLine),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                        write_subprocess_log(whippetquant,logger)



    logger.info("====> Gunzip, Parse , Annotate")

    hashFiles = {}

    for analyse_key,analyse_values in config.parameters["analysis"].items():
        logger.info(analyse_key)
        logger.debug("TEST :"+analyse_values.get("SAMPLE"))

        samples = []

        for sample_groupTest in config.parameters["files"][analyse_values.get("SAMPLE")] :
            for hash_keyTest,hash_valuesTest in sample_groupTest.items() :
                samples.append(hash_keyTest)

                if(os.path.isfile(config.parameters['path_to_output']+hash_keyTest+".psi")) :
                    logger.info("No need to process again. File "+hash_key+".psi already exist")
                else :

                    logger.info("====> gunzip ")

                    gunzip_command = "gunzip  -k "+config.parameters['path_to_output']+hash_keyTest+".psi.gz"
                    logger.info(gunzip_command)
                    gunzip = subprocess.run((gunzip_command),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(gunzip,logger)

                logger.info(config.parameters['path_to_output']+hash_keyTest+".psi")

                for event in ["CE","AA","AD","RI","AF","AL","TE","TS"] : #

                    if(os.path.isfile(config.parameters['path_to_output']+hash_keyTest+"."+event+".psi")) :
                        continue

                    logger.info("====> parse ")
                    filter_command = config.parameters["path_to_cleaner"]+" "+ event +" "+hash_keyTest+" "+config.parameters['path_to_output']
                    logger.info("EVENT : "+event)
                    logger.info(filter_command)
                    filter = subprocess.run((filter_command),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(filter,logger)
                    logger.info("====> annotate ")
                    Rcommand ="Rscript "+parameters.Rpath+"annotSymbol.R --organism="+config.parameters['organism']+" --file="+config.parameters['path_to_output']+hash_keyTest+"."+event+".psi"
                    logger.info(Rcommand)
                    R = subprocess.run((Rcommand),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(R,logger)
                    logger.info("====> remove Intermediates ")
                    #subprocess.run(("rm "+config.parameters['path_to_output']+hash_keyTest+"."+event+".psi"),shell=True)
                    #subprocess.run(("rm "+config.parameters['path_to_output']+hash_keyTest+".psi"),shell=True)
                    
                    
                    