# splicing_project_moreau

---


Scripts available here are in Python3.
It's not required but advised to install Conda if python3 is not set up on your computer.
It will make things easier then for installing tools or switching to older python version if needed.

Look to INSTALL, every steps are documented.
First you need to install Whippet [here](https://github.com/timbitz/Whippet.jl)  

Here you will find the project manager [here](https://trello.com/b/XFuccCgE/splicingprojectcolab)
Just in case if we want a clean follow-up of what we are doing.

## 1. Create json config for each fastq pair


Each time, you need to change python3 path accordingly to your home user path.

---

It creates all the json files used by the script.
The json files are used to define path to input,output directories and name of each samples.
The script will create a file called listing.json.per.sample.txt with all links to json configs to process sequencially in the next step.

```shell
     /home/jean-philippe.villemin/anaconda3/bin/python3 /home/jean-philippe.villemin/splicing_project_moreau/src/createJsonConfForPsi.py \
     -l /home/jean-philippe.villemin/samples.txt \
     -i /eqmoreaux/data/global_data/MMC/RNAseq/FASTQ/ \
     -o /home/jean-philippe.villemin/moreau_splicing/ \
     -s /home/jean-philippe.villemin/splicing_project_moreau/bash/whippet_filter_eventType_wrapped_for_psiOnly.sh 

```
The output directory is set to moreau_splicing. Inside you will find _config_ dir with json files and _output_ dir where splicing results will be written. 


## 2. Process sequentially for splicing analysis

---
```shell
	cat listing.json.per.sample.txt| xargs -n 1  -I %  wrapper.sh % 
```

It will read each line of listing.txt to access each json config.
It's a loop in fact...  

Wrapper.sh calls the main script as follows using the config path received :

###### VERY IMPORTANT :

Now you need to set paramerers -j, -i and -r.
-j : is the path to julia executable
-i : path to whippet jls index
-r : path to dir where you find annotSymbol.R

So you need to modify accordingly to your own conf.

You can test for one sample if it works.

```shell
nohup /home/jean-philippe.villemin/anaconda3/bin/python3 /home/jean-philippe.villemin/splicing_project_moreau/src/splicingWhippetPSI.py --config /home/jean-philippe.villemin/moreau_splicing/config/E10061.json  -i /home/jean-philippe.villemin/index_whippet_human.jls  -r /home/jean-philippe.villemin/splicing_project_moreau/Rscript/ -j /home/jean-philippe.villemin/julia-d55cadc350/bin/julia &
```
Then for all

```shell
nohup cat listing.json.per.sample.test.txt | xargs -n 1  -I %  ../splicing_project_moreau/bash/wrapper.sh %  & 
```

Read whippet output to get an idea of what we have at the end.

Here I show you a screenshot of what we have in the directory output for one experiment treated.

![alt text](https://github.com/LucoLab/splicing_project_moreau/blob/master/img/main_output.png "Outputs")

## 3. Create PSI Matrice by concatening all Samples 

This will create a matrice named output.tsv.

```shell
/home/luco/localLib/anaconda3/bin/python3 /home/luco/code/python/prepareDataForHeatmap.py  -l /home/luco/PROJECT/BEAUTY/BEAUTY_ID2GROUP.tsv -d /home/luco/PROJECT/BEAUTY/output/ -e /home/luco/PROJECT/TCGA/TCGA/CE.whippet.tsv
```

Examples are given in config dir for the files used by the script.
_BEAUTY_ID2GROUP.tsv_ is just a two columns file with nameOfPatient & group.  

*I strongly advise to regenerate CE.whippet.tsv yourself or at least check it correspond to file.CE.psiannoted.csv (same number of lines, orderded the same way...)* It depends of the version gtf you used. Exons can be ordered differently, removed or added.

Take any file.CE.psiannoted.csv  and apply the following command to regenerate it properly from your dataset.

```shell
/usr/bin/gawk -F ","  'BEGIN {OFS="\t";}  {  if ( match($4, "^(chr.*):([0-9]+)-([0-9]+)", ary) ) print ary[1],ary[2]-1,ary[3],NR-1,0,$5 ;}' /home/luco/PROJECT/BEAUTY/output/EX173639.CE.psiannoted.csv 
```

## 4. Apply filters to this matrice

Now you get the whole matrice, you need to apply fitlters.
You can do it by yourself with your own scripts or I provide some tools that will create a output_filtered.tsv file

```shell
python3 /home/luco/code/python/filterHeatmap.py  -m /home/luco/PROJECT/BEAUTY/output.tsv -g /home/luco/PROJECT/BEAUTY/genes.txt
 ```

This script work with a provided list of genes symbol. It also use a cluster algorithm (scikit-learn : DBSCAN, harcoded eps=0.3, min_samples  10) to detect some potential splicing changes between groups.  So your are not force to use a gene lists. It's optional but recommanded because the cluster algorithm is not optimal yet.  


## 5. Use Morpheus [here](https://software.broadinstitute.org/morpheus/) to visualise the matrice in a heatmap.

---

At the end a quick way to visualize your matrice is this tool.
Try it, play with it. It's cool.


