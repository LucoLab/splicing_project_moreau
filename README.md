# splicing_project_moreau

---

Scripts available here are in Python3.
It's not required but advised to install Conda if python3 is not set up on your computer.
It will make things easier then for installing tools or switching to older python version if needed.

First you need to install Whippet [here](https://github.com/timbitz/Whippet.jl)  

Here you will find the project manager [here](https://trello.com/b/XFuccCgE/splicingprojectcolab)
Just in case if we want a clean follow-up of what we are doing.

## 1. Create json config for each fastq pair

---

```shell
    python3 createJsonConfForPsi.py -l samples.csv 
```

The script need to be changed depending on samples.csv format.  
It should be something like ID_Project,FASQT1_NAME,FASTQ2_NAME.
Hardcoded paths in the script need to be modified also.
SAMPLE1 is a keyword. Do not change.
example1 must be replace by the name of your sample.
You need to set path to output and input directories.

## 2. Create whole listing of json configs

---

If you set all your json configs in a dir. Use a command like 'find' to print all the paths in file to use then in the nex step.

```shell
	find path2directory -name *.json > listing.txt
```

So you will have something like following : 

pathToConf/example1.json
pathToConf/example2.json

## 3. Process sequentially for splicing analysis

---
```shell
	cat listing.txt | xargs -n 1  -I %  wrapper.sh 
```

It will read each line of listing.txt to access each json config.

Wrapper.sh call the main script as follows using the config path received :

```shell
python3 splicingWhippetPSI.py --config $1
```

Read whippet output to get an idea of what we have at the end.

Here I show you a screenshot of what we have in the dir for experiment treated.

![alt text](https://github.com/LucoLab/splicing_project_moreau/blob/master/img/main_output.png "Outputs")


## 4. Use Morpheus [here](https://software.broadinstitute.org/morpheus/) to visualise the matrice in a heatmap.

---

I was using this script to parse all outputs to grap the psi values from different experiments an concat them in one matrice. 
Also I was using an another bed to get a subset of the matrice for exon of interest.
Need to be modified.

```shell
    prepareDataForHeatmap.py  -l listing.csv  -d /home/jean-philippe.villemin/lakitu/PROJECT/BEAUTY/output/ -e specific.bed.tsv
```
