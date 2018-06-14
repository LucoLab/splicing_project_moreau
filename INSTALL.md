# splicing_project_moreau

---

## 1. Install anaconda to have python3

 ```shell
	 wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh
 		bash Anaconda3-5.2.0-Linux-x86_64.sh
  ```
  
 no add to bashrc
 no visualstudio
 
## 2. Install julia to install whippet next
 
 ```shell
 	wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.3-linux-x86_64.tar.gz
 	tar -zxvf julia-0.6.3-linux-x86_64.tar.gz
 ```
 
Now you can call both soft with these links :

/home/jean-philippe.villemin/julia-d55cadc350/bin/julia  
/home/jean-philippe.villemin/anaconda3/bin/python3

## 3. Install whippet into julia

/home/jean-philippe.villemin/julia-d55cadc350/bin/julia

```shell
	Pkg.update()
	Pkg.add("Whippet")
	using Whippet
```

Whippet should be install here ~/.julia/v0.6/Whippet

## 4. Copy file
 
 ```shell
	git clone https://github.com/LucoLab/splicing_project_moreau.git
```
## 5. Annotation GenomeRef

[here](https://www.gencodegenes.org/releases/current.html)

```shell
	mkdir ref
	cd ref
```

```shell
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.primary_assembly.annotation.gtf.gz
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz
```

## 6. Create Index for whippet

```shell
	/home/jean-philippe.villemin/julia-d55cadc350/bin/julia ~/.julia/v0.6/Whippet/bin/whippet-index.jl --fasta /home/jean-philippe.villemin/ref/GRCh38.primary_assembly.genome.fa.gz --gtf /home/jean-philippe.villemin/ref/gencode.v28.primary_assembly.annotation.gtf.gz --suppress-low-tsl -x index_whippet_human
```

## 7. Install R packages

```shell
	R
	source("http://bioconductor.org/biocLite.R")
	biocLite("data.table")
	biocLite("GenomicFeatures")
	biocLite("AnnotationDbi")
	biocLite("biomaRt")
	biocLite("optparse")
```

yes in 
~/R/x86_64-pc-linux-gnu-library/3.1
