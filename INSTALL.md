# splicing_project_moreau

---

## 1. Install anaconda to have python3

 wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh
 
 bash Anaconda3-5.2.0-Linux-x86_64.sh
 no add to bashrc
 no visualstudio
 
## 2. Install julia to install whippet next
 
 wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.3-linux-x86_64.tar.gz
 tar -zxvf julia-0.6.3-linux-x86_64.tar.gz
 
Now you can call both soft with these links :

/home/jean-philippe.villemin/julia-d55cadc350/bin/julia
/home/jean-philippe.villemin/anaconda3/bin/python3

## 3. Install whippet into julia

/home/jean-philippe.villemin/julia-d55cadc350/bin/julia

Pkg.update()
Pkg.add("Whippet")
using Whippet