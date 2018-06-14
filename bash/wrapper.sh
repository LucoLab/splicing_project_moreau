#!/bin/bash

echo "Here we go..."
echo $1
/home/jean-philippe.villemin/anaconda3/bin/python3 /home/jean-philippe.villemin/splicing_project_moreau/src/splicingWhippetPSI.py -i /home/jean-philippe.villemin/index_whippet_human.jls  -r /home/jean-philippe.villemin/splicing_project_moreau/Rscript/ -j /home/jean-philippe.villemin/julia-d55cadc350/bin/julia  --config $1
