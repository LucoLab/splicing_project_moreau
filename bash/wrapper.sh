#!/bin/bash

/home/jean-philippe.villemin/anaconda3/bin/python3 /home/luco/code/python/splicingWhippetPSI.py --config $1 -i /home/jean-philippe.villemin/index_whippet_human.jls  -r /home/jean-philippe.villemin/splicing_project_moreau/Rscript/annotSymbol.R -j /home/jean-philippe.villemin/julia-d55cadc350/bin/julia 