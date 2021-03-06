#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset
#set -o xtrace

# Set magic variables for current file & dir
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
__file="${__dir}/$(basename "${BASH_SOURCE[0]}")"
__base="$(basename ${__file} .sh)"

#################################################################
#
#date: Fev 22, 2017
#platform: Ubuntu 16.04
#author: Villemin Jean-Philippe
#team: Epigenetic Component of Alternative Splicing - IGH
#
# whippet.sh
# Usage : 
# whippet.sh 
# Purpose :
# Filter supposed reliable events.
#
#################################################################


EVENT=$1
FILE=$2
PATH=$3


TYPE=""
echo ${PATH}${FILE}.psi


#################################################################
#Gene	Node	Coord	Strand	Type	Psi	ges psi
#ENSG00000162851.7	1	chr1:246565826-246566324	-	NA
#ENSG00000162851.7	2	chr1:246564346-246564434	-	CE

#CI_Width	CI_Lo,Hi	 Total_Reads	Complexity	Entropy	Inc_Paths	Exc_Paths	Edges
#1.0	0.341	0.652,0.993	6.0	       K0	     0.0     1-2-3:1.0	NA	1-2
#NA	NA	NA	NA	NA	NA	NA	NA	NA

if [ ${EVENT} == "CE" ] || [ ${EVENT} == "SE"  ] ; then

TYPE="CE"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="CE"  ) print ary[1],$3,$4,$5,$6,$9,$10,$11 ;}' ${PATH}${FILE}.psi  > ${PATH}${FILE}.${TYPE}.psi

fi
#################################################################

if [ ${EVENT} == "A3SS" ] || [ ${EVENT} == "AA"  ] ; then

TYPE="AA"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AA"  ) print ary[1],$3,$4,$5,$6,$9,$10,$11 ;}' ${PATH}${FILE}.psi  > ${PATH}${FILE}.${TYPE}.psi


fi

#################################################################

if [ ${EVENT} == "A5SS" ] || [ ${EVENT} == "AD"  ] ; then

TYPE="AD"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AD" ) print ary[1],$3,$4,$5,$6,$9,$10,$11 ;}' ${PATH}${FILE}.psi  > ${PATH}${FILE}.${TYPE}.psi

fi

#################################################################

if [ ${EVENT} == "RI" ] || [ ${EVENT} == "IR"  ] ; then


TYPE="RI"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="RI" ) print ary[1],$3,$4,$5,$6,$9,$10,$11 ;}' ${PATH}${FILE}.psi  > ${PATH}${FILE}.${TYPE}.psi


fi

#################################################################
 
if [ ${EVENT} == "TS" ] ; then

TYPE="TS"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="TS" ) print ary[1],$3,$4,$5,$6,$9,$10,$11 ;}' ${PATH}${FILE}.psi  > ${PATH}${FILE}.${TYPE}.psi


fi
#################################################################
#################################################################
 
if [ ${EVENT} == "TE" ] ; then

TYPE="TE"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="TE" ) print ary[1],$3,$4,$5,$6,$9,$10,$11 ;}' ${PATH}${FILE}.psi  > ${PATH}${FILE}.${TYPE}.psi


fi

#################################################################
 
if [ ${EVENT} == "AF" ] ; then

TYPE="AF"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AF" ) print ary[1],$3,$4,$5,$6,$9,$10,$11 ;}' ${PATH}${FILE}.psi  > ${PATH}${FILE}.${TYPE}.psi


fi
#################################################################
 
if [ ${EVENT} == "AL" ] ; then

TYPE="AL"
/usr/bin/gawk -F "\t"  'BEGIN {OFS="\t";}  {  if ( match($1, "^(\\w+)\\.([0-9]+)", ary) && $5=="AL" ) print ary[1],$3,$4,$5,$6,$9,$10,$11 ;}' ${PATH}${FILE}.psi  > ${PATH}${FILE}.${TYPE}.psi


fi


/bin/sed -i $'1 i\\\ngene\tcoordinates\tstrand\tevent\tpsi\ttotalReads\tcomplexity\tentropy' ${PATH}${FILE}.${TYPE}.psi



