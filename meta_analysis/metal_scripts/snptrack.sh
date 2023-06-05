#! /bin/bash
# Script to generate hg19 rsID from input list of CHR and POS.

java -jar /Volumes/scratch/merrimanlab/CoreTools/CurrentVersions/snptracker/snptracker.jar --merge-file /Volumes/scratch/merrimanlab/CoreTools/CurrentVersions/snptracker/resources/b142_GRCh19.RsMergeArch.bcp.gz --coor-file /Volumes/scratch/merrimanlab/CoreTools/CurrentVersions/snptracker/resources/b142_SNPChrPosOnRef_GRCh19p105.bcp.gz --hist-file /Volumes/scratch/merrimanlab/CoreTools/CurrentVersions/snptracker/resources/b142_GRCh19.SNPHistory.bcp.gz --in ${1} --chr 1 --pos 2 --ref hg19 --by-pos --out ${1%.*} --no-web

awk '{print $1"_"$3, $2}' ${1%.*}.result.txt > ${1%.*}.tmp && mv ${1%.*}.tmp ${1%.*}.result.txt
