#!/bin/bash

i=$1
j=$2
ibdfile=$3
mapfile=$4
idfile=$5

echo "chr=$i"
echo "split=$j"

echo "copy database"
mkdir /scratch/chr${i}_split${j}
cp chr${i}/chr${i}.maf0.01snp.gt.recode.db /scratch/chr${i}_split${j}
cp chr${i}/chr${i}.raresnp.gt.recode.db /scratch/chr${i}_split${j}
export PATH=/python/anaconda2/bin:$PATH

echo "get mutation count"
python sql_gtcount_split.py ibdtrios/chr${i}.maf0.25.2.5-6cM.txt /scratch/chr${i}_split${j}/chr${i}.raresnp.gt.recode.db /scratch/chr${i}_split${j}/chr${i}.maf0.01snp.gt.recode.db ${ibdfile} ${mapfile} ${idfile} $j > chr${i}/triosgt.maf0.25.2.5-6cM.${j}.txt

echo "clear repository in the end"
rm -r /scratch/chr${i}_split${j}  
echo ""