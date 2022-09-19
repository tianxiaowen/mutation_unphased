#!/bin/bash
i=$1
Nhaplotype=$2
export PATH=/python/anaconda2/bin:$PATH 
echo "python sql_load.py chr${i}/chr${i}.maf0.01snp.gt.recode.gz chr${i}/chr${i}.maf0.01snp.gtstats chr${i}/chr${i}.maf0.01snp.gt.recode.db $Nhaplotype"
echo "python sql_load.py chr${i}/chr${i}.raresnp.gt.recode.gz chr${i}/chr${i}.raresnp.gtstats chr${i}/chr${i}.raresnp.gt.recode.db $Nhaplotype"