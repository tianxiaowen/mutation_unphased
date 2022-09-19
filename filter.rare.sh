#!/bin/bash
i=$1
(awk '{if ($4!="A" && $4!="C" && $4!="T" && $4!="G") print $2}' chr${i}/chr${i}.gtstats
 awk '{if ($5!="A" && $5!="C" && $5!="T" && $5!="G") print $2}' chr${i}/chr${i}.gtstats
 awk '{if ($11>=0.01||$11==0) print $2}' chr${i}/chr${i}.gtstats
 cut -f1,2 chr${i}/chr${i}.gtstats | sort -n | uniq -d | awk '{print $2}') | sort | uniq > chr${i}/chr${i}.exclude.for.snp.rare
