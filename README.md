This repository provides resources for mutation rate estimation using the methods developed in the manuscript *Estimating the genome-wide mutation rate from thousands of unrelated individuals* by Tian et al.

The detailed analysis pipeline is outlined in example_analysis_pipeline.md.
Additional scripts called in example_analysis_pipeline.md can be found in the following links.

- IBD detection with [hap-ibd](https://github.com/browning-lab/hap-ibd)
- IBD endpoints estimation with [ibd-ends](https://github.com/browning-lab/ibd-ends)
- Effective population estimtation with [IBDNe](https://faculty.washington.edu/browning/ibdne.html)
- Three-way IBD construction [findibd.R](https://github.com/tianxiaowen/mutation_unphased/blob/main/findibd.R)
- Recoding each 0/1 field in vcf file to represent major/minor alleles [changevcf.py](https://github.com/tianxiaowen/mutation_unphased/blob/main/changevcf.py)
- Get mutation counts [sql_load.py](https://github.com/tianxiaowen/mutation_unphased/blob/main/sql_load.py), [sql_gtcount.py](https://github.com/tianxiaowen/mutation_unphased/blob/main/sql_gtcount.py), [sql_gtcount_split.py](https://github.com/tianxiaowen/mutation_unphased/blob/main/sql_gtcount_split.py)
- Heterozygosity estimation:
- Likelihood calculation for mutation rate [addGC.v3.jar](https://github.com/tianxiaowen/mutation_unphased/blob/main/addGC.v3.jar)
- Get MLE from likelihood outputs
