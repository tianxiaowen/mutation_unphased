# Mutation rate estimation

This is an example pipeline for estimating mutation rate from phased genotype data using the methods developed in the manuscript *Estimating the genome-wide mutation rate from thousands of unrelated individuals* by Tian et al. The pipeline includes the following steps:   
1. IBD detection     
2. Estimate effective population size      
3. Search for three-way IBD sharing       
4. Get mutation count       
5. Estimate mutation rate        

## STEP 1: IBD detection
We use hap-ibd and ibd-ends to detect pairwise IBD segments from the phased haplotypes
When detecting IBD segments using hap-ibd, we consider only markers with MAF $\geq$ 0.25 by setting the `minmac` parameter of hap-ibd to 
$$minmac=2N*0.25$$ 
where $N$ is the sample size. 

We used ibd-ends for a second step of IBD segment inference, which estimates the posterior medians of the endpoints of IBD segments as the adjusted endpoints of the IBD segments. 

```      
mkdir ibdseg     
for i in `seq 1 22`; do       
    java -jar hap-ibd.jar gt=chr${i}/chr${i}.vcf.gz map=${map_file} min-seed=0.5 min-extend=0.1 max-gap=5000 min-output=1 min-mac=${minmac} out=ibdseg/chr${i}.maf.0.25 nthreads=12
done      
      
mkdir ibdends     
for i in `seq 1 22`; do     
    java -jar -Xmx180g ibd-ends.jar gt=chr${i}/chr${i}.vcf.gz ibd=ibdseg/chr${i}.maf.0.25.ibd.gz map=${map_file} out=ibdends/chr${i}.maf.0.25.ibdends nsamples=1     
done     

for i in `seq 1 22`; do     
    zcat ibdends/chr${i}.maf.0.25.ibdends.ibd.gz | cut -f1-5,10-12 | sed '1 d' | gzip > ibdends/chr${i}_maf.0.25.ibd.gz     
done     
```     
As a final step of IBD detection, from the output IBD segments by ibd-ends, we need to removed any IBD segments shared between pairs of individuals who are duplicated samples, monozygotic twins, or parent and offspring. 

## STEP 2: Estimate effective population size
We inferred the recent effective populations sizes of the study population using the default settings of IBDNe with IBD segments output by ibd-ends.
```     
zcat ibdends/chr*_maf.0.25.ibd.gz | java -Xmx90g -jar ibdne.23Apr20.ae9.jar map=${map_file} out=${IBDNe_output} nthreads=12     
```     

## STEP 3. Search for three-way IBD sharing
When searching three-way IBD sharing, we restricted to only IBD segments with length between 2.5 cM to 6 cM. We use IBD segments output by ibd-ends excluding any first-degree relatives. We provided an R program findibd.R to search for three-way IBD sharing. This program takes the following arguments:    

- `chr`: chromosome number 1-22
- `ibdname`: path to the IBD file
- `minlen`: minimum length of the three-way IBD sharing
- `outname`: path to save the output file 
- `idfile`: text file that has the ID names of samples in the VCF file as a single column
- `mapfile`: path to the genetic maps, examples can be found https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/

```
mkdir ibdtrios     
for i in `seq 1 22`; do     
    Rscript findibd.R $i ${ibdname} 2.5 ibdtrios/chr${i}.maf0.25.2.5-6cM.txt ${idfile} ${mapfile}     
done     
```

## STEP 4. Get mutation count
To prepare data for getting mutation counts, the following steps are needed: 

- separate VCF files to data on SNP variants with MAF $\geq 1\%$ and data on rare SNP variants with MAF $< 1\%$
- recode 0/1 field in the VCF files so that 0/1 represents major/minor allele instead of REF/ALT

To compute MAF on each marker, we can use the software [`gtstats.jar`](https://faculty.washington.edu/browning/beagle_utilities/utilities.html#gtstats) to compute genotype statistics on VCF files. 
```
for i in `seq 1 22`; do
    zcat chr${i}/chr${i}.vcf.gz | java -jar gtstats.jar > chr${i}/chr${i}.gtstats
done
```

Next, we filter variants to be excluded from the VCF files on SNPs with MAF $\geq 1\%$. Note that we also need to exclude non-biallelic markers, INDELs, and duplicated positions. This can be done using the following shell script: 
```
#!/bin/bash
i=$1
(awk '{if ($4!="A" && $4!="C" && $4!="T" && $4!="G") print $2}' chr${i}/chr${i}.gtstats
 awk '{if ($5!="A" && $5!="C" && $5!="T" && $5!="G") print $2}' chr${i}/chr${i}.gtstats
 awk '{if ($11<0.01) print $2}' chr${i}/chr${i}.gtstats
 cut -f1,2 chr${i}/chr${i}.gtstats | sort -n | uniq -d | awk '{print $2}') | sort | uniq > chr${i}/chr${i}.exclude.for.snp.maf0.01
```
Suppose the above script is saved under the name `filter.maf.01.sh`: 
```
for i in `seq 1 22`; do
    filter.maf.01.sh $i
done
```

Similarly, for filtering variants to be excluded from the VCF files on SNPs with MAF $< 1\%$: 
```
#!/bin/bash
i=$1
(awk '{if ($4!="A" && $4!="C" && $4!="T" && $4!="G") print $2}' chr${i}/chr${i}.gtstats
 awk '{if ($5!="A" && $5!="C" && $5!="T" && $5!="G") print $2}' chr${i}/chr${i}.gtstats
 awk '{if ($11>=0.01||$11==0) print $2}' chr${i}/chr${i}.gtstats
 cut -f1,2 chr${i}/chr${i}.gtstats | sort -n | uniq -d | awk '{print $2}') | sort | uniq > chr${i}/chr${i}.exclude.for.snp.rare
```
Suppose the above script is saved under the name `filter.rare.sh`: 
```
for i in `seq 1 22`; do
    filter.rare.sh $i 
done
```

We then make VCF files for SNP markers with MAF $\geq 1\%$ or MAF $< 1\%$ by filtering out lines in the VCF files that corresponds to these variants for exclusion. We can use the software [`filterlines.jar`](https://faculty.washington.edu/browning/beagle_utilities/utilities.html#filterlines). 
```
for i in `seq 1 22`; do
    zcat chr${i}/chr${i}.vcf.gz | java -jar filter-lines.jar '#' -2 chr${i}/chr${i}.exclude.for.snp.maf0.01 | gzip > chr${i}/chr${i}_snp_maf0.01.vcf.gz

    zcat chr${i}/chr${i}.vcf.gz | java -jar filter-lines.jar '#' -2 chr${i}/chr${i}.exclude.for.snp.rare | gzip > chr${i}/chr${i}_snp_rare.vcf.gz
done
```

For the next step, we need to recode the 0/1 field in the VCF files so that 0/1 represents major/minor allele instead of REF/ALT. 

We need the genotype statistics on the new VCF files for SNPs with MAF $\geq 1\%$ or MAF $< 1\%$ to know the major/minor allele at each marker: 
```
for i in `seq 1 22`; do
    zcat chr${i}/chr${i}_snp_maf0.01.vcf.gz | java -jar gtstats.jar > chr${i}/chr${i}_maf0.01snp.gtstats
done

for i in `seq 1 22`; do
    zcat chr${i}/chr${i}_snp_rare.vcf.gz | java -jar gtstats.jar > chr${i}/chr${i}_raresnp.gtstats
done
```

Next, we identify markers for which we need to recode the 0/1 fields (i.e., the non-REF allele is not the minor allele): 
```
for i in `seq 1 22`; do
    cat chr${i}/chr${i}_maf0.01snp.gtstats | awk -v OFS='\t' '{print \$0,\"nochange\"}' | awk -v OFS='\t' '{if (\$11!=\$13) \$16=\"change\";print \$0}' > chr${i}/chr${i}.maf0.01snp.markers2recode

    cat chr${i}/chr${i}_raresnp.gtstats | awk -v OFS='\t' '{print \$0,\"nochange\"}' | awk -v OFS='\t' '{if (\$11!=\$13) \$16=\"change\";print \$0}' > chr${i}/chr${i}.raresnp.markers2recode"
done
```

And we change the format of the VCF files, so that each line now all has 0's and 1's separated by tabs representing haplotypes:
```
for i in `seq 1 22`; do
    zcat chr${i}/chr${i}_snp_maf0.01.vcf.gz | grep -v '^#' | awk -v OFS='\t' '{for(i=10; i<=NF; i++) \$i=substr(\$i,1,3); print \$0}' | sed 's|\||\t|g' > chr${i}/chr${i}.snp.maf0.01.gt

    zcat chr${i}/chr${i}_snp_rare.vcf.gz | grep -v '^#' | awk -v OFS='\t' '{for(i=10; i<=NF; i++) \$i=substr(\$i,1,3); print \$0}' | sed 's|\||\t|g' > chr${i}/chr${i}.snp.rare.gt
done
```

We provide a python program `changevcf.py` for recoding each 0/1 field to represent major/minor alleles. `Nhaplotype` is the number of haplotypes in the sample, which is twice of the sample size. 

```
for i in `seq 1 22`; do
    python changevcf.py chr${i}/chr${i}.maf0.01snp.markers2recode chr${i}/chr${i}.snp.maf0.01.gt ${Nhaplotype} | gzip > chr${i}/chr${i}.maf0.01snp.gt.recode.gz
done
for i in `seq 1 22`; do
    rm chr${i}/chr${i}.snp.maf0.01.gt
done

for i in `seq 1 22`; do
    python changevcf.py chr${i}/chr${i}.raresnp.markers2recode chr${i}/chr${i}.snp.rare.gt ${Nhaplotype} | gzip > chr${i}/chr${i}.raresnp.gt.recode.gz
done
for i in `seq 1 22`; do
    rm chr${i}/chr${i}.snp.rare.gt 
done
```
We need to add a header line to the output file of haplotypes. The 1-9 entries in the header line will be the 1-9 column headers of the original VCF file. Each subject ID in the original VCF file will be changed to ID of the two haplotypes of that subject, i.e., ID of the haplotype of the first sample will be HAP1 \t HAP2. Say we save this header line in a text file named `HAPids`: 

```
for i in `seq 1 22`; do
    (cat HAPids
    zcat chr${i}/chr${i}.maf0.01snp.gt.recode.gz) | gzip > chr${i}/chr${i}.maf0.01snp.gt.recode.tmp.gz

    (cat HAPids
    zcat chr${i}/chr${i}.raresnp.gt.recode.gz) | gzip > chr${i}/chr${i}.raresnp.gt.recode.tmp.gz
done

for i in `seq 1 22`; do
    mv chr${i}/chr${i}.maf0.01snp.gt.recode.tmp.gz chr${i}/chr${i}.maf0.01snp.gt.recode.gz
    mv chr${i}/chr${i}.raresnp.gt.recode.tmp.gz chr${i}/chr${i}.raresnp.gt.recode.gz
done
```

One last step before getting mutation counts is to extract the genotype statistics required for getting mutation counts from the gtstats output on the SNP VCF files: 
```
for i in `seq 1 22`; do
    cat chr${i}/chr${i}_maf0.01snp.gtstats | awk -v OFS='\t' '{print \$2\"_\"\$4\"_\"\$5,\$0}' | cut -f2,3,5,6,14,15 > chr${i}/chr${i}.maf0.01snp.gtstats

    cat chr${i}/chr${i}_raresnp.gtstats | awk -v OFS='\t' '{print \$2\"_\"\$4\"_\"\$5,\$0}' | cut -f2,3,5,6,14,15 > chr${i}/chr${i}.raresnp.gtstats
done
```

We provide python programs `sql_gtcount.py` and `sql_gtcount_split.py` to get mutation counts from files that record haplotypes following the above format. The later program split the input data into 10 or 11 separate files for parallel computing.  

For computation efficiency, we can first load the required data to database using the provided helper python program `sql_load.py`, for example, using the following shell script: 

```
#!/bin/bash
i=$1
Nhaplotype=$2
export PATH=/python/anaconda2/bin:$PATH 
echo "python sql_load.py chr${i}/chr${i}.maf0.01snp.gt.recode.gz chr${i}/chr${i}.maf0.01snp.gtstats chr${i}/chr${i}.maf0.01snp.gt.recode.db $Nhaplotype"
echo "python sql_load.py chr${i}/chr${i}.raresnp.gt.recode.gz chr${i}/chr${i}.raresnp.gtstats chr${i}/chr${i}.raresnp.gt.recode.db $Nhaplotype"
```
And then to execute the shell script (saved as `loadpy.sh`): 
```
for i in `seq 1 22`; do
    loadpy.sh $i $Nhaplotype 
done
```

Below is an example shell script to get mutation counts from input data on three-way IBD sharing using the parallel version of the program: 
```
#!/bin/bash

i=$1
j=$2
ibdfile=$3
mapfile=$4
idfile=$5
usephase=$6

echo "chr=$i"
echo "split=$j"

echo "copy database"
mkdir /scratch/chr${i}_split${j}
cp chr${i}/chr${i}.maf0.01snp.gt.recode.db /scratch/chr${i}_split${j}
cp chr${i}/chr${i}.raresnp.gt.recode.db /scratch/chr${i}_split${j}
export PATH=/python/anaconda2/bin:$PATH

echo "get mutation count"
python sql_gtcount_split.py ibdtrios/chr${i}.maf0.25.2.5-6cM.txt /scratch/chr${i}_split${j}/chr${i}.raresnp.gt.recode.db /scratch/chr${i}_split${j}/chr${i}.maf0.01snp.gt.recode.db ${ibdfile} ${mapfile} ${idfile} $j $usephase> chr${i}/triosgt.maf0.25.2.5-6cM.${j}.txt

echo "clear repository in the end"
rm -r /scratch/chr${i}_split${j}  
echo ""
```
When executing the above shell script, the required command-line inputs include

- `chr`: chromosome number 1-22
- `j`: part number of the separated data
- `ibdfile`: path to the IBD file
- `mapfile`: path to the genetic maps
- `idfile`: path to a text file that has the ID names of all samples as a single column
- `usephase`: if usephase=true, the phase will be used as given without adjustment

To execute the shell script (named `gtcountpy_split.sh`)
```
for chr in `seq 1 22`; do
    for j in `seq 1 10`; do
        gtcountpy_split.sh $chr $j ${ibdfile} ${mapfile} ${idfile} false
    done
done
```

Finally, we collect the output from each parallel process together: 
```
for i in `seq 1 22`; do
    cat chr${i}/triosgt.maf0.25.2.5-6cM.[1-9]*.txt > chr${i}/triosgt.maf0.25.2.5-6cM.txt
done

for i in `seq 1 22`; do
    rm chr${i}/triosgt.maf0.25.2.5-6cM.[1-9]*.txt
done
```

## STEP 5. Estimate mutation rate
Before finally estimating the mutation rate, we need two more steps. First, we remove mutation count greater and equal to 50 using the R program `checkgt.R`. 
```
Rscript checkgt.R
```
Second, we need to estimate heterozygosity rates using the R program `het.R` based on the Hardy-Weinberg exact test p-value reports generated by [`plink`](https://www.cog-genomics.org/plink/) on all 22 chromosomes: 
```
#!/bin/bash
for i in {1..22}; do
    plink --vcf chr${i}/chr${i}.vcf.gz --hardy
    mv plink.hwe chr${i}/chr${i}.plink.hwe 
done

Rscript het.R
```

Now we are ready to conduct a two-stage likelihood calculation to estimate the mutation rate. We use the software `addGC.v3.jar`. In the first stage, we find MLEs for the two error rates $\epsilon_1$ and $\epsilon_2$ seperately. In the code below, `mapfile` denotes path to the genetic map file, `Nefile` denotes path to the output .ne file by IBDNe, and `hetrate1` and `hetrate2` denote the estimated heterozygosity rates using `het.R` above.

```
for i in `seq 1 22`;do
    cat chr${i}/triosgt.maf0.25.2.5-6cM.txt | java -jar ../addGC.v3.jar map=${mapfile} ne=${Nefile} out=chr${i}/mle_epsilon1_maf0.25 mu.start=1.0E-8 mu.end=1.7E-8 mu.step=1.0E-9 gc.start=5E-7 gc.end=5E-6 gc.step=5.0E-8 theta1=${hetrate1} theta2={hetrate2} stage=1 err1.start=1.0E-13 err1.ratio=10
done

for i in `seq 1 22`;do
    cat chr${i}/triosgt.maf0.25.2.5-6cM.txt | java -jar ../addGC.v3.jar map=${mapfile} ne=${Nefile} out=chr${i}/mle_epsilon2_maf0.25 mu.start=1.0E-8 mu.end=1.7E-8 mu.step=1.0E-10 gc.start=5E-7 gc.end=5E-6 gc.step=5.0E-8 theta1=${hetrate1} theta2={hetrate2} stage=2 err2.start=1.0E-13 err2.ratio=10
done
```

The software `addGC.v3.jar` will produce a output file that contains the log-likelihood calculated at each parameter value. We provide R programs `search.gc.eps1.R` and `search.gc.eps2.R` to find the value of $\epsilon_1$ and $\epsilon_2$ that maximize the likelihood. 
```
Rscript search.gc.eps1.R ${path_to_likelihood_file}
Rscript search.gc.eps2.R ${path_to_likelihood_file}
```

Finally, we plug in the MLEs of $\epsilon_1$ (`e1_mle` below) and $\epsilon_2$ (`e2_mle` below) to calculate the mutation rate $\mu$, with the R program `search.gc.mu.R`. 
```
for i in `seq 1 22`;do
    echo "cat chr${i}/triosgt.maf0.25.2.5-6cM.txt | java -jar ../addGC.v3.jar map=${mapfile} ne=${Nefile} out=chr${i}/mle_mu_maf0.25 mu.start=1.0E-8 mu.end=1.7E-8 mu.step=1.0E-10 gc.start=5.0E-7 gc.end=5.0E-6 gc.step=5.0E-8 theta1=${hetrate1} theta2={hetrate2} err1.start=${e1_mle} err1.end=${e1_mle} err2.start=${e2_mle} err2.end=${e2_mle} stage=3
done

Rscript search.gc.mu.R ${path_to_likelihood_file}
``` 
