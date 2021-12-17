# Spurp_SexDetermination_RAD
Bioinformatic WorkFlow used to find sex specific genotypes within Spurp using RAD seq

Purple Sea Urchin: 10 M, 18 F

Extracted ID's from meta data file. 

## Process Reads

Reads were processed using fastp and found to already have had been trimmed for barcodes.

## Reference Alignment

Spur_5.0 scaffold assembly: GCA_000002235.4 
	Mapping Rate ~98%
Lvar_3.0 chromosomal assembly: GCA_018143015.1 
	Mapping Rate ~8%

## Variant Calling

Switched to locally installing via anaconda ``conda install -c bioconda freebayes/1.3.5``

Succesfully genotyped ~30,000 sites

## Variant Filtering

Sites were removed based on genotype quality, mapping quality, and missingness within Sex group (Male, Female).

Differences in sequencing depth makes thresholds hard to determine, min/maxDP filters produces high levels of missingness in Female samples.

<p align="center">
<img src="03-VARIANT/PLOTS/Individual_Depth-1.png" width = "45%">
<img src="03-VARIANT/PLOTS/Missing_Individuals-1.png" width = "45%">
</p>

**Figure 1.** Distribution of data for individual missingness and mean-depth from Spurp.minQ20.minGQ20.mac3.miss1.snps.ri.vcf (~7,000 SNPs)

<p align="center">
<img src="04-PCA/PLOTS/PCA-1-1.png" width = "55%">
</p>

**Figure 2.** principal component analysis via plink ~6% variance explained on both axis (PC1 vs. PC2) (~2103 SNPs)

## Presence/Absence

Performed a logistical regression on individual genotypes using Sambada

## Case/Control

''plink --assoc fisher'' using sex as the encoded phenotype.

## Radsex

<p align="center">
<img src="02-RADSEX/distribution.png" width = "45%">
<img src="02-RADSEX/mapping_circos.png" width = "45%">
</p>

**Figure 3.** heatmap showing association (presence/absence) of marker (read based on minimum depth threshold). R 
circlize package showing significant levels in sex bias on the top track and association on the bottom track. 

> #### Fine Tuning
> - distrib --min-depth <threshold>
> - extract a subset and cluster markers based on depth -> radsex_markers_depth() *subset* command
> - mapping quality threshold  

## De-Novo

Constructed a sex specific reference genome and aligned reads based on opposing sex. 

```bash
#!/bin/local/env bash
PROG=/usr/local/bin
DIR=/home/Shared_Data/Spurp_RAD

rm avg_depth.txt
for i in *bam;
do
        AVG_DEPTH=$($PROG/samtools depth -Q 40 $i | awk '{sum+=$3} END {print sum/NR}')
        SD_DEPTH=$($PROG/samtools depth -Q 40 $i | awk '{sum+=$3;a[NR]=$3} END {for (i in a)y+=(a[i]-(sum/NR))^2;print sqrt(y/(NR-1))}')
        echo -e "$i\t$AVG_DEPTH\t$SD_DEPTH" >> avg_depth.txt
done
```

|ID |avg depth |sd depth |
|---|----------|---------|
|Spf12.bam|1.85434|23.9615|
|Spf13.bam|1.51986|1.47941|
|Spf14.bam|1.43037|2.53296|
|Spf17.bam|5.69743|8.33223|
|Spf18.bam|2.18523|6.11296|
|Spf24.bam|3.23958|5.97413|
|Spf25.bam|2.30345|4.6854|
|Spf26.bam|1.76948|4.92149|
|Spf28.bam|13.9681|23.0103|
|Spf29-3.bam|3.13654|4.26893|
|Spf30.bam|9.09218|19.4664|
|Spf32.bam|15.6657|26.2012|
|Spf35.bam|4.76147|6.8294|
|Spf42.bam|2.51205|55.2533|
|Spf43.bam|2.40871|69.5291|
|Spf44.bam|2.55186|59.204|
|Spf45.bam|2.52544|80.0525|
|Spf6.bam|2.27803|14.7091|
|Spm10.bam|52.8826|142.069|
|Spm11.bam|50.7519|144.781|
|Spm1-3.bam|43.3796|127.697|
|Spm14.bam|44.7545|129.224|
|Spm16.bam|48.8474|133.306|
|Spm2.bam|44.0387|124.847|
|Spm3.bam|58.2594|161.314|
|Spm4.bam|55.4492|154.455|
|Spm7.bam|41.4106|120.693|
|Spm8.bam|50.9566|144.762|


> #### Thoughts/Concerns
> - Inconsistent read depth between sexes makes accurate comparison for sex bias reads difficult. 
> - subset female individuals: Spf28, Spf30, and Spf32 have depths highest average depth
>> - Explore genotypes like reference pipepline
> - resequence female individuals to match male libraries (or include 5 males vs 5 females in one lane with high coverage for confident sex bias calling)

