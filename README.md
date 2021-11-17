# Spurp_SexDetermination_RAD
Bioinformatic WorkFlow used to find sex specific genotypes within Spurp using RAD seq
10 M, 18 F

## Process Reads

Reads were processed using fastp and found to already have had been trimmed for barcodes.

## Alignment

Spur_5.0 scaffold assembly: GCA_000002235.4
Lvar_3.0 chromosomal assembly: GCA_018143015.1


**Spur Alignment Statistics** are stored **here**: /02-ALIGN/bam_stats.txt 

## Variant Calling

error message from /usr/local/bin/freebayes : 
/usr/local/bin/freebayes: error while loading shared libraries: libbz2.so.1.0: cannot open shared object file: No such file or directory

Switched to locally installing via anaconda ``conda install -c bioconda freebayes/1.3.5``

Succesfully genotyped ~30,000 sites

## Variant Filtering

differences in sequencing depth makes thresholds hard to determine, min/maxDP filters produces high levels of missingness in Female samples.

<p align="center">
<img src="03-VARIANT/PLOTS/Individual_Depth-1.png" width = "45%">
<img src="03-VARIANT/PLOTS/Missing_Individuals-1.png" width = "45%">
</p>

**Figure 1.** Distribution of data for individual missingness and mean-depth from Spurp.minQ20.minGQ20.mac3.miss1.snps.ri.vcf (~10,000 SNPs)

<p align="center">
<img src="04-PCA/PLOTS/PCA-1-1.png" width = "55%">
</p>

**Figure 2.** principal component analysis via plink ~6% variance explained on both axis (PC1 vs. PC3)

## Presence & Absence w/ Sex



## Bedtools 
``bedtools coverage``

## Mapping

map short reads to a genome with assigned chromosomes.
