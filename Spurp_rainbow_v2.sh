#!/bin/local/env bash

# ~~~~~~~~~~~~~~~
# DE-NOVO CONSTRUCTION ON PROCESSED READS
# ~~~~~~~~~~~~~~~

cd /home/Shared_Data/Spurp_RAD/test

# Cluster Reads
bigfile1=../01-PROCESS/bigfile.1.fq.gz
bigfile2=../01-PROCESS/bigfile.2.fq.gz
rainbow cluster -1 $bigfile1 -2 $bigfile2 -m 3 -e 2000 > rbcluster.out 2> rbcluster.log

rainbow div -i rbcluster.out -o rbdiv.out -k 7 -K 10 -f .5 2> rbdiv.log

rainbow merge -i rbdiv.out -o rbmerge.out -a -r 5 -N300 -R300 -l 20 -f .9 2> rbmerge.log

# Consensus De-Novo Reference
/usr/local/bin/select_best_rbcontig_plus_read1.pl rbmerge.out rbdiv.out > rainbow.fasta

# ~~~~~~~~~~~~~~~~~~
# ALIGNING INDV TO DE-NOVO REFERENCE
# ~~~~~~~~~~~~~~~~~~

# Index De-Novo Reference Genome
bwa index -a bwtsw rainbow.fasta -p rainbow

FQ1=../01-PROCESS/*.out_1.fq.gz
for fq1 in $FQ1;
do
	PREFIX=$(basename "$fq1" .out_1.fq.gz)
	fq2=../01-PROCESS/"${PREFIX}.out_2.fq.gz"
	BAM="${PREFIX}.bam"
	RG=$(echo \@RG\\tID:$PREFIX\\tSM:$PREFIX)

	# Align Reads
	bwa mem -t 15 -R $RG rainbow $fq1 $fq2 | \
	samtools view -S -h -q 30 - | \
	samtools sort - > $BAM
	samtools index $BAM

	samtools coverage $BAM > $PREFIX.cov
done


#  
# Calculate Average Coverage
#

OUT=FINAL_RESULTS
mkdir ./${OUT}

echo "read coverage meandepth meanbaseq meanmapq" > $OUT/avg_female.cov
	# skip first line of female bam files

	# stache coverage from each bam file and calculate average for each row
awk '{
	c[FNR]=$1;a[FNR]+=$6;d[FNR]+=$7;e[FNR]+=$8;f[FNR]+=$9;b[FNR]++;
}END{
	# skip the first line by setting i = 2
	for(i=2;i<=FNR;i++)
	print c[i],a[i]/b[i],d[i]/b[i],e[i]/b[i],f[i]/b[i]
}' *f*.cov >> $OUT/avg_female.cov

	# repeat for males
echo "read coverage meandepth meanbaseq meanmapq" > $OUT/avg_male.cov
awk '{
        c[FNR]=$1;a[FNR]+=$6;d[FNR]+=$7;e[FNR]+=$8;f[FNR]+=$9;b[FNR]++;
}END{
        for(i=2;i<=FNR;i++)
        print c[i],a[i]/b[i],d[i]/b[i],e[i]/b[i],f[i]/b[i]
}' *m*.cov >> $OUT/avg_male.cov

# Write reads with 0 coverage to FINAL_RESULTS respective of sex

awk '{if ($3 == 0) print}' $OUT/avg_male.cov > $OUT/avg_male_0.cov

awk '{if ($3 == 0) print}' $OUT/avg_female.cov > $OUT/avg_female_0.cov


#
# Retain only unique calls between sexes 
#

touch full_0.cov
cat avg_male_0.cov avg_female_0.cov > full_0.cov
# Gather duplicates into file
awk 'seen[$1]++' full_0.cov > duplicate_readID.txt
# Remove READS present in file1 from file2
awk -F '\t' 'FNR==NR{a[$1];next};!($1 in a)' duplicate_readID.txt avg_male_0.cov > avg_male_0_uniq.cov

awk -F '\t' 'FNR==NR{a[$1];next};!($1 in a)' duplicate_readID.txt avg_female_0.cov > avg_female_0_uniq.cov


