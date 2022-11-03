#!/bin/local/env bash
PROG=/usr/local/bin
DIR=/home/Shared_Data/Spurp_RAD
cd /home/Shared_Data/Spurp_RAD/test2_g
#rm ./avg_depth.txt
rm *.cov
for i in *bam;
do
	#AVG_DEPTH=$($PROG/samtools depth -Q 40 $i | awk '{sum+=$3} END {print sum/NR}') 
	#SD_DEPTH=$($PROG/samtools depth -Q 40 $i | awk '{sum+=$3;a[NR]=$3} END {for (i in a)y+=(a[i]-(sum/NR))^2;print sqrt(y/(NR-1))}')
	#echo -e "$i\t$AVG_DEPTH\t$SD_DEPTH" >> avg_depth.txt
	echo $i
	COV=$(echo ${i} | sed 's/.bam/.cov/')
	$PROG/samtools coverage ${i} > ${COV}
done

# Calculate read coverage across sex
# write headers
echo "read coverage meandepth meanbaseq meanmapq" > avg_female.cov
	# skip first line of female bam files

	# stache coverage from each bam file and calculate average for each row
awk '{
	c[FNR]=$1;a[FNR]+=$6;d[FNR]+=$7;e[FNR]+=$8;f[FNR]+=$9;b[FNR]++;
}END{
	# skip the first line by setting i = 2
	for(i=2;i<=FNR;i++)
	print c[i],a[i]/b[i],d[i]/b[i],e[i]/b[i],f[i]/b[i]
}' *F*.cov >> avg_female.cov

	# repeat for males
echo "read coverage meandepth meanbaseq meanmapq" > avg_male.cov
awk '{
        c[FNR]=$1;a[FNR]+=$6;d[FNR]+=$7;e[FNR]+=$8;f[FNR]+=$9;b[FNR]++;
}END{
        for(i=2;i<=FNR;i++)
        print c[i],a[i]/b[i],d[i]/b[i],e[i]/b[i],f[i]/b[i]
}' *M*.cov >> avg_male.cov

awk '{if ($3 == 0) print}' avg_male.cov > avg_male_0.cov
awk '{if ($3 == 0) print}' avg_female.cov > avg_female_0.cov

touch full_0.cov
cat avg_male_0.cov avg_female_0.cov > full_0.cov
# Gather duplicates into file
awk 'seen[$1]++' full_0.cov > duplicate_readID.txt
# Remove READS present in file1 from file2
awk -F '\t' 'FNR==NR{a[$1];next};!($1 in a)' duplicate_readID.txt avg_male_0.cov > avg_male_0_uniq.cov

awk -F '\t' 'FNR==NR{a[$1];next};!($1 in a)' duplicate_readID.txt avg_female_0.cov > avg_female_0_uniq.cov


# Low Coverage in the opposite sex
# Female
awk '{if ($6 < 50) print $1}' *F*.cov | awk 'seen[$1]++' > female_low_coverage.txt
# Male
awk '{if ($6 < 50) print $1}' *M*.cov | awk 'seen[$1]++' > male_low_coverage.txt

# Remove Reads from candidate list
# Female
awk -F '\t' 'FNR==NR{a[$1];next};!($1 in a)' male_low_coverage.txt avg_female_0.cov > avg_female_0_uniq_oppositelowcov.cov
# Male
awk -F '\t' 'FNR==NR{a[$1];next};!($1 in a)' female_low_coverage.txt avg_male_0.cov > avg_male_0_uniq_oppositelowcov.cov