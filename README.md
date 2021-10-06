# Spurp_SexDetermination_RAD
Bioinformatic WorkFlow of sex determination within Spurp using RAD seq

## Process Reads

```shell
#!/usr/bin/env bash

INDIR=/home/Shared_Data/Spurp_RAD/raw_reads_3
OUTDIR=/home/Shared_Data/Spurp_RAD/01-PROCESS
for i in ls $INDIR/*_1.fq.gz;
do
        fq1=${i}
        fq2=$(echo ${i} | sed 's/_1./_2./')

        FQ1=$(echo ${fq1} | sed -e 's/_1./.out_1./' -e 's+raw_reads_3+01-PROCESS+')
        FQ2=$(echo ${fq2} | sed -e 's/_2./.out_2./' -e 's+raw_reads_3+01-PROCESS+')

        OUTPUT=$(basename $fq1 _1.fq.gz)

        /usr/local/bin/fastp -i $fq1 -I $fq2 -o $FQ1 -O $FQ2 --correction -h ${OUTDIR}/$OUTPUT.html &> ${OUTDIR}/$OUTPUT.log


done
```
