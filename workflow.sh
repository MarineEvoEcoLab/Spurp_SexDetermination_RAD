#!/usr/bin/env bash
PROG=/usr/local/bin

########################
#### Copy files W/ meta data associated ---------------------------------------------------------
########################
# write unique values from c-4 into txt
awk '!a[$4]++  {print $4}' ./meta/Sp_Radseq_Files_Guide.txt > sample_list.txt

# start from line 2
DATA=$(tail --line=+2 sample_list.txt)

# copy into raw_reads_3
for i in $DATA;
do
        for pair in {1..2};
        do
                cp ./archive/"${i}"/*L3_**"${pair}".fq.gz ./raw_reads/"${i}"_"${pair}".fq.gz

        done
done
rm sample_list.txt

########################
####### Process Raw Reads -------------------------------------------------------------------------
########################
INDIR=/home/Shared_Data/Spurp_RAD/raw_reads_3
OUTDIR=/home/Shared_Data/Spurp_RAD/01-PROCESS
for i in ls $INDIR/*_1.fq.gz;
do
        # Input Reads
        fq1=${i}
        fq2=$(echo ${i} | sed 's/_1./_2./')

        # Output Reads
        FQ1=$(echo ${fq1} | sed -e 's/_1./.out_1./' -e 's+raw_reads_3+01-PROCESS+')
        FQ2=$(echo ${fq2} | sed -e 's/_2./.out_2./' -e 's+raw_reads_3+01-PROCESS+')

        # Output Reports
        OUTPUT=$(basename $fq1 _1.fq.gz)

        /usr/local/bin/fastp -i $fq1 -I $fq2 -o $FQ1 -O $FQ2 --correction -h ${OUTDIR}/$OUTPUT.html &> ${OUTDIR}/$OUTPUT.log

done
##################################
# Reference Genome Pipeline      # -----------------------------------------------------------
# lva chromosomal assembly       # -----------------------------------------------------------
##################################

# lva = ~8% mapping rate
# sp5 = ~95% 

REFDIR=/home/Shared_Data/Spurp_RAD/ref_genome/L_variegatus
REF=${REFDIR}/Lva_genome.FINAL.from_draft.fasta
#INDEX=$(echo $REF | cut -d "/" -f 7 | sed 's/.fasta.gz//')
INDEX=${REFDIR}/Spurp
# Write Index Spurp from reference

bwa index -a bwtsw -p ${INDEX} ${REF}

INDIR=/home/Shared_Data/Spurp_RAD/01-PROCESS

#OUTDIR=/home/Shared_Data/Spurp_RAD/02-ALIGN/sp5
OUTDIR=/home/Shared_Data/Spurp_RAD/02-ALIGN/lva
mkdir -p $OUTDIR

rm $OUTDIR/bam_stats.txt


for i in ${INDIR}/*1.fq.gz;
do
        # RAD Paired Files
        FQ1=${i}
        FQ2=$(echo $FQ1 | sed 's/1.fq.gz/2.fq.gz/')
        # Prefix (ID)
        SAM=$(basename $FQ1 .out_1.fq.gz)
        BAM=$OUTDIR/${SAM}.bam
        #echo $FQ1 $FQ2 $SAM

        RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)

        # Align Reads
        $PROG/bwa mem -t 20 -R $RG ${INDEX} $FQ1 $FQ2 | \
        $PROG/samtools view -S -h -u - | \
        $PROG/samtools sort - > $BAM
        $PROG/samtools index $BAM

        # Get Alignment Statistics
        MM=$($PROG/samtools flagstat $BAM | grep -E 'mapped \(|properly' | cut -f1 -d '+' | tr -d '\n')
        CM=$($PROG/samtools idxstats $BAM | mawk '$3 > 0' | wc -l)
        CC=$($PROG/samtools idxstats $BAM | mawk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
        DD=$($PROG/samtools idxstats $BAM | mawk '$3 >0' | mawk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
        BM=$($PROG/samtools flagstat $BAM | grep mapQ | cut -f1 -d ' ')
        MR=$($PROG/samtools flagstat $BAM |  grep "mapped (" | awk '{print $5}' | cut -b 2-5)
        echo -e "$SAM\t$MM\t$CM\t$CC\t$DD\t$BM\t$MR" >> $OUTDIR/bam_stats.txt

done


######################
###### Variant Calling --------------------------------------------------------------------------------
# Freebayes

######### input, output, directories 
#INDIR=/home/Shared_Data/Spurp_RAD/02-ALIGN/L_variegatus
#OUTDIR=/home/Shared_Data/Spurp_RAD/03-VARIANT
#mkdir -p $OUTDIR

# reference genome
sp5_REF=/home/Shared_Data/Spurp_RAD/ref_genome/sp5_0_GCF_genomic.fa
lva_REF=/home/Shared_Data/Spurp_RAD/ref_genome/L_variegatus/Lva_genome.FINAL.from_draft.fasta
# list of bam files
#SAMPLES_FILE=/home/Shared_Data/Spurp_RAD/sample_list.txt
#BAMLIST=$OUTDIR/bam.list
#tail -n +2 $SAMPLES_FILE | sed 's/$/.bam/' | sed "s,^,$INDIR/," > $BAMLIST


        # Directories
        INDIR=/home/Shared_Data/Spurp_RAD/02-ALIGN/$i
        OUTDIR=/home/Shared_Data/Spurp_RAD/03-VARIANT/$i
        mkdir -p $OUTDIR

        # list of bam files
        SAMPLES_FILE=/home/Shared_Data/Spurp_RAD/sample_list.txt
        BAMLIST=$OUTDIR/bam.list
        tail -n +2 $SAMPLES_FILE | sed 's/$/.bam/' | sed "s,^,$INDIR/," > $BAMLIST

        ######### run freebayes -------------------------------------
        freebayes -f $lva_REF --bam-list $BAMLIST \
                -m 30 -q 20 \
                --min-coverage 100 --skip-coverage 50000 > ${OUTDIR}/Spurp_.vcf


########################
###### Variant Filtering -----------------------------------------------------------------------------
########################

# Setup FILES
FILE=${DIR}/Spurp.vcf
CLEANED_FILE=${DIR}/Spurp.F.vcf
F1=$(echo $FILE | sed 's/vcf/minQ20.minGQ20.mac4.F.vcf/')
F2=$(echo $F1 | sed 's/vcf/miss95.F.vcf/')
F3=$(echo $F2 | sed 's/vcf/snps.F.vcf/')
F4=$(echo $F3 | sed 's/vcf/ri.F.vcf/')

# Remove Problematic Indv.
$PROG/bcftools --vcf $FILE --samples-file ${DIR}/problematic.indv.txt --out

# Filter Sites
$PROG/vcftools --vcf $CLEANED_FILE --minQ 20 --minGQ 20 --mac 4 --recode --recode-INFO-all --stdout > $F1

########## per-pop miss --------------------------------------------

# write popmap.txt
awk '{print $4"\t"$3}' /home/Shared_Data/Spurp_RAD/meta/Sp_Radseq_Files_Guide.txt | tail -n +2 | sort | uniq > popmap.txt

sh filter_miss_by_pop.sh $FILTER1 popmap.txt 1 .95 $FILTER2


#####################
###### Principal component analysis ------------------------------------------------------------------
#####################
DIR=/home/Shared_Data/Spurp_RAD

# v1.07 ------
plink=$PROG/plink
# IN ---------
VCF=$DIR/03-VARIANT/Spurp.minQ20.minGQ20.mac4.miss99.F.ann.vcf
# OUT --------
OUTDIR=$DIR/04-PCA
mkdir $OUTDIR
OUTFILE=$OUTDIR/Spurp.minQ20.minGQ20.mac4.miss99.F

########## plink PCA -------------------------------------------------------

$plink --vcf $VCF --indep-pairwise 100 25 0.1 \
--double-id --allow-extra-chr --set-missing-var-ids @:# --out $OUTFILE

$plink --file $OUTFILE --extract ${OUTFILE}.prune.in \
--double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out $OUTFILE

####################
###### Plink case/control association ------------------------------------------------------------
####################
SDIR=${DIR}/04-PLINK
mkdir $SDIR


# input files
PED=$DIR/04-PCA/Spurp.minQ20.minGQ20.mac4.miss99
OUT=$SDIR/Spurp.minQ20.minGQ20.mac4.miss99
# write covariate file
cat ${DIR}/meta/Sp_Radseq_Files_Guide.txt | awk -F "\t" '{print $3"\t"$3"\t"$2}' | uniq | \
sed -e 's/Female/1/' -e 's/Male/2/' -e '1,/Urchin/ s/Urchin/FID/' -e 's/Urchin/IID/' > $SDIR/pheno.dat

#plink --file $PED --make-pheno $pheno.dat --allow-extra-chr

# logistic
plink --file $PED --pfilter .01 --logistic \
--pheno $SDIR/pheno.dat --pheno-name Sex --out $OUT --noweb --allow-extra-chr --allow-no-sex

# filter sig. and output ID
#awk '{if( $9 < .01) {print $2}}' | cut -d : $OUT.assoc > sig_snps.txt
awk '{if( $9 < .01) {print $2}}' $OUT.assoc.logistic | sed 's/:/\t/' > sig_snps.txt

# sig SNP df w/ annotations
bcftools query -T sig_snps.txt --format '%CHROM\t%POS\t%ANN\t[%SAMPLE=%GT ]\n' ../03-VARIANT/Spurp.minQ20.minGQ20.mac4.miss95.snps.ri.F.ann.vcf > sig_snps_df.txt
#NUMCOL=$(head -1 fisher_snps.txt | grep -o "|" | wc -l)


#######################
# De-Nove Assembly    # -----------------------------------------------------------------------
# Rainbow 2.0.4       # -----------------------------------------------------------------------
#######################

# Assemble based on Sex (M,F)

DIR=/home/Shared_Data/Spurp_RAD
OUT=/home/Shared_Data/Spurp_RAD/02-RAINBOW
META=/home/Shared_Data/Spurp_RAD/meta/Sp_Radseq_Files_Guide.txt
mkdir $OUT

# two files with individuals grouped by sex
awk -F"\t" '{if ($2 == "Male") print $3}' $META | uniq | sed "s,^,$DIR/01-PROCESS/," > Males.txt
awk -F"\t" '{if ($2 == "Female") print $3}' $META | uniq | sed "s,^,$DIR/01-PROCESS/," > Females.txt


for sex in Males Females;
do
        cat ${sex}.txt | while read ind;
        do
                FQ1=${ind}.out_1.fq.gz
                FQ2=$(echo $FQ1 | sed -e 's/1.fq.gz/2.fq.gz/')
                echo "Read 1: $FQ1 Read 2: $FQ2"
                rainbow cluster -1 $FQ1 -2 $FQ2 -m 4 -e 2000 > $OUT/rbcluster_$sex.out 2> $OUT/rbcluster_$sex.log
        done

        rainbow div -i $OUT/rbcluster_$sex.out -o $OUT/rbdiv_$sex.out -k 7 -K 10 -f .5 2> $OUT/rbdiv_$sex.log

        rainbow merge -i $OUT/rbdiv_$sex.out -o $OUT/rbmerge_$sex.out -a -r 2 -N1000 -R1000 -l 20 -f .9 2> $OUT/rbmerge_$sex.log
        /usr/local/bin/select_best_rbcontig_plus_read1.pl
        PROG=/usr/local/bin
        
	# write best contigs into fasta format
	$PROG/select_best_rbcontig_plus_read1.pl $OUT/rbmerge_$sex.out $OUT/rbdiv_$sex.out > $OUT/rainbow_$sex.fasta

done


############################################
# Align reads to reference of opposite sex # ------------------------------------------------
############################################

for sex in Males Females;
do

done

