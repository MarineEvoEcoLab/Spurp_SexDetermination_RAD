#!/bin/local/env bash

#awk '{print $1}' avg_female_0_uniq.cov > avg_female_0_uniq.txt

#conda activate seqkit
seqkit grep -n -f FemaleCandidateID.txt reference.fasta > female_0.fasta
#conda deactivate 

#conda activate diamond
diamond blastx -p 2 -k 1 -e 0.00001 -d uniprot_sprot_v2.dmnd -q female_0.fasta -o female_0.annot.fa




#conda activate seqkit
seqkit grep -n -f MaleCandidateID.txt reference.fasta > male_0.fasta
#conda deactivate 

#conda activate diamond
diamond blastx -p 2 -k 1 -e 0.00001 -d uniprot_sprot_v2.dmnd -q male_0.fasta -o male_0.annot.fa
#conda deactivate

#awk '{print $1}' avg_female_0_uniq.cov > avg_female_0_uniq.txt

#conda activate seqkit
#seqkit grep -n -f avg_female_0_uniq.txt rainbow_Males.fasta > avg_female_0_uniq.fasta
#conda deactivate 

#conda activate diamond
#diamond blastx -p 2 -k 1 -e 0.00001 -d /home/jgreen/databases/uniprot_sprot.dmnd -q avg_male_0.fasta -o avg_male_0.annot.fa
#conda deactivate

# optimizing the assembly 

# Average coverage and is found in more than 5 individuals minimum coverage in opposite sex

