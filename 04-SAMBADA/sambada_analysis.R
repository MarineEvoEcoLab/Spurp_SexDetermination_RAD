library(tidyverse)
library(dplyr)
library(vcfR)
# Read result file for univariate models
## Insert correct name here
res.1=read.table("Spurp-Out-1.geno", T, " ", stringsAsFactors=F)

# Read the result file for the null models
## Insert correct name here
res.0=read.table("Spurp-Out-0.geno", T, " ", stringsAsFactors=F)

# Drop the last column in both files since they are null
res.1=res.1[,1:15]
res.0=res.0[,1:5]

# Remove unconverted models (i.e. models with an error)
res.1=res.1[res.1$NumError == 0,]

# Check number of monomorphic markers
table(res.0$NumError)

# Keep the number of polymorphic markers as reference
numMark=table(res.0$NumError)[[1]]

# Number of environmental variables
## Insert correct number here
numEnv=1

# Choice of alpha (could be 0.05)
alpha=0.01

# Score threshold for model significance (same value for G and Wald)
# Models with a score higher than this value are considered as significant.
# When using the option "SIGNIF" or "BEST", Sambada enforce this threshold for both G and Wald scores. One can also use a single test (G or Wald) for assessing significance.
# Note: The number of degree of freedom (df) is 1 since we compare models involving one env var with models without any env var
score.threshold=qchisq(alpha/(numEnv*numMark), 1, lower.tail=F)

# Score threshold for model significance (same value for G and Wald)
# Models with a p-value lower than this value are considered as significant.
# When using the option "SIGNIF" or "BEST", Sambada enforce this threshold for both G and Wald p-values. One can also use a single test (G or Wald) for assessing significance.
pVal.threshold=alpha/(numEnv*numMark)
print(pVal.threshold)
# Computing p-values:
res.1=cbind(res.1, pvalG=pchisq(res.1$Gscore, 1, lower.tail=F), pvalWald=pchisq(res.1$WaldScore, 1, lower.tail=F))
#write.table(res.1,file="res.1.txt")
# Selecting models passing the G test (p-value for G score lower than the threshold)
res.1[res.1$pvalG<pVal.threshold,]
# Check if there are matches here

# If not, check if there are a couple of significant models if setting alpha to 0.05
alpha=0.05
pVal.threshold=alpha/(numEnv*numMark)
res.1[res.1$pvalG<pVal.threshold,]

res.1.a <- res.1 %>% select(Marker,Gscore,WaldScore,pvalG,pvalWald) %>% arrange(pvalG)
write.table(x = res.1.a,file = "sambada_results.tsv",quote = F,sep = "\t",row.names = F,col.names = T)


entap <- read.table("../05-ENTAP/final_annotations_lvl3.tsv",sep="\t",fill = T)
vcfR_obj <- read.vcfR("../05-SNPEFF/Spurp.minQ20.minGQ20.mac4.miss95.snps.ri.F.ann.vcf")

vcfR_df <- vcfR2tidy(vcfR_obj)

vcf_df_fix <- vcfR_df$fix %>% 
  select(chrom_id=CHROM,pos=POS,ANN) %>% 
  separate(col = ANN,into = c("allele","annotation","annotation_impact","gene_name",
                              "gene_id","feature_type","feature_id","transcript_biotype",
                              "rank","HGVS.c","HGVS.p","cDNA.pos/cDNA.length","cDNS.pos/cDNS.length",
                              "AA.pos/AA.length","distance","err"), sep = "\\|",convert=T,remove=F,fill="right") %>%
  select(!(ANN)) 
View(vcf_df_fix)






