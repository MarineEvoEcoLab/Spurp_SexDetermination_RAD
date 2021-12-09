Landscape Genomics PCA Analysis
================
Gabriel Barrett
9/17/2020

## PCA Analysis From Plink

``` r
#load eigenvectors/values 
pca <- read.delim("Spurp.minQ20.minGQ20.mac4.miss99.eigenvec", header = FALSE, sep = " ")
eigenval <- scan("Spurp.minQ20.minGQ20.mac4.miss99.eigenval")
#remove the first column since it's repeated twice
pca <- pca[,-1]
#name the first column indv
names(pca)[1] <- "ind"
#name 2nd to number of columns in pca as PC 1 to 20
names(pca)[2:ncol(pca)] <- paste0("PC",1:(ncol(pca)-1))
meta <- read.delim("../meta/Sp_Radseq_Files_Guide.txt") %>% select(Urchin, Sex)

df <- left_join(pca,meta,by=c("ind"="Urchin"))
```

``` r
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_classic()
```

![](D:/Puritz_Lab/Spurp_SexDetermination_RAD/04-PCA/PLOTS/Variance%20Explained-1.png)<!-- -->

``` r
library(ggsci)
library(viridis)
```

    ## Warning: package 'viridis' was built under R version 4.0.5

    ## Loading required package: viridisLite

    ## Warning: package 'viridisLite' was built under R version 4.0.5

``` r
ggplot(df, aes(PC2, PC4, color=Sex)) + geom_point(size = 3, alpha = .875) + 
  coord_equal() + 
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  labs(colour = "") + #stat_ellipse(aes(group = reg)) +
  #annotate("text", x = -0.05, y = .0045, label = paste0("PC1 (", signif(pve$pve[1], 3), "%)"),size=4.95) +
  #annotate("text", x = -.0045, y = -.045, label = paste0("PC2 (", signif(pve$pve[2], 3), "%)"), angle = 90,size=4.95) +
  #scale_color_igv(breaks = c("up","down","lower")) + 
  geom_hline(yintercept=0, linetype="longdash") + geom_vline(xintercept=0, linetype="dotdash") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    legend.position = c(.7,.9)
  ) + 
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))
```

![](D:/Puritz_Lab/Spurp_SexDetermination_RAD/04-PCA/PLOTS/PCA-1-1.png)<!-- -->

``` r
#ggsave("LG.snps1.FFP5.MAF.01.HWE.01.MISS.65.MD.69.PCA.png",last_plot(),device="png")
```
