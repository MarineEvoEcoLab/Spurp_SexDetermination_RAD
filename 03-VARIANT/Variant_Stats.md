Variant Statisticsr
================
Gabriel Barrett
10/13/2021

``` r
knitr::opts_knit$set(root.dir = "D:/Puritz_Lab/Spurp_SexDetermination_RAD/03-VARIANT/")

knitr::opts_chunk$set(fig.path = "PLOTS/")

library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.5     v purrr   0.3.4
    ## v tibble  3.1.3     v dplyr   1.0.7
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.4.0     v forcats 0.5.0

    ## Warning: package 'ggplot2' was built under R version 4.0.5

    ## Warning: package 'dplyr' was built under R version 4.0.5

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(ggpubr)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(readxl)
```

``` r
var_qual <- read_delim("STATS/fb.lqual", delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)

q <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_classic()
q 
q + coord_cartesian(xlim = c(-0.01,50000))
```

``` r
var_depth <- read_delim("STATS/fb.ldepth", delim = "\t",
                        col_names = c("chr", "pos", "sum_depth", "sumsq_depth"), skip = 1)

sd <- ggplot(var_depth, aes(sum_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_classic() + xlab("Sum Site Depth") + ylab("Density") + theme(plot.margin=grid::unit(c(0,0,0,0), "cm")) + coord_cartesian(expand=FALSE)

sd100 <- sd + xlab("Sum Site Depth (0-100)") + coord_cartesian(xlim = c(-0.01,100), ylim = c(0,.25), expand=F)
ggarrange(sd + labs(tag="A"), sd100 + labs(tag="B"), ncol=2)

summary(var_depth$sum_depth)
```

``` r
var_mean_depth <- read_delim("STATS/Spurp.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   chr = col_character(),
    ##   pos = col_double(),
    ##   mean_depth = col_double(),
    ##   var_depth = col_double()
    ## )

``` r
msd <- ggplot(var_mean_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_classic() + xlab("Mean Site Depth") + ylab("Density") + theme(plot.margin=grid::unit(c(0,0,0,0), "cm")) + coord_cartesian(expand=FALSE)

msd100 <- msd + xlab("Mean Site Depth (0-150)") + coord_cartesian(xlim = c(-0.01,150), ylim = c(0,.25), expand=F)
```

    ## Coordinate system already present. Adding new coordinate system, which will replace the existing one.

``` r
summary(var_mean_depth$mean_depth)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    1.50   46.33   52.26   59.39   62.81 2235.67

``` r
mean <- mean(var_mean_depth$mean_depth)
std <- sd(var_mean_depth$mean_depth)
cutoff <- sum(mean + (2*std))

ggarrange(msd + labs(tag="A"), msd100+ labs(tag="B"), ncol=2)
```

![](PLOTS/MEAN%20SITE%20DEPTH-1.png)<!-- -->

``` r
var_miss <- read_delim("STATS/Spurp.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   chr = col_character(),
    ##   pos = col_double(),
    ##   nchr = col_double(),
    ##   nfiltered = col_double(),
    ##   nmiss = col_double(),
    ##   fmiss = col_double()
    ## )

``` r
summary(var_miss$fmiss)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 0.00000 0.05660 0.07692 0.08949 0.12000 0.80645

``` r
ms <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue2", colour = "black", alpha = 0.3) + 
  theme_classic() + xlab("Fraction Missing Per Site") + ylab("Density") +  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"))  + coord_cartesian(xlim = c(-0.01,1),expand=F)
ms
```

![](PLOTS/Missing%20Sites-1.png)<!-- -->

``` r
var_freq <- read_delim("STATS/Spurp.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# calculate the minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

freq_sum<-var_freq%>%
  summarise(avg=mean(maf), med=median(maf))

# Distribution of allele frequencies
af <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue3", colour = "black", alpha = 0.3) + 
  theme_classic() + xlab("Minor Allele Frequency (maf)") + ylab("Density") +
  geom_text(x = 0.3, y = 40,aes(label=paste0("Mean:",round(avg,3)," Median:",round(med,3))),data = freq_sum) +  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"))  + coord_cartesian(xlim = c(-0.009,.51),expand=F)
af
```

``` r
var_hwe <- read_delim("STATS/Spurp.hwe", delim = "\t",
                      col_names = c("chr", "pos", "obs", "exp", "chi_sq", "p_hwe", "p_het_def", "p_het_exc"), skip=1)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   chr = col_character(),
    ##   pos = col_double(),
    ##   obs = col_character(),
    ##   exp = col_character(),
    ##   chi_sq = col_double(),
    ##   p_hwe = col_double(),
    ##   p_het_def = col_double(),
    ##   p_het_exc = col_double()
    ## )

``` r
summary(var_hwe$p_hwe)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   0.000   1.000   1.000   0.856   1.000   1.000

``` r
hwe <- ggplot(var_hwe, aes(p_hwe)) + geom_density(fill = "dodgerblue3", colour = "black", alpha = 0.3) + theme_classic()
hwe
```

![](PLOTS/Hardy-Weinberg%20Equilibrium-1.png)<!-- --> \#\# Individual
Stats

``` r
ind_miss  <- read_delim("STATS/Spurp.minQ20.minGQ20.mac3.miss1.snps.ri.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   ind = col_character(),
    ##   ndata = col_double(),
    ##   nfiltered = col_double(),
    ##   nmiss = col_double(),
    ##   fmiss = col_double()
    ## )

``` r
meta <- read.delim("../meta/Sp_Radseq_Files_Guide.txt") %>% select(Urchin, Sex)
ind_miss_join <- left_join(ind_miss,meta,by=c("ind"="Urchin"))

ind_miss %>% arrange(desc(fmiss))
```

    ## # A tibble: 28 x 5
    ##    ind   ndata nfiltered nmiss fmiss
    ##    <chr> <dbl>     <dbl> <dbl> <dbl>
    ##  1 Spf42 10401         0  9113 0.876
    ##  2 Spf14 10401         0  6618 0.636
    ##  3 Spf45 10401         0  6222 0.598
    ##  4 Spf43 10401         0  4941 0.475
    ##  5 Spf13 10401         0  4629 0.445
    ##  6 Spf12 10401         0  3665 0.352
    ##  7 Spf26 10401         0  2775 0.267
    ##  8 Spf44 10401         0  2482 0.239
    ##  9 Spf6  10401         0  1217 0.117
    ## 10 Spf18 10401         0  1184 0.114
    ## # ... with 18 more rows

``` r
# Focus on the fmiss column telling us the proportion missing
s <- ggplot(ind_miss_join, aes(fmiss,fill=Sex)) + geom_histogram(colour = "black", alpha = 0.3) + theme_classic() + 
  ylab("# of Individuals") + xlab("Fraction Genotypes Missing")  + coord_cartesian(expand=F)

high_missing<-ind_miss %>%
  filter(fmiss > .5)

knitr::kable(
  high_missing, format = "html",caption = "Individuals Missing Greater than .4 Genotypes",
  booktabs = TRUE
)
```

<table>

<caption>

Individuals Missing Greater than .4 Genotypes

</caption>

<thead>

<tr>

<th style="text-align:left;">

ind

</th>

<th style="text-align:right;">

ndata

</th>

<th style="text-align:right;">

nfiltered

</th>

<th style="text-align:right;">

nmiss

</th>

<th style="text-align:right;">

fmiss

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Spf45

</td>

<td style="text-align:right;">

10401

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

6222

</td>

<td style="text-align:right;">

0.598212

</td>

</tr>

<tr>

<td style="text-align:left;">

Spf14

</td>

<td style="text-align:right;">

10401

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

6618

</td>

<td style="text-align:right;">

0.636285

</td>

</tr>

<tr>

<td style="text-align:left;">

Spf42

</td>

<td style="text-align:right;">

10401

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

9113

</td>

<td style="text-align:right;">

0.876166

</td>

</tr>

</tbody>

</table>

``` r
s
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](PLOTS/Missing_Individuals-1.png)<!-- -->

``` r
ind_depth <- read_delim("STATS/Spurp.minQ20.minGQ20.mac3.miss1.snps.ri.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1) %>% left_join(.,meta,by=c("ind"="Urchin"))
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   ind = col_character(),
    ##   nsites = col_double(),
    ##   depth = col_double()
    ## )

``` r
ind_depth %>% arrange(depth)
```

    ## # A tibble: 48 x 4
    ##    ind   nsites depth Sex   
    ##    <chr>  <dbl> <dbl> <chr> 
    ##  1 Spf14   3783  1.95 Female
    ##  2 Spf13   5772  2.18 Female
    ##  3 Spf26   7626  3.07 Female
    ##  4 Spf12   6736  4.06 Female
    ##  5 Spf18   9217  4.81 Female
    ##  6 Spf18   9217  4.81 Female
    ##  7 Spf25   9329  5.17 Female
    ##  8 Spf25   9329  5.17 Female
    ##  9 Spf6    9184  5.47 Female
    ## 10 Spf6    9184  5.47 Female
    ## # ... with 38 more rows

``` r
ind_depth %>% arrange(desc(depth))
```

    ## # A tibble: 48 x 4
    ##    ind   nsites depth Sex  
    ##    <chr>  <dbl> <dbl> <chr>
    ##  1 Spm3   10275  142. Male 
    ##  2 Spm3   10275  142. Male 
    ##  3 Spm4   10253  130. Male 
    ##  4 Spm4   10253  130. Male 
    ##  5 Spm8   10271  124. Male 
    ##  6 Spm8   10271  124. Male 
    ##  7 Spm16  10286  123. Male 
    ##  8 Spm16  10286  123. Male 
    ##  9 Spm11  10327  117. Male 
    ## 10 Spm11  10327  117. Male 
    ## # ... with 38 more rows

``` r
d <- ggplot(ind_depth, aes(depth,fill=Sex)) + geom_histogram(colour = "black", alpha = .3) + ylab("# of Individuals") + xlab("Mean Depth") + theme_classic() + coord_cartesian(expand=F)

d
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](PLOTS/Individual_Depth-1.png)<!-- -->

``` r
ind_het <- read_delim("STATS/Spurp_F2.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1) 

ind.het.plot <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_classic() +
  ylab("# of Individuals")+ xlab("Heterozygosity")  + coord_cartesian(expand=F)
ind_het %>% arrange(desc(f))
ind_het %>%
  filter(f < -.5 | f > .5)
```
