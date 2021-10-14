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
ind_miss  <- read_delim("STATS/Spurp.minQ20.minGQ20.minDP40.maxDP166.miss1.snps.ri.imiss", delim = "\t",
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
    ##  1 Spf42 19186         0 16800 0.876
    ##  2 Spf14 19186         0 11734 0.612
    ##  3 Spf45 19186         0 10885 0.567
    ##  4 Spf43 19186         0  8593 0.448
    ##  5 Spf13 19186         0  7980 0.416
    ##  6 Spf12 19186         0  6313 0.329
    ##  7 Spf26 19186         0  4672 0.244
    ##  8 Spf44 19186         0  4386 0.229
    ##  9 Spf6  19186         0  1980 0.103
    ## 10 Spf18 19186         0  1920 0.100
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

19186

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

10885

</td>

<td style="text-align:right;">

0.567341

</td>

</tr>

<tr>

<td style="text-align:left;">

Spf14

</td>

<td style="text-align:right;">

19186

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

11734

</td>

<td style="text-align:right;">

0.611592

</td>

</tr>

<tr>

<td style="text-align:left;">

Spf42

</td>

<td style="text-align:right;">

19186

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

16800

</td>

<td style="text-align:right;">

0.875638

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
ind_depth <- read_delim("STATS/Spurp.minQ20.minGQ20.minDP40.maxDP166.miss1.snps.ri.idepth", delim = "\t",
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
    ##  1 Spf42   2386  1.75 Female
    ##  2 Spf14   7452  1.85 Female
    ##  3 Spf13  11206  2.20 Female
    ##  4 Spf45   8301  2.54 Female
    ##  5 Spf43  10593  2.82 Female
    ##  6 Spf12  12873  2.97 Female
    ##  7 Spf26  14514  2.98 Female
    ##  8 Spf44  14800  4.45 Female
    ##  9 Spf18  17266  4.84 Female
    ## 10 Spf18  17266  4.84 Female
    ## # ... with 38 more rows

``` r
ind_depth %>% arrange(desc(depth))
```

    ## # A tibble: 48 x 4
    ##    ind   nsites depth Sex  
    ##    <chr>  <dbl> <dbl> <chr>
    ##  1 Spm3   19007  157. Male 
    ##  2 Spm3   19007  157. Male 
    ##  3 Spm4   19011  144. Male 
    ##  4 Spm4   19011  144. Male 
    ##  5 Spm8   19004  139. Male 
    ##  6 Spm8   19004  139. Male 
    ##  7 Spm16  18981  136. Male 
    ##  8 Spm16  18981  136. Male 
    ##  9 Spm11  19080  129. Male 
    ## 10 Spm11  19080  129. Male 
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
