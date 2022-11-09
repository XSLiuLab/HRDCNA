# HRD<sub>CNA</sub>: Homologous recombination deficiency prediction by copy number alteration features

## Introduction

HRD<sub>CNA</sub> is a robust HRD predictor based on copy number alteration (CNA) features. CNA information can be obtained from a diverse type of data, such as shallow WGS, WES, SNP array, and panel sequencing, and could represent a cost-effective type of biomarker for cancer diagnosis and clinical response prediction. HRD<sub>CNA</sub> can precisely predict HR status across cancer types using CNA features data derived from different platforms and it provides a robust tool for cost-effective HRD prediction and also demonstrates the applicability of CNA features in cancer precision medicine. We made HRD<sub>CNA</sub> into an R package `HRDCNA` for use.

The package provides you with the functions and example data to automise this process of extracting CNA features and calculating HRD<sub>CNA</sub> scores.

## Getting started

### Requirements

-   Software: R
-   Operating system: Linux, OS X, Windows
-   R version: 4.1.0

### Installation
You can install the development version of `HRDCNA` from Github with:
``` r
# From GitHub
install.packages("remotes")
remotes::install_github("XSLiuLab/HRDCNA")
library(HRDCNA)
```

### Extracting CNA features

#### Input file
The `HRDCNA` input requires absolute copy number profile with following information:

-   Segment chromosome.
-   Segment start.
-   Segment end.
-   Absolute copy number value for this segment.
-   Sample ID.

The input data can be result from any software which provides information above.
Useful softwares are listed below:

- [ABSOLUTE](https://software.broadinstitute.org/cancer/cga/absolute)
- [Sequenza](https://cran.r-project.org/web/packages/sequenza/index.html)
- [FACETS](https://github.com/mskcc/facets)
- [PennCNV](https://penncnv.openbioinformatics.org/en/latest/) & [ASCAT](https://www.crick.ac.uk/research/labs/peter-van-loo/software)
- [CNVkit](https://github.com/etal/cnvkit)

``` r
cn_wgs <- readRDS("./data/test/cn_60_wgs.rds")
head(cn_wgs)
```
```
#      sample chromosome    start      end segVal
# 1 FBC020030          1    13116  1598432      2
# 2 FBC020030          1  1599547  1661844      1
# 3 FBC020030          1  1662895 11043464      2
# 4 FBC020030          1 11044593 28695187      3
# 5 FBC020030          1 28696392 28747637      2
# 6 FBC020030          1 28749014 29497250      3
```

#### Extracting CNA features
```r
nmfcn_wgs <- sigminercopy(data = cn_wgs, "hg19")
```
```
# ℹ [2022-11-09 05:06:01]: Started.
# ℹ [2022-11-09 05:06:01]: Genome build  : hg19.
# ℹ [2022-11-09 05:06:01]: Genome measure: called.
# ✔ [2022-11-09 05:06:01]: Chromosome size database for build obtained.
# ℹ [2022-11-09 05:06:01]: Reading input.
# ✔ [2022-11-09 05:06:01]: A data frame as input detected.
# ✔ [2022-11-09 05:06:01]: Column names checked.
# ✔ [2022-11-09 05:06:01]: Column order set.
# ✔ [2022-11-09 05:06:02]: Chromosomes unified.
# ✔ [2022-11-09 05:06:02]: Data imported.
# ℹ [2022-11-09 05:06:02]: Segments info:
# ℹ [2022-11-09 05:06:02]:     Keep - 25915
# ℹ [2022-11-09 05:06:02]:     Drop - 0
# ✔ [2022-11-09 05:06:02]: Segments sorted.
# ℹ [2022-11-09 05:06:02]: Joining adjacent segments with same copy number value. Be patient...
# ✔ [2022-11-09 05:06:06]: 24809 segments left after joining.
# ✔ [2022-11-09 05:06:06]: Segmental table cleaned.
# ℹ [2022-11-09 05:06:06]: Annotating.
# ✔ [2022-11-09 05:06:07]: Annotation done.
# ℹ [2022-11-09 05:06:07]: Summarizing per sample.
# ✔ [2022-11-09 05:06:07]: Summarized.
# ℹ [2022-11-09 05:06:07]: Generating CopyNumber object.
# ✔ [2022-11-09 05:06:07]: Generated.
# ℹ [2022-11-09 05:06:07]: Validating object.
# ✔ [2022-11-09 05:06:07]: Done.
# ℹ [2022-11-09 05:06:07]: 6.438 secs elapsed.
# ℹ [2022-11-09 05:06:07]: Started.
# ℹ [2022-11-09 05:06:08]: Step: getting copy number features.
# ℹ [2022-11-09 05:06:08]: Getting breakpoint count per 10 Mb...
# ℹ [2022-11-09 05:06:14]: Getting breakpoint count per chromosome arm...
# ℹ [2022-11-09 05:06:20]: Getting copy number...
# ℹ [2022-11-09 05:06:20]: Getting change-point copy number change...
# ℹ [2022-11-09 05:06:25]: Getting length of chains of oscillating copy number...
# ℹ [2022-11-09 05:06:29]: Getting (log10 based) segment size...
# ℹ [2022-11-09 05:06:29]: Getting the minimal number of chromosome with 50% CNV...
# ℹ [2022-11-09 05:06:32]: Getting burden of chromosome...
# ✔ [2022-11-09 05:06:33]: Gotten.
# ℹ [2022-11-09 05:06:33]: Step: generating copy number components.
# ✔ [2022-11-09 05:06:33]: `feature_setting` checked.
# ℹ [2022-11-09 05:06:33]: Step: counting components.
# ✔ [2022-11-09 05:06:35]: Counted.
# ℹ [2022-11-09 05:06:35]: Step: generating components by sample matrix.
# ✔ [2022-11-09 05:06:35]: Matrix generated.
# ℹ [2022-11-09 05:06:35]: 27.582 secs elapsed.
```

This step returns an NMF matrix containing information about CNA features and their counts. Below we can see what a part of the NMF matrix looks like.
```r
head(nmfcn_wgs)
```
```
  BP10MB[0] BP10MB[1] BP10MB[2] BP10MB[3] BP10MB[4] BP10MB[5] BP10MB[>5] BPArm[0]
1       246        29        31         4         3         1          2       11
2       262        30        19         4         0         1          0       19
3       147        54        50        23        21        12          9        4
4       222        59        26         7         2         0          0        8
5       175        61        46        18         9         4          3        5
6       127        72        49        32        14         9         13        4
```

### Calculating HRD<sub>CNA</sub> score
Once we have the NMF matrix containing information about CNA features and their counts, HRD<sub>CNA</sub> score can be calculated and we can use it for predicting HRD.
```r
score_wgs <- HRDprediction(data = nmfcn_wgs)
head(score_wgs)
```
```
#       HRDCNAScore    sample
# 1 0.06710088 FBC013587
# 2 0.05277765 FBC016006
# 3 0.97387076 FBC016026
# 4 0.98887292 FBC016050
# 5 0.97191659 FBC020021
# 6 0.96064644 FBC020030
```

The higher the HRD<sub>CNA</sub> score, the greater the probability that the sample is HRD.

The development process of HRD<sub>CNA</sub> model, its applications in biology, and generated data and figures can be achieved can be read online at [InterpretationAnalysisHRDCNA](https://github.com/XSLiuLab/InterpretationAnalysisHRDCNA).


