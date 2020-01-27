---
title: "Clonal maintenance of transcriptional and epigenetic states in cancer cells"
author: "Zohar Meir"
date: "January 27, 2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Summary
This reposiroty contains three parts of analysis performed in Meir et al. - Expression, Dna methylation and Hic.
Each part consists of a hierarchy of smaller modules and steps.


## HTML page
In order to view the html go to https://tanaylab.github.io/Meir_et_al_nat_gen_2020_clonemem/ 


## Compilation
Each of the three R scripts (expression.r, methylation.r and hic.r) can run indepedently.
However:

Each of the three R scripts (expression.r, methylation.r and hic.r) can run indepedently.
However:
  - at each new session, start with sourcing the main functions script:

```r
source("Meir_et_al_2020_nat_gen_functions.r")
```

  - also, make sure you run the scripts from the same parent directory, maintaining folders of processed table (sub-directories expression_data/, methylation_data/ and HiC_data) just below it.

  - Only in part 3 (HiC.r), You will have to download large files in order to perform parts of the analysis (raw contacts map, in order to generate normalized SHAMAN maps of specific regions, or to compute insulation score yourself).,




a link to vignette:
https://tanaylab.github.io/Meir_et_al_nat_gen_2020_clonemem/scRNA_cellCycle_Vignette.html

