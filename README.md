# Clonal maintenance of transcriptional and epigenetic states in cancer cells

## Summary
This reposiroty contains three parts of analysis performed in [Meir et al.](https://www.nature.com/articles/s41588-020-0645-y) - Expression, Dna methylation and Hic.
Each part consists of a hierarchy of smaller modules and steps.

## Compilation

Each of the three R scripts (expression.r, methylation.r and hic.r) can run indepedently.
However:

  - at each new session, start with sourcing the main functions script:

```r
source("Meir_et_al_2020_nat_gen_functions.r")
```

  - also, make sure you run the scripts from the same parent directory, maintaining folders of processed table (sub-directories expression_data/, methylation_data/ and HiC_data/) just below it.

  - Some of the analysis rely on tables you need to download. All functions that requires downloading will put the files in the proper subdirctory. Only expection in part 3 (HiC.r), where you will have to generate a "misha" database and download large files in order to perform parts of the analysis (raw contacts map, in order to generate normalized SHAMAN maps of specific regions, or to compute insulation score yourself). In this cases you will have to make sure you put the tracks in the misha-database subfolder after you download them.




a link to vignette examplifying logic of isolating cell-cycle independent variation in scRNA datasets:
https://tanaylab.github.io/Meir_et_al_nat_gen_2020_clonemem/scRNA_cellCycle_Normalization.html

