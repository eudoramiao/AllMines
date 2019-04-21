# Analysis of diversity of Mine sites


Contents
========
The repository contains prelimary analyses for mine sites, all the analyses can be found in [All mines analysis](AllMines_Graph.md)


Session Info
========
```
R version 3.5.3 (2019-03-11)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] stats4    parallel  grid      stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] ampvis_1.27.0               scales_1.0.0               
 [3] gridExtra_2.3               magrittr_1.5               
 [5] ggdendro_0.1-20             DESeq2_1.22.2              
 [7] SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
 [9] BiocParallel_1.16.6         matrixStats_0.54.0         
[11] Biobase_2.42.0              GenomicRanges_1.34.0       
[13] GenomeInfoDb_1.18.2         data.table_1.12.2          
[15] Biostrings_2.50.2           XVector_0.22.0             
[17] IRanges_2.16.0              S4Vectors_0.20.1           
[19] BiocGenerics_0.28.0         igraph_1.2.4               
[21] reshape2_1.4.3              jsonlite_1.6               
[23] rvest_0.3.3                 xml2_1.2.0                 
[25] ggrepel_0.8.0               kableExtra_1.1.0           
[27] knitr_1.22                  vegan_2.5-4                
[29] lattice_0.20-38             permute_0.9-5              
[31] metacoder_0.3.2             taxa_0.3.2                 
[33] cowplot_0.9.4               RColorBrewer_1.1-2         
[35] forcats_0.4.0               stringr_1.4.0              
[37] dplyr_0.8.0.1               purrr_0.3.2                
[39] readr_1.3.1                 tidyr_0.8.3                
[41] tibble_2.1.1                tidyverse_1.2.1            
[43] ggplot2_3.1.1               phyloseq_1.26.1            

loaded via a namespace (and not attached):
  [1] readxl_1.3.1           backports_1.1.4        Hmisc_4.2-0           
  [4] plyr_1.8.4             lazyeval_0.2.2         splines_3.5.3         
  [7] usethis_1.5.0          digest_0.6.18          foreach_1.4.4         
 [10] htmltools_0.3.6        fansi_0.4.0            checkmate_1.9.1       
 [13] memoise_1.1.0          cluster_2.0.8          remotes_2.0.4         
 [16] annotate_1.60.1        ggfittext_0.6.0        modelr_0.1.4          
 [19] prettyunits_1.0.2      colorspace_1.4-1       blob_1.1.1            
 [22] haven_2.1.0            xfun_0.6               callr_3.2.0           
 [25] crayon_1.3.4           RCurl_1.95-4.12        genefilter_1.64.0     
 [28] survival_2.44-1.1      iterators_1.0.10       ape_5.3               
 [31] glue_1.3.1             gtable_0.3.0           zlibbioc_1.28.0       
 [34] webshot_0.5.1          pkgbuild_1.0.3         Rhdf5lib_1.4.3        
 [37] DBI_1.0.0              Rcpp_1.0.1             viridisLite_0.3.0     
 [40] xtable_1.8-3           htmlTable_1.13.1       bit_1.1-14            
 [43] foreign_0.8-71         Formula_1.2-3          htmlwidgets_1.3       
 [46] httr_1.4.0             acepack_1.4.1          pkgconfig_2.0.2       
 [49] XML_3.98-1.19          nnet_7.3-12            locfit_1.5-9.1        
 [52] utf8_1.1.4             tidyselect_0.2.5       labeling_0.3          
 [55] rlang_0.3.4            AnnotationDbi_1.44.0   munsell_0.5.0         
 [58] cellranger_1.1.0       tools_3.5.3            cli_1.1.0             
 [61] RSQLite_2.1.1          generics_0.0.2         ade4_1.7-13           
 [64] devtools_2.0.2         broom_0.5.2            evaluate_0.13         
 [67] biomformat_1.10.1      yaml_2.2.0             bit64_0.9-7           
 [70] processx_3.3.0         fs_1.2.7               nlme_3.1-139          
 [73] GA_3.2                 compiler_3.5.3         rstudioapi_0.10       
 [76] curl_3.3               geneplotter_1.60.0     stringi_1.4.3         
 [79] ps_1.3.0               desc_1.2.0             Matrix_1.2-17         
 [82] multtest_2.38.0        pillar_1.3.1           BiocManager_1.30.4    
 [85] bitops_1.0-6           R6_2.4.0               latticeExtra_0.6-28   
 [88] sessioninfo_1.1.1      codetools_0.2-16       MASS_7.3-51.3         
 [91] assertthat_0.2.1       pkgload_1.0.2          rhdf5_2.26.2          
 [94] rprojroot_1.3-2        withr_2.1.2            GenomeInfoDbData_1.2.0
 [97] mgcv_1.8-28            hms_0.4.2              rpart_4.1-15          
[100] rmarkdown_1.12         lubridate_1.7.4        base64enc_0.1-3
```