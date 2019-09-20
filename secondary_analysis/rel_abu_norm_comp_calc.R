library(here)
library(metagenomeSeq)
library(tidyverse)

mrexp_obj <- readRDS(here("data","mothur_mrexp.rds"))

normalize_counts <- function(method, mrexp){
    require(matrixStats)
    ## Extract count matrix from MRexperiment
    count_mat <- mrexp@assayData$counts 
    
    ## extract normalizaed counts
    if (method == "RAW") {
        ## Raw counts no normalization applied
        norm_factors = 1
    } else if (method == "UQ") {
        ## Upper quartile normalization
        count_mat[count_mat == 0] <- NA
        norm_factors <- colQuantiles(count_mat, p = 0.75 ,na.rm = TRUE)
    } else if (method == "CSS") {
        ## Cumulative sum scaling Paulson et al. 2013
        norm_factors <- metagenomeSeq::calcNormFactors(count_mat, p = 0.75)
        norm_factors <- norm_factors$normFactors
    } else if ( method == "TSS") {
        ## Total sum scaling 
        norm_factors <- colSums(count_mat)
    } else if (method %in% c("RLE","TMM")) {
        ## EdgeR RLE and TMM normalization methods
        norm_factors <- edgeR::calcNormFactors(count_mat, method = method) 
        norm_factors <- norm_factors * colSums(count_mat)
    } else {
        warning("Normalization method not defined")
    }
    
    ## Normalizing counts
    sweep(count_mat, 2, norm_factors,'/')
    
}

tidy_norm_counts <- function(norm_mat, mrexp){
    ## Tidy counts
    count_df <- norm_mat %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "feature_id") %>% 
        gather("id", "count", -feature_id)
    
    ## Adding sample data
    count_df <- pData(mrexp) %>% 
        select(biosample_id, id, t_fctr) %>% 
        right_join(count_df)
    
    ## Summarizing count values across PCR replicates 
    count_df %>% 
        group_by(biosample_id, feature_id, t_fctr) %>% 
        summarise(mean_count = mean(count),
                  var_count = var(count),
                  sd_count = sd(count))
}

## Normalizing counts using different methods/ normlization factors
norm_methods <- list(RAW = "RAW", ## No normalization
                     RLE = "RLE", ## EdgeR - relative log expression
                     TMM = "TMM", ## EdgeR - weighted trim mean of M-values
                     UQ = "UQ",   ## EdgeR - upperqurtile
                     CSS = "CSS", ## metagenomeSeq - cumulative sum scaling, with p = 0.75
                     TSS = "TSS") ## Total sum scaling (proportions)

norm_count_mats <- norm_methods %>% 
    map(normalize_counts, mrexp = mrexp_obj)

saveRDS(norm_count_mats, here("data","norm_count_mats.RDS"))

## Summarizing across PCR replicates and converting to a tidy data frame
norm_count_df <- norm_count_mats %>% 
    map(tidy_norm_counts, mrexp = mrexp_obj) %>% 
    map_df(bind_rows, .id = "norm_method") %>%
    mutate(pipe = "mothur")


saveRDS(norm_count_df, here("data","norm_count_df.RDS"))