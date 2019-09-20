library(here)
require(tidyverse)
require(metagenomeSeq)
mrexp_obj <- here("data","mothur_mrexp.rds") %>% readRDS()
norm_count_mats <- here("data","norm_count_mats.RDS") %>% readRDS()

make_titration_comp_df <- function(mrexp_obj){
    pdat <- pData(mrexp_obj)  %>% 
        mutate(t_fctr = factor(t_fctr, level = c(0:5,10,15,20)))
    t1_pdat <- pdat %>% dplyr::select(biosample_id, id, t_fctr) %>% 
        dplyr::rename(T1 = t_fctr, T1_id = id) 
    t2_pdat <- pdat %>% dplyr::select(biosample_id, id, t_fctr) %>% 
        dplyr::rename(T2 = t_fctr, T2_id = id) 
    
    ## Generate titration comparison sets and reduce redundancy
    left_join(t1_pdat, t2_pdat) %>% 
        filter(as.numeric(T1) < as.numeric(T2)) %>% 
        group_by(biosample_id, T1, T2) %>% 
        summarise(sam_names = c(T1_id,T2_id) %>% unique() %>% list(),
                  n_sams = c(T1_id,T2_id) %>% unique() %>% length())
}

make_titration_comp_subset_df <- function(norm_mat, titration_comp_df, mrexp){
    
    ## Subsetting mrexp for titrations
    subset_mat <- function(sam_names){
        sub_mat <- norm_mat[,which(colnames(norm_mat) %in% sam_names)]
        ## Remove features with no counts for all samples 
        sub_mat <- sub_mat[rowSums(sub_mat) > 0,]
    }
    
    get_design_mat <- function(sam_names){
        mrexp_obj[,which(colnames(mrexp_obj) %in% sam_names)] %>% 
            pData() %>% model.matrix(~t_fctr, data = .)
    }
    
    titration_comp_df %>% 
        ## Dataframe with list of subsetted mrexp
        mutate(norm_mat_sub = map(sam_names, subset_mat),
               design_mat = map(sam_names, get_design_mat))
}

calc_logFC <- function(subset_mat, design_mat){
    require(limma)
    fit <- lmFit(object = subset_mat,design = design_mat)
    fit <- eBayes(fit)
    
    ## logFC and test results
    topTable(fit,number = Inf) %>% 
        rownames_to_column(var = "feature_id")
}

calc_titration_comp_logFC <- function(norm_count_mat, titation_comp_df){
    titration_comp_df <- make_titration_comp_subset_df(norm_count_mat, titration_comp_df)  
    
    titration_comp_logFC_df <- titration_comp_df %>% 
        mutate(logFC_df = map2(norm_mat_sub, design_mat, calc_logFC))
}

## Calculate logFC for titration comparisons for multiple normalization methods 
titration_comp_df <- make_titration_comp_df(mrexp_obj)

norm_logFC <- norm_count_mats %>% map(calc_titration_comp_logFC, titration_comp_df)
saveRDS(norm_logFC, here("data","norm_logFC.RDS"))

## Generate tidy data frame with logFC results
norm_logFC_df <- norm_logFC %>% 
    bind_rows(.id = "norm_method") %>%
    select(biosample_id, norm_method, T1, T2, logFC_df) %>%
    mutate(pipe = "mothur") %>% 
    unnest()

saveRDS(norm_logFC_df, here("data","norm_logFC_df.RDS"))