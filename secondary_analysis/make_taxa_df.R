## Creating taxa_df for set and relative abundance plots
library(here)
library(phyloseq)
library(tidyverse)

ps_list <- list(
    dada2 = here("data","dada2_ps.rds"),
    mothur =  here("data","mothur_ps.rds"),
    qiime =  here("data","qiime_ps.rds")
) %>% 
    map(readRDS)

get_taxa_df <- function(ps){
    taxa_df <- tax_table(ps) %>% 
        as.data.frame() %>%
        tibble::rownames_to_column(var = "feature_id") %>% 
        mutate_all(funs(str_remove(.,".__")))
    
    taxa_df <- tibble(feature_id = taxa_names(ps), 
                      f_counts = taxa_sums(ps)) %>% 
        left_join(taxa_df) 
    
    taxa_df %>% 
        group_by(Rank1, Rank2, Rank3, Rank4, Rank5, Rank6) %>% 
        summarise(total_count = sum(f_counts))
}

taxa_df <- ps_list %>% map_dfr(get_taxa_df, .id = "pipe")

taxa_df <- taxa_df %>% 
    group_by(pipe) %>% 
    mutate(rel_abu = total_count / sum(total_count))

saveRDS(taxa_df, here("data","taxa_df.RDS"))