## Extract ERCC qPCR data
library(here)
library(tidyverse)
library(readxl)

## Standard Curve 2016-09-19 ---------------------------------------------------
get_ercc <- function(xcl_file){
      cols <- c("ercc_plate","ercc_std",
                paste0("std_rep",rep(1:3, each = 2),"_",rep(c("Ct","quant"),3)),
                "sampleID",
                paste0("sam_rep",rep(1:3, each = 2),"_",rep(c("Ct","quant"),3)))
      
      ## ERCC plasmid id
      ercc_assays <- rep(c(84,12,34,157,57,108,130,2,92,35), each = 8)
      
      read_excel(path = xcl_file, sheet = "ERCC_Quant_20161208",
                 skip = 3, col_names = cols, na = "Undetermined") %>% 
         filter(!is.na(ercc_plate)) %>% 
            add_column(ercc = ercc_assays) %>% 
            gather("key","value",-ercc_plate,-ercc,-ercc_std, -sampleID) %>% 
            mutate(sampleID = if_else(grepl("std",key), ercc_std, sampleID)) %>% 
            select(-ercc_std) %>% separate(key, c("sample_type","rep","key"), sep = "_") %>% 
            spread(key,value) %>% 
            mutate(Ct = as.numeric(Ct), quant = as.numeric(quant))
}

here("data","raw", "MixStudy_Nate_20161209.xls") %>% 
   get_ercc() %>% 
   saveRDS(here("data","qpcrERCC.RDS"))


# ERCC Spike-in ----------------------------------------------------------------

## loading source data files
get_erccMeta <- function(){
   ercc_assay <- read_tsv(here("data","raw","thermo_fisher_ercc_qPCR_assays.txt"))
   ercc_spike <- read_tsv(here("data","raw","ERCC-Spike-ins.tsv"),comment = "#") %>%
      rename(`Gene Symbol` = Control)
   
   
   ## Assays selected to minimize differences in amplicon length between the pre
   ## and post qPCR assays.
   selected_assays <- c("Ac03459877_a1", "Ac03459922_a1", "Ac03459958_a1",
                        "Ac03459987_a1", "Ac03460000_a1","Ac03460028_a1",
                        "Ac03460000_a1","Ac03460028_a1","Ac03459872_a1",
                        "Ac03460039_a1","Ac03459892_a1","Ac03459925_a1")
   
   ercc_spike_ins <- ercc_assay %>%
      select(`Gene Symbol`, `Assay ID`,`Amplicon Length`) %>%
      filter(`Assay ID` %in% selected_assays) %>% full_join(ercc_spike, .) %>%
      rename(ercc_id = `Gene Symbol`,
             biosample_id = Sample_spike,
             treatment = Treatment,
             assay_id = `Assay ID`,
             GC = `GC (excluding polyA tail)`,
             amplicon_length = `Amplicon Length`)
   
   write_csv(ercc_spike_ins,here("data","raw","ercc_spike_ins.csv"))
   ercc_spike_ins
}

get_erccMeta() %>% saveRDS(here("data","erccMeta.RDS"))
