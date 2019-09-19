library(here)
library(tidyverse)
library(readxl)

## Extract bacterial qPCR standard curves

## Standard Curve 2016-09-19 ---------------------------------------------------
get_std_0919 <- function(xcl_file){
      cols <- c("well","sample_name","conc",
                "plate1", "skip5", "plate2", "skip7", "plate3",
                paste0("skip", 9:13))
      
      read_excel(path = xcl_file, sheet = "QDNA_20160919", 
                 skip = 3, col_names = cols, na = "Undetermined",) %>% 
            select(-starts_with("skip"))  %>%
            filter(sample_name %in% paste0("Std",1:7)) %>% 
            gather("plate","Ct",-well, -sample_name, -conc) %>%
            mutate(conc = as.numeric(conc), Ct = as.numeric(Ct),
                   std = "zymo", date = "2016-09-19")
}

## Standard Curve 2016-09-19 ---------------------------------------------------
get_std_1209 <- function(xcl_file){
      basic_cols <- c("well", "sample_name", "conc",
                      "plate1","plate2", "plate3")
      cols <- c("well", "sample_name", "conc",
                "plate1","skip5", "plate2", "skip7", "plate3", "skip9")

      full_cols <- c(paste("shan",cols,sep = "_"), "skip10",
                     paste("zymo",cols,sep = "_"))
      
      
      bac_std <- 
         read_excel(path = xcl_file,
                            sheet = "ReDo_QDNA_20161209",skip = 3, 
                            na = "Undetermined", col_names = full_cols) %>% 
            filter(shan_sample_name %in% paste0("Std",1:7)) %>% 
            select(-contains("skip"))
      
      shan_std <- bac_std %>% select(starts_with("shan")) %>% 
            set_colnames(basic_cols) %>% mutate(std = "shan")
      
      bac_std <- bac_std %>% select(starts_with("zymo")) %>% 
            set_colnames(basic_cols) %>% mutate(std = "zymo") %>% 
            bind_rows(shan_std)
      
      bac_std %>% gather("plate","Ct",-well, -sample_name, -conc, -std) %>% 
            mutate(conc = as.numeric(conc), Ct = as.numeric(Ct), date = "2016-12-09")
}

## Generating full dataset and caching
qpcrBacStd <- here("data","raw","MixStudy_Nate_20161209.xls") %>% 
      {full_join(get_std_0919(.), get_std_1209(.))}

saveRDS(qpcrBacStd, here("data","qpcrBacStd.RDS"))