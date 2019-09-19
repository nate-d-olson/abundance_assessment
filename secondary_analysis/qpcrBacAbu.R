library(here)
library(tidyverse)
library(readxl)

## Munge Bacterial Abundance qPCR data
get_abu_0909 <- function(xcl_file){
      cols <- c("well","sample_name","conc",
                "plate1_Ct","plate1_quant",
                "plate2_Ct","plate2_quant",
                "plate3_Ct","plate3_quant")
      
      read_excel(path = xcl_file, sheet = "QDNA_20160919", 
                 skip = 11, na = "Undetermined", col_names = cols) %>% 
            mutate(std = "zymo", date = "2016-09-19") 
}

get_abu_1209 <- function(xcl_file){
      cols <- c("well","sample_name","conc",
                "plate1_Ct","plate1_quant",
                "plate2_Ct","plate2_quant",
                "plate3_Ct","plate3_quant")
      
      std_cols <- cols %>% {c(paste("shan", ., sep = "_"), 
                              "blank",
                              paste("zymo", ., sep = "_"))}
      
      bac_con_raw <- read_excel(path = xcl_file,
                                sheet = "ReDo_QDNA_20161209",
                                skip = 11, na = "Undetermined", 
                                col_names = std_cols)
      
      shan_con_raw <- bac_con_raw %>%
         select(starts_with("shan")) %>%
         set_colnames(cols) %>%
         mutate(std = "shan")
      
      bac_con_raw %>%
         select(starts_with("zymo")) %>%
         set_colnames(cols) %>%
         mutate(std = "zymo") %>%
         bind_rows(shan_con_raw) %>%
         mutate(date = "2016-12-09")
}

tidy_abu_con <- function(bac_con_raw){
      bac_con <- bac_con_raw %>% 
            gather("id","value", -well, -sample_name, -std, -date) %>% 
            separate(id, c("plate","var"), sep = "_") %>%
            spread(var,value) %>%
            mutate(sam_type = if_else(grepl('\\(', sample_name), "unmixed","titration"),
                   sam_type = if_else(sample_name == "NTC", "NTC",sam_type)) %>%
            rename(stine_quant = quant) %>%
            filter(!(well %in% c("Plate", paste0("P",1:3))))
      
      bac_unmixed <- bac_con %>% filter(sam_type == "unmixed") %>% 
            mutate(sample_name = str_replace(sample_name, ".*_", ""), 
                   sample_name = str_replace(sample_name, '\\('," "),
                   sample_name = str_replace(sample_name, '\\)',""))
      
      bac_con %>% filter(sam_type != "unmixed") %>% bind_rows(bac_unmixed)
}

## QPCR results file
xcl_file <- here("data","raw","MixStudy_Nate_20161209.xls")

## Parsing results for two dates
abu_909 <- get_abu_0909(xcl_file)
abu_129 <- get_abu_1209(xcl_file)

## Combining and tidying data frames
bac_con_raw <- bind_rows(abu_909, abu_129) %>% select(-conc)
qpcrBacAbu <- tidy_abu_con(bac_con_raw)

saveRDS(qpcrBacAbu, here("data", "qpcrBacAbu.RDS"))