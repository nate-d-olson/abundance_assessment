## Functions for validating and fixing sample sheet
library(tidyverse)

check_samplesheet <- function(sample_df, fix = TRUE){
    ## Add columns with correct biosample (b_check) and titration (t_check)
    t20_pos <- paste0("A", c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11))
    t0_pos <- paste0(rep(c("B", "C", "E", "F", "G"), each = 2), 
                     rep(c(6, 12), 5))
    ## Setting titration factor as character to prevent comparison errors
    sample_df <- sample_df %>% 
        mutate(t_fctr = as.character(t_fctr),
               biosample_id = as.character(biosample_id))
    
    sam_check_df <- sample_df %>%
        mutate(t_fctr = as.character(t_fctr),
            t_check = case_when(pos %in% t20_pos ~ "20",
                                pos %in% t0_pos ~ "0",
                                TRUE ~ t_fctr
            ),
            b_check = case_when(
                pos %in% c("B6", "B12") ~ "E01JH0004",
                pos %in% c("C6", "C12") ~ "E01JH0011",
                pos %in% c("E6", "E12") ~ "E01JH0016",
                pos %in% c("F6", "F12") ~ "E01JH0017",
                pos %in% c("G6", "G12") ~ "E01JH0038",
                TRUE ~ biosample_id
            )
        )
    if (all(sam_check_df$t_fctr == sam_check_df$t_check,na.rm = TRUE) |
        all(sam_check_df$biosample_id == sam_check_df$b_check, na.rm = TRUE)) {
        if (fix) {
            msg <- str_c(
                "Titration or biosample metadata is not corrected. ",
                "Returning data frame with corrected metadata"
            )
            message(msg)
            
            corrected_sample_df <- sam_check_df %>%
                mutate(t_fctr = t_check,
                       biosample_id = b_check) %>%
                select(-t_check, -b_check)
            
            return(corrected_sample_df)
        } else {
            message(str_c(
                "Titration or biosample metadata is not correct ",
                "and has not been corrected. Returning input data frame."
            ))
            return(sample_df)
        }
    }
    
    message(str_c("Titration and biosample metadata is correct", 
            ", returning input data frame"))
    return(sample_df)
}

unfix_sample_sheet <- function(sample_df){
        ## Set biosample_id to experimental design values
       sample_df %>%
        mutate(
            biosample_id = case_when(
                pos %in% c("B6", "B12") ~ "E01JH0011",
                pos %in% c("C6", "C12") ~ "E01JH0016",
                pos %in% c("E6", "E12") ~ "E01JH0017",
                pos %in% c("F6", "F12") ~ "E01JH0004",
                pos %in% c("G6", "G12") ~ "E01JH0038",
                TRUE ~ biosample_id
            )
        )
}

## Sanity checks for sample check functions    
# ps <- readRDS("data/dada2_ps.rds")
# sample_df <- sample_data(ps) %>% as.data.frame()
# df <- check_samplesheet(sample_df, fix = FALSE)
# df <- check_samplesheet(sample_df, fix = TRUE)
# incorrect_df <- unfix_sample_sheet(sample_df)
# uncorrected_df <- check_samplesheet(incorrect_df, fix = FALSE)
# corrected_df <- check_samplesheet(incorrect_df, fix = TRUE)
# 
# all(corrected_df == sample_df, na.rm = TRUE)
