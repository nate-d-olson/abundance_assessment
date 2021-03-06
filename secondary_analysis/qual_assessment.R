## ------------------------------------------------------------------------
library(here)
library(tidyverse)
library(broom)
library(ggpubr)


## Initiating parallelization 
library(furrr)
plan(multiprocess, workers = 6, gc = TRUE)

## input data files
raw_counts_file <- here("data","raw_counts.RDS")
nb_counts_file <- here("data","nb_counts_titrations.RDS")
pa_summary_anno_file <- here("data","pa_summary_anno_df.RDS")
theta_est_file <- here("data","bootstrap_theta_estimates.RDS")


## ----qual_tidy-----------------------------------------------------------
## Feature-PCR rep level tidy data frame
## Relevant raw_count factors
count_df <- readRDS(raw_counts_file) %>%
      select(biosample_id, t_fctr, id, pipe, feature_id, count) %>% 
      # grouping rows by `id` column and then adding column of total_counts
      group_by(id) %>%
      mutate(total_count = sum(count))

## Sanity Check
count_df$total_count %>% unique() %>% summary()

## Loading negative binomial count data
nb_counts <- readRDS(nb_counts_file)

## Annotating with feature type
# essentially just reordered cols of count_df and changed what is being "grouped" by
count_df <- readRDS(pa_summary_anno_file) %>%
      select(biosample_id, pipe, feature_id, pa_specific) %>%
      left_join(count_df)

## Inferred Feature Proportions
### Same code used in Rel_Abu_Error
pre_post_prop <- nb_counts %>%
      ungroup() %>%
      filter(t_fctr %in% c(0,20)) %>% #this grabs all of the unmixed samples
      mutate(end_point = if_else(t_fctr == 0 , "post", "pre")) %>%
      select(-t_fctr) %>%
      ## setting values to 0 when one or more of the PCR replicates are 0 for titration end-points
      spread(end_point, nb_prop, fill = 0)

prop_inferred <- readRDS(theta_est_file) %>%
      filter(pipe == "unclustered") %>%
      ungroup() %>%
      mutate(t_fctr = factor(t_fctr, levels = c(0:5, 10, 15, 20))) %>%
      select(biosample_id, theta_hat_mean, t_fctr) %>%
      right_join(nb_counts) %>% 
      right_join(pre_post_prop) %>%
      filter(t_fctr %in% c(1:5,10,15)) %>%
      ## Setting negative inferred theta values to 0
      mutate(theta_hat_mean = if_else(theta_hat_mean < 0, 
                                      0, theta_hat_mean)) %>%
      ## Using inferred theta estimates to calculate expected values
      mutate(inferred_prop = post * theta_hat_mean + pre * (1-theta_hat_mean))

count_df <- prop_inferred %>%
      select(biosample_id, t_fctr, pipe, feature_id,
             pre, post, inferred_prop, theta_hat_mean) %>%
      left_join(count_df)

unmix_count_df <- filter(count_df,
                    t_fctr %in% c(1:5, 10, 15),
                    pa_specific == "unmixed",
                    # should not have to filter by inferred_prop unmix should
                    # have non-zero inferred prop or at least non-zero pre +
                    # post
                    pre + post != 0, count == 0)

mix_count_df <- filter(count_df,
                  t_fctr %in% c(1:5, 10, 15),
                  count != 0, pre == 0, post == 0)


## ------------------------------------------------------------------------
unmix_binom_test <- function(unmix_count_df, output_file, confidence){
    if(file.exists(output_file)){
        return(readRDS(output_file))
    } else {
        print("could not find results file in 'data' ... running code now")
    }
    ## Binomial test - Is the true proportion of a feature within a mixture less
    ## than the inferred proportion assuming no observed counts for a titration
    ## given the number of reads sequenced for a PCR replicate.
      
    unmix_binom_test_df <- unmix_count_df %>%
        group_by(pipe, biosample_id, id, feature_id, t_fctr) %>%
        nest() %>%                                                
        mutate(binom_test = future_map(data, 
                                       ~with(., 
                                             binom.test(
                                               x = 0,  
                                               n = total_count , 
                                               p = inferred_prop, 
                                               alternative = "less",
                                               conf.level = confidence
                                               )
                                             ),
                                       .progress = TRUE
                                       ),
               binom_tidy = map(binom_test, tidy))

    ## tidy-ing up results
    unmix_binom_tidy <- unmix_binom_test_df %>%
        select(pipe, biosample_id, id, feature_id, t_fctr, binom_tidy) %>%
        unnest() 
        ## multiple test correction
        mutate(adj.pvalue = p.adjust(p.value, method = "BH"))

    ## Saving data frame as RDS to prevent having to rerun
    saveRDS(unmix_binom_tidy, output_file)
}


## ------------------------------------------------------------------------
unmix_binom_0.95 <- unmix_binom_test(
  unmix_count_df,
  output_file = here("data", "unmix_binom_test.RDS"),
  confidence = 0.95
)

## TODO - figure out which one to use
unmix_binom_0.975 <- unmix_binom_test(
  unmix_count_df,
  output_file = here("data", "unmix_binom_0.975.RDS"),
  confidence = 0.975
)


## ----echo=TRUE, message = FALSE, warning=FALSE, include=FALSE------------
calc_bayes_mc <- function(max_prop, total_count, obs_count = count, ALPHA, BETA){
      ## Using sapply to simulate join beta/uniform and binomal distributions
      
      # changing from simulating samples from a uniform distribution to
      # a beta distribution while varrying parameter selection on the
      # beta distribution
      
      # for pi_obs >= max_prop, switching to a beta distribution and
      # iterating over parameters alpha (shape1) and beta (shape2)
      prop_sim <- rbeta(n = 10000, shape1 = ALPHA, shape2 = BETA)
      mc_gt_count <- sapply(prop_sim, rbinom, size = total_count, n = 1)
      p_gt_max <- sum(mc_gt_count >= obs_count)/length(mc_gt_count) 
      
      ## sim distribution for P(Pi <= max_prop)
      ## TODO - ask if this should also be beta distributed???
      ## - I think it should be prop_sim for 0 to 1 distribution then split sims by < count and >= count
      prop_sim <- runif(n = 10000, min = 0, max = max_prop)
      mc_lt_count <- sapply(prop_sim, rbinom, size = total_count, n = 1)
      p_lt_max <- sum(mc_lt_count >= obs_count)/length(mc_lt_count)
      
      ## Bayesian Hypothesis test
      (p_lt_max * 0.5)/sum((p_lt_max * 0.5), (p_gt_max * 0.5))
}

mix_binom_bayes <- function(output_file, count_df, mix_count_df){
      ALPHA <- str_extract(outfiles, "(?<=test_).*(?=_)") %>% 
        as.numeric()
      
      BETA <- str_extract(outfiles, "(?<=[:digit:]_).*(?=.RDS)") %>% 
        as.numeric()
      
      ## Minimum expected proportion
      min_prop <- count_df %>% 
            ## Only allowing for non-zero minimum proportions
            ## This potentially results in a higher minimum expected proprotion than
            ## allowing for one of the two unmixed samples to have a 0 minimum
            ## proportion.
            mutate(pre = if_else(pre == 0, 1, pre),
                   post = if_else(post == 0, 1, post)) %>% 
            group_by(pipe, biosample_id) %>% 
            ## For some samples with the minimum expected prop is 0
            mutate(min_pre = min(pre),
                   min_post = min(post)) %>% 
            select(pipe, biosample_id, theta_hat_mean, min_pre, min_post) %>% 
            unique() %>% 
            mutate(min_exp_prop = min_post * theta_hat_mean + min_pre * (1 - theta_hat_mean))
      
      mix_count<- mix_count_df %>% 
            filter(pa_specific == "mixed") %>% 
            left_join(min_prop) %>% 
            ## excluding Unclustered decreases the number of tests from 136k to 28k
            filter(pipe != "unclustered")
      
      mix_binom_test <- mix_count %>% 
            ## Minimizing the number of test performed
            select(min_exp_prop, total_count, count) %>% unique() %>% 
            rowwise() %>%
            mutate(bayes_test = calc_bayes_mc(max_prop = min_exp_prop, 
                                              total_count = total_count, 
                                              obs_count = count,
                                              ALPHA, BETA))
      
      mix_binom_test_df <- mix_count %>% 
        left_join(mix_binom_test) %>% 
        mutate(adj.p.value = p.adjust(bayes_test, method = "BH"))

      saveRDS(mix_binom_test_df, here("data", "mix_bayes_res", output_file))
}


## ------------------------------------------------------------------------
## alpha and beta parameter values
model_values <- str_pad(c(1, 2, 10, 25), width = 2, side = "left", pad = "0")
model_params <- cross2(model_values, model_values)

outfiles <- model_params %>% 
  map_chr(paste, collapse = "_") %>% 
  {paste0("test_", ., ".RDS")}

## Defining output files
model_params %>% 
  map_chr(paste, collapse = "_") %>% 
  {paste0("test_", ., ".RDS")} %>% 
  ## performing MC tests
  future_map(mix_binom_bayes, 
             count_df, 
             mix_count_df,
             .progress = TRUE)

