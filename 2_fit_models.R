# dependencies
library(rethinking)
library(rstan)
library(viridis)
library(cmdstanr)
library(parallel)

d <- read.csv("s_person_year_df.csv", stringsAsFactors = FALSE)

person_ids <- sort(unique(d$person_id))

# select subset
set.seed(1)
n_rp <- 100
# N = 36595 in the full sample

rp_sub <- sample(person_ids, size = n_rp)

dm <- subset(d, d$person_id %in% rp_sub)
# re-index person_id for the stats model
person_ids <- sort(unique(dm$person_id))
dm$person_id <- match(dm$person_id, person_ids)

age_list <- sort(unique(dm$age))

# re-index age for model; real age 0 (last birthday) is now age_bin[1]
dm$age_bin <- dm$age + 1
# note that this sort of assumes that the *minimum* age is age 0, but is that always true??

# distance matrix for ages
d_mat <- as.matrix(dist(age_list, upper = TRUE, diag = TRUE))

# data list
data <- list(N_ages = length(age_list),
             N_ind = length(person_ids),
             N = nrow(dm), 
             y = dm$n_moves,
             age = dm$age_bin,
             person_id = dm$person_id,
             d_mat = d_mat,
             run_estimation = 0)

m_pois <- cmdstan_model("m_pois.stan")

# simulate from the stan mdoel
d_pois_sim <- m_pois$sample(
  data = data,
  chains = 1,
  iter_sampling = 2,
  fixed_param = TRUE)

stansim_pois <- rstan::read_stan_csv(d_pois_sim$output_files())
dat_pois <- extract.samples(stansim_pois)

data <- list(N_ages = length(age_list),
             N_ind = length(person_ids),
             N = nrow(dm), 
             y = dat_pois$y_sim[1,], # pass a simulation as input into the stats model
             age = dm$age_bin,
             person_id = dm$person_id,
             d_mat = d_mat,
             run_estimation = 1)

# fit poisson
stanfit_pois <- cstan(file = "m_pois.stan", 
                      data = data,
                      chains = 4, 
                      cores = 60, 
                      control = list(adapt_delta = 0.8))

# extract samples 
post_pois <- extract.samples(stanfit_pois) 

# save post_pois
save(post_pois, file = "post_pois.RData")

#-----------------------------------------------------------------------------------------------------------------------

# Fit negative binomial (gamma) (to account for over-dispersion)

stanfit_negbin <- cstan(file = "m_negbin.stan", 
                        data = data, 
                        chains = 4, 
                        cores = 60, 
                        control = list(adapt_delta = 0.8))

# extract samples
post_negbin <- extract.samples(stanfit_negbin) 

# save post_negbin
save(post_negbin, file = "post_negbin.RData")

#-----------------------------------------------------------------------------------------------------------------------

# cohort model

cohorts <- seq(from = 1860, to = 1910, by = 1)

d_split <- split(d, d$b_y)
d_split <- d_split[as.character(cohorts)]

m_output_list <- mclapply(d_split, function(dm) {

  person_ids <- sort(unique(dm$person_id))
  n_rp <- length(person_ids)
  rp_sub <- sample(person_ids, size = n_rp)
  
  dm$person_id <- match(dm$person_id, person_ids)

  age_list <- sort(unique(dm$age))
  
  dm$age_bin <- dm$age + 1
  d_mat <- as.matrix(dist(age_list, upper = TRUE, diag = TRUE))
  
  # data list
  data <- list(N_ages = length(age_list),
               N_ind = length(person_ids),
               N = nrow(dm), 
               y = dm$n_moves,
               age = dm$age_bin,
               person_id = dm$person_id,
               d_mat = d_mat,
               run_estimation = 1)
  
  # run model
  m_pois_coh <- cstan(file = "m_pois.stan", 
                      data = data, 
                      chains = 1, 
                      cores = 1, 
                      control = list(adapt_delta = 0.8))

  return(m_pois_coh)

}, mc.cores = 3)
# where mc.cores times stan chains will occupy that many cores

m_coh_samples <- lapply(m_output_list, extract.samples)
