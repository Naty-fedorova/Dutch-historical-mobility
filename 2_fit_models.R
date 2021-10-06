# dependencies
library(rethinking)
library(rstan)
library(viridis)
library(cmdstanr)
library(parallel)

# load in data for model
# when working with real data:
d <- read.csv("s_person_year_df.csv", stringsAsFactors = FALSE)

# when working with simulated data
d <- read.csv("s_person_year_sim.csv", stringsAsFactors = FALSE)

# select subset
set.seed(1)
person_ids <- sort(unique(d$person_id))
# if selecting subset of rps
n_rp <- 1000

# if wanting to work with entire set
n_rp <- length(person_ids) 

rp_sub <- sample(person_ids, size = n_rp)
dm <- subset(d, d$person_id %in% rp_sub)
person_ids <- sort(unique(dm$person_id))
dm$person_id <- match(dm$person_id, person_ids)

age_list <- sort(unique(dm$age))

# re-index age for model
dm$age_bin <- dm$age + 1

# distance matrix for ages
d_mat <- as.matrix(dist(age_list, upper = TRUE, diag = TRUE))

# data list
data <- list(N_ages = length(age_list),
             N_ind = length(person_ids),
             N = nrow(dm), 
             y = dm$n_moves,
             age = dm$age_bin,
             person_id = dm$person_id,
             d_mat = d_mat)

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
cohorts <- c(1860, 1870, 1900)

d_split <- split(d, d$b_y)
d_split <- d_split[as.character(cohorts)]

m_pois_rstan <- stan_model("stan/m_pois.stan")

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
               d_mat = d_mat)
  
  
  # run model
  m_pois_coh <- cstan(file = "m_pois.stan", 
                      data = data, 
                      chains = 4, 
                      cores = 60, 
                      control = list(adapt_delta = 0.8))

  return(m_pois_coh)

}, num_cores = 3)

m_coh_samples <- lapply(m_output_list, extract.samples)

