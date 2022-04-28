
# dependencies
library(rethinking)
library(rstan)
library(cmdstanr)
library(parallel)
library(posterior)

# stan code for models
m_pois <- cmdstan_model("m_pois.stan")
m_negbin <- cmdstan_model("m_negbin.stan")

# simulate the person-year table
n_rp <- 36595     # 36,595 in the full analysis dataset

birth_years <- sample(1850:1922, n_rp, replace = TRUE)

genders <- sample(c(1:2), n_rp, prob = c(0.5, 0.5), replace = TRUE)


d <- expand.grid(age = 0:50, person_id = 1:n_rp)
d$b_y <- birth_years[d$person_id]
d$address_start_y <- birth_years[d$person_id] + d$age
d$obs_end <- 60 + birth_years[d$person_id]
d$gender <- genders[d$person_id]
d$n_moves <- (-999) # this is what we want to sim!

person_ids <- sort(unique(d$person_id))

age_list <- sort(unique(d$age))

# re-index age for model; real age 0 (last birthday) is now age_bin[1]
d$age_bin <- d$age + 1

# distance matrix for ages
d_mat <- as.matrix(dist(age_list, upper = TRUE, diag = TRUE))

# data list
data <- list(N_ages = length(age_list),
             N_ind = length(person_ids),
             N = nrow(d), 
             y = d$n_moves, # dummy values, because we want to simulate this
             age = d$age_bin,
             gender = d$gender,  # female = 1, male = 2
             person_id = d$person_id,
             d_mat = d_mat,
             run_estimation = 0)

# simulate from the stan mdoel
d_pois_sim <- m_pois$sample(
  data = data,
  chains = 1,
  iter_sampling = 2,
  fixed_param = TRUE)

draws <- d_pois_sim$draws()
samples <- as.data.frame(as_draws_df(draws))

stansim_pois <- rstan::read_stan_csv(d_pois_sim$output_files())
dat_pois <- extract.samples(stansim_pois)

# outcome has been simulated from our prior:
d$n_moves <- dat_pois$y_sim[1,]

# can save for future use if need be
write.csv(d, "d_sim.csv", row.names = FALSE)

###############

# if you have access to HSN, can otherwise load `s_person_year_df.csv` here
#d <- read.csv("s_person_year_df.csv", stringsAsFactors = FALSE)

# select subset
set.seed(1)
person_ids <- sort(unique(d$person_id))
n_rp <- 36595
# N is 36595 in the full sample

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
d_mat <- as.matrix(dist(age_list, upper = TRUE, diag = TRUE)) # note this takes length of age_list, so if ages are not a complete sequence, it will cause problems with indexing

# normalize distance matrix
d_mat <- (d_mat - min(d_mat))/(max(d_mat) - min(d_mat)) 

# data list
data <- list(N_ages = length(age_list),
             N_ind = length(person_ids),
             N = nrow(dm), 
             y = dm$n_moves,
             age = dm$age_bin,
             gender = dm$gender,  # female = 1, male = 2
             person_id = dm$person_id,
             d_mat = d_mat,
             run_estimation = 1)

# fit poisson
stanfit_pois <- cstan(file = "m_pois.stan", 
                      data = data,
                      seed = 101,
                      chains = 4, 
                      cores = 60, 
                      iter = 1000,
                      warmup = 500,
                      control = list(adapt_delta = 0.8, max_treedepth = 15))

# if needing to save fit
saveRDS(stanfit_pois, file = "stanfit_pois.rds")

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
                        iter = 1000,
                        seed = 101,
                        control = list(adapt_delta = 0.8, max_treedepth = 15)) # recommended increasing adapt_delta

# extract samples
post_negbin <- extract.samples(stanfit_negbin) 

# save post_negbin
# save(post_negbin, file = "post_negbin.RData")

#-----------------------------------------------------------------------------------------------------------------------

# cohort model
num_cores <- 20

cohorts <- seq(from = 1850, to = 1922, by = 1)

d_split <- split(d, d$b_y)
d_split <- d_split[as.character(cohorts)]

m_output_list <- mclapply(d_split, function(dm) {
  
  person_ids <- sort(unique(dm$person_id))
  n_rp <- length(person_ids)
  rp_sub <- sample(person_ids, size = n_rp)
  
  dm$person_id <- match(dm$person_id, person_ids)
  
  dm$age_bin <- dm$age + 1
  
  age_list <- sort(unique(dm$age_bin))
  
  d_mat <- as.matrix(dist(age_list, upper = TRUE, diag = TRUE))
  
  # data list
  data <- list(N_ages = length(age_list),
               N_ind = length(person_ids),
               N = nrow(dm), 
               y = dm$n_moves,
               age = dm$age_bin,
               gender = dm$gender,  # female = 1, male = 2
               person_id = dm$person_id,
               d_mat = d_mat,
               run_estimation = 1)
  
  # run model
  m_pois_coh <- cstan(file = "m_pois.stan", 
                      data = data, 
                      chains = 4, 
                      cores = 4, 
                      control = list(adapt_delta = 0.8))
  
  return(m_pois_coh)
  
}, mc.cores = num_cores)
# where `mc.cores` times stan's `chains` will occupy that many cores

m_coh_samples <- lapply(m_output_list, extract.samples, pars = c("beta", "a", "mu", "sd_id"))
