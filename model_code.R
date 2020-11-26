# model code

# dependencies
library(rethinking)
library(rstan)
library(viridis)

# create folder to contain plots
dir.create("Figures")

# load in data for model
# when working with real data:
# d <- read.csv("s_person_year_df.csv", stringsAsFactors = FALSE)

# when working with simulated data
d <- read.csv("s_person_year_sim.csv", stringsAsFactors = FALSE)

# select subset
set.seed(1)
person_ids <- sort(unique(d$person_id))
n_rp <- 100
# n_rp <- length(person_ids)    # if wanting to work with entire set
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

# stan code for non-centered poisson model
m_pois_gaus_nc_full <-
  "functions {
    matrix cov_GPL2(matrix x, real eta, real rho, real delta) {
      int N = dims(x)[1];
      matrix[N, N] K;
      for (i in 1:(N-1)) {
        K[i, i] = eta + delta;
        for (j in (i + 1):N) {
          K[i, j] = eta * exp(-rho * square(x[i,j]) );
          K[j, i] = K[i, j];
        }
      }
      K[N, N] = eta + delta;
      return K;
    }
  }
data {
  int N_ages;
  int N_ind;
  int N;
  int y[N];
  int age[N];
  int person_id[N];
  matrix[N_ages, N_ages] d_mat;
}
parameters {
  vector[N_ages] z;
  vector[N_ind] z_id;
  real mu;
  real<lower=0> sd_id;
  real<lower=0> eta;
  real<lower=0> rho;
}
transformed parameters{
  vector[N_ages] beta;
  vector[N] lambda;
  matrix[N_ages, N_ages] L_SIGMA;
  matrix[N_ages, N_ages] SIGMA;
  vector[N_ind] a = sd_id * z_id;
  // Calculate the covariance for the GP
  SIGMA = cov_GPL2(d_mat, eta, rho, 0.01);
  // cholesky factor of a covariance
  L_SIGMA = cholesky_decompose(SIGMA);
  // covariance matrix = Cholesky covariance factor * z_scores
  beta = L_SIGMA * z;
    // Compute lambda
    for (i in 1:N) {
    lambda[i] = mu + a[person_id[i]] + beta[age[i]];
    }
}
model {
  rho ~ gamma(2, 2);
  eta ~ normal(0, 1);
  mu ~ normal(0, 1);
  z ~ normal(0, 1);
  z_id ~ normal(0, 1);
  sd_id ~ normal(0, 1);
  y ~ poisson_log(lambda);
}"

# compile
m_nc_real_bin <- stan_model(model_code = m_pois_gaus_nc_full)

# run model
m_nc_real <- sampling(m_nc_real_bin,
                      data = data,
                      iter = 2e3, 
                      chains = 1, 
                      cores = 1, 
                      control = list(adapt_delta = 0.8))


# extract samples and save
post <- extract.samples(m_nc_real) 

# if wanting to save samples, uncomment line below
# save(post, file = "samples_m_nc_real.RData")

# plot for reporting model estimates
png("Figures/model_estimates.png", res = 300, height = 10, width = 8, units = "cm")
plot(precis(m_nc_real, pars = c("mu", "rho", "eta"), depth = 2), ylab = "Parameter")
dev.off()

# lambda
lambda_mod <- apply(post$lambda, 2, mean)
exp(mean(lambda_mod))
HPDI(exp(lambda_mod))

# beta estimates
plot(precis(m_nc_real, depth = 2, pars = "beta"))

#-----------------------------------------------------------------------------------------------------------------------

# stan code for nc model with negative binomial (gamma) (to account for over-dispersion)
m_negbin_gaus_nc_full <-
  "functions {
    matrix cov_GPL2(matrix x, real eta, real rho, real delta) {
      int N = dims(x)[1];
      matrix[N, N] K;
      for (i in 1:(N-1)) {
        K[i, i] = eta + delta;
        for (j in (i + 1):N) {
          K[i, j] = eta * exp(-rho * square(x[i,j]) );
          K[j, i] = K[i, j];
        }
      }
      K[N, N] = eta + delta;
      return K;
    }
  }
data {
  int N_ages;
  int N_ind;
  int N;
  int y[N];
  int age[N];
  int person_id[N];
  matrix[N_ages, N_ages] d_mat;
}
parameters {
  vector[N_ages] z;
  vector[N_ind] z_id;
  real mu;
  real<lower=0> sd_id;
  real<lower=0> eta;
  real<lower=0> rho;
  real<lower=0> phi;
}
transformed parameters{
  vector[N_ages] beta;
  matrix[N_ages, N_ages] L_SIGMA;
  matrix[N_ages, N_ages] SIGMA;
  vector[N_ind] a = sd_id * z_id;
  // Calculate the covariance for the GP
  SIGMA = cov_GPL2(d_mat, eta, rho, 0.01);
  // cholesky factor of a covariance
  L_SIGMA = cholesky_decompose(SIGMA);
  // covariance matrix = Cholesky covariance factor * z_scores
  beta = L_SIGMA * z;
    // Compute lambda
  for (i in 1:N) {
    lambda[i] = exp(mu + a[person_id[i]] + beta[age[i]]);
  }
}
model {
  vector[N] lambda;
  rho ~ gamma(2, 2);
  eta ~ normal(0, 1);
  mu ~ normal(0, 1);
  z ~ normal(0, 1);
  z_id ~ normal(0, 1);
  sd_id ~ normal(0, 1);
  phi ~ exponential(1);
  y ~ neg_binomial_2(lambda, phi);
}"

# compile
m_nc_negbin_s <- stan_model(model_code = m_negbin_gaus_nc_full)

# run model
m_nc_negbin <- sampling(m_nc_negbin_s,
                      data = data,
                      iter = 2e3, 
                      chains = 1, 
                      cores = 1, 
                      control = list(adapt_delta = 0.8))


# extract samples and save
post_g <- extract.samples(m_nc_negbin) 

# plots for reporting gamma model
png("Figures/negbin_estimates.png", res = 300, height = 10, width = 8, units = "cm")
plot(precis(m_nc_negbin, pars = c("mu", "rho", "eta", "phi"), depth = 2), ylab = "Parameter")    # TODO add lambda to this list - model code should produce it....
dev.off()

# define colors for plotting
vir <- viridis(20)
vir_int <- viridis(20, alpha = 0.3)

# betas
beta_mod <- apply(post$beta, 2, mean)
beta_mod_int <- apply(post$beta, 2, HPDI)
mu_mod <- mean(post$mu)

png("Figures/beta_estimates_negbin.png", res = 300, height = 15, width = 25, units = "cm")
par(mar = c(5.1,4.1,4.1,0.5))
plot(y = exp(beta_mod + mu_mod), 
     x = 0:(length(beta_mod)-1), 
     ylim = c(0, 0.6),
     xlim = c(0,(length(beta_mod)-1)),
     xlab = "Age", 
     pch = 19, 
     ylab = "Estimated number of moves per year", 
     col = vir[7], 
     main = "Gamma-poisson model estimate of number of moves per year per age", 
     font.main = 1,
     bty = "n")
shade(exp(beta_mod_int + mu_mod), 0:(length(beta_mod)-1), col = vir_int[1])
dev.off()

# alphas
a_mod <- apply(post$a, 2, mean)
a_mod_int <- apply(post$a, 2, HPDI)

a <- t(rbind(a_mod_int, a_mod))
a_o <- a[order(a[,"a_mod"]), ]
a_int <- t(a_o[,1:2])

png("Figures/a_est_negbin.png", res = 300, height = 15, width = 20, units = "cm")
plot(exp(a_o[,"a_mod"] + mu_mod), 
     xlab = "RP index", 
     ylab = "Moves per year", 
     main = "Individual variation in moves per year", 
     font.main = 1,
     col = vir[7],
     bty = "n", 
     ylim = c(0,2))
shade(exp(a_int + mu_mod), 1:length(a_mod), col = vir_int[1])
dev.off()


#-----------------------------------------------------------------------------------------------------------------------

# cohort model
cohorts <- seq(from = 1860, to = 1910, by = 1)
m_output_list <- as.list(cohorts)
m_coh_samples <- as.list(cohorts)

for(i in 1:length(cohorts)){
  
  b_seq <- seq(from = cohorts[i], to = cohorts[i], by = 1)
  d_coh <- d[which(d$b_y %in% b_seq), ]
  
  set.seed(1)
  person_ids <- sort(unique(d_coh$person_id))
  n_rp <- 100
  rp_sub <- sample(person_ids, size = n_rp)
  
  dm <- subset(d, d$person_id %in% rp_sub)
  
  person_ids <- sort(unique(dm$person_id))
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
  m_nc_coh <- sampling(m_nc_real_bin,
                       data = data,
                       iter = 2e3, 
                       chains = 1, 
                       cores = 1, 
                       control = list(adapt_delta = 0.8))
  
  
  m_output_list[i] <- list(m_nc_coh)
  m_coh_samples[i] <- list(extract.samples(m_nc_coh))
  
}

