functions {
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
  int sex[N];                               // NEW input
  int person_id[N];
  matrix[N_ages, N_ages] d_mat;
  int<lower=0, upper=1> run_estimation;
}
parameters {
  array[2] vector[N_ages] z;
  vector[N_ind] z_id;
  real<lower=0> sd_id;
  array[2] real<lower=0> eta;               // eta for each sex
  array[2] real<lower=0> rho;               // rho for each sex
  real mu;
}
transformed parameters{
  array[2] vector[N_ages] beta;             // beta vector for each sex
  vector[N] lambda;
  vector[N_ind] a;
  {
    matrix[N_ages, N_ages] L_SIGMA;
    matrix[N_ages, N_ages] SIGMA;
    // Calculate the covariance for the GP
    for ( s in 1:2 ) {
      SIGMA = cov_GPL2(d_mat, eta[s], rho[s], 0.01);
      // cholesky factor of a covariance
      L_SIGMA = cholesky_decompose(SIGMA);
      // covariance matrix = Cholesky covariance factor * z_scores
      beta[s] = L_SIGMA * z[s];
    }
    a = sd_id * z_id;
    // Compute lambda
    for (i in 1:N) {
                                         // access correct beta vector
      lambda[i] = mu + a[person_id[i]] + beta[sex[i],age[i]];
      //lambda[i] = mu + beta[age[i]];
    }
  }
}
model {
  rho ~ normal(3, 3);
  eta ~ exponential(1);
  mu ~ normal(0, 0.05);
  for(s in 1:2) z[s] ~ normal(0, 1);    // subset z prior
  z_id ~ normal(0, 1);
  sd_id ~ normal(0, 1);
  if (run_estimation == 1){
    y ~ poisson_log(lambda);
  }
}
generated quantities {
  int<lower = 0> y_sim[N] = poisson_log_rng(lambda);
}
