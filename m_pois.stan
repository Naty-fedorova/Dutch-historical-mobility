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
  int person_id[N];
  matrix[N_ages, N_ages] d_mat;
  int<lower=0, upper=1> run_estimation;
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
  vector[N_ind] a = sd_id * z_id;
  {
    matrix[N_ages, N_ages] L_SIGMA;
    matrix[N_ages, N_ages] SIGMA;
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
}
model {
  rho ~ gamma(2, 2);
  eta ~ normal(0, 1);
  mu ~ normal(0, 1);
  z ~ normal(0, 1);
  z_id ~ normal(0, 1);
  sd_id ~ normal(0, 1);
  if (run_estimation == 1){
    y ~ poisson_log(lambda);
  }
}
generated quantities {
  int<lower = 0> y_sim[N] = poisson_log_rng(lambda);
}
