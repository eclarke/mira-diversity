data {
  int<lower=1> n_specimens;   // Number of timepoints
  int<lower=1> n_subjects;    // Number of subjects
  int<lower=1> n_abx;         // Number of antibiotics
  int<lower=0> n_continuous;  // Days of continuous administration
  int subjects[n_specimens];  // Index of subject for each specimen
  real<lower=0> D1[n_specimens];
  real<lower=0> D1_lag[n_specimens];

  // abx[i,j,k]: for each specimen i, if abx j was administered for k continuous days
  int<lower=0, upper=1> abx[n_specimens, n_abx, n_continuous];  
}

parameters {
  real<lower=0> sigma;
  real a;     // Global intercept
  real b_lag; // Lag term

  real a_subj[n_subjects];  // Subject intercepts
  real<lower=0> s_subj;     // Pooled subject variance

  // Non-centered parameterization for b_abx2
  // b_abx2 = za_abx2 + zb_abx2 * s_abx2
  real zb_abx2[n_subjects, n_abx, n_continuous];
  real za_abx2[n_abx, n_continuous];
  real<lower=0> s_abx2[n_abx, n_continuous];
}
transformed parameters {
  real mu[n_specimens];

  real b_abx2[n_subjects, n_abx, n_continuous];

  for (i in 1:n_specimens) {
    int subject_id = subjects[i];
    mu[i] = a + a_subj[subject_id] + b_lag * D1_lag[i];
    // Add each antibiotic-specific coefficient separately
    for (j in 1:n_abx) {
      for (k in 1:n_continuous) {
        b_abx2[subject_id, j, k] = za_abx2[j, k] + zb_abx2[subject_id, j, k] * s_abx2[j, k];
        mu[i] += (b_abx2[subject_id, j, k]) * abx[i, j, k];
      }
    }
  }
}
model {
  a ~ normal(0, 1);
  b_lag ~ normal(0, 5);
  
  a_subj ~ normal(0, s_subj);
  s_subj ~ normal(0, 5);
  
  sigma ~ normal(0, 5);

  for (k in 1:n_continuous) {
    for (j in 1:n_abx) {
      za_abx2[j, k] ~ normal(0, 1);
      s_abx2[j, k] ~ normal(0.5, 0.5);
      for (i in 1:n_subjects) {
        zb_abx2[i, j, k] ~ normal(0, 1);
      }
    }
  }
  
  D1 ~ lognormal(mu, sigma);
}

