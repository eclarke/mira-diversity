/* 
* Autoregressive binomial model, v3.0
* Changes from v2.x.x
*  - Multiple abx lag terms corresponding to today (AR(0)), yesterday (AR(1)), AR(2), etc
* 
* D1_delta ~ LogNormal(mu, sigma)
* mu[i,j,k] = a + a_spec[i] + a_subj[j] * ((b_abx1[k] * abx[i,k]) + (b_abx2[j,k] * abx[i,k])) * (1-prev[i]) + (b_lag * prev[i])
* 
* a = global intercept
* a_spec[i] = specimen i intercept (partially pooled between specimens)
* a_subj[j] = subject j intercept (partially pooled btwn subjects)
* b_abx1[k] = antibiotic k slope (partially pooled btwn abx)
* abx[i,k] = indicator variable for antibiotic k at specimen i
* b_abx2[j,k] = subject j + antibiotic k slope (partially pooled between subjects)
* b_lag = coefficient for previous specimen's reads
* prev[i] = previous specimen's reads (as proportion of total)
*
*/
data {
  int<lower=1> n_specimens;   // Number of timepoints
  int<lower=1> n_subjects;    // Number of subjects
  int<lower=1> n_abx;         // Number of antibiotics
  int<lower=0> n_lag;         // Days of lag 
  int subjects[n_specimens];  // Index of subject for each specimen
  real<lower=0> D1_delta[n_specimens];

  int<lower=0, upper=1> abx[n_specimens, n_abx, n_lag];  // Antibiotics on/off for each specimen/abx at each lag 
}

parameters {
  real<lower=0> sigma;
  real a;     // Global intercept
  real b_lag; // Global lag coefficient

  real a_spec[n_specimens]; // Specimen intercepts
  real<lower=0> s_spec;     // Pooled specimen variance
  
  real a_subj[n_subjects];  // Subject intercepts
  real<lower=0> s_subj;     // Pooled subject variance

  // Non-centered parameterization for b_abx2
  // b_abx2 = za_abx2 + zb_abx2 * s_abx2
  real zb_abx2[n_subjects, n_abx, n_lag];
  real za_abx2[n_abx, n_lag];
  real<lower=0> s_abx2[n_abx, n_lag];
}
transformed parameters {
  real mu[n_specimens];

  real b_abx2[n_subjects, n_abx, n_lag];

  for (i in 1:n_specimens) {
    int subject_id = subjects[i];
    mu[i] = a + a_spec[i] + a_subj[subject_id];
    // Add each antibiotic-specific coefficient separately
    for (j in 1:n_abx) {
      for (k in 1:n_lag) {
        b_abx2[subject_id, j, k] = za_abx2[j, k] + zb_abx2[subject_id, j, k] * s_abx2[j, k];
        mu[i] += (b_abx2[subject_id, j, k]) * abx[i, j, k];
      }
    }
  }
}
model {
  a ~ normal(0, 1);
  b_lag ~ normal(0, 1);
  
  a_subj ~ normal(0, s_subj);
  s_subj ~ normal(0, 5);
  
  a_spec ~ normal(0, s_spec);
  s_spec ~ normal(0, 5);
  
  sigma ~ normal(0, 5);

  for (k in 1:n_lag) {
    for (j in 1:n_abx) {
      za_abx2[j, k] ~ normal(0, 1);
      s_abx2[j, k] ~ normal(0, 1);
      for (i in 1:n_subjects) {
        zb_abx2[i, j, k] ~ normal(0, 1);
      }
    }
  }
  
  D1_delta ~ lognormal(mu, sigma);
}

