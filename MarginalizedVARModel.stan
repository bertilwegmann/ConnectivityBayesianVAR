data {
  int<lower=0> p; // number of brain regions
  int<lower=0> q; // number of covariates
  int<lower=0> S; // number of subjects
  int<lower=0> qp; // number of VAR coefficients, qp = q*p
  int<lower=0> T; // number of time points
  
  array[S] matrix[p, p] R_s; // array with matrices R_s for all subjects
  array[S] matrix[q, p] E_s; // array with matrices E_s for all subjects
  array[S] matrix[q, q] Q_s_inv; // array with matrices Q_s_inv for all subjects
  
  // Prior settings
  vector[qp] B_0_spec; // Prior Mean Special with ONLY zeros
  cov_matrix[q] Chol_Cov_B; // cholesky decomposition of prior covariance matrix for B 
  int nu_0; // degrees of freedom in Inverse-Wishart prior for Sigma
  cov_matrix[p] Psi_0; // Scale matrix in the prior covariance matrix Sigma
  cov_matrix[qp] I_Mat; // Identity matrix to form the prior covariance matrix for B, ie diag(qp)
}
parameters {
  cov_matrix[p] Sigma; // covariance matrix
  matrix[q, p] B_spec; // matrix of VAR coefficients
  real<lower=p + 2, upper=10000> nu; // Degree of freedom in inverse-Wishart matrix for Sigma_s
}
transformed parameters {
  matrix[q, p] B;
  B = Chol_Cov_B * B_spec * cholesky_decompose(Sigma); // B_0 = 0 
}
model {
  real Sum_logdet;
  matrix[p, p] Part_s;
  
  // priors
  Sigma ~ inv_wishart(nu_0, Psi_0); // prior for the covariance matrix Sigma
  to_vector(B_spec) ~ multi_normal(B_0_spec, I_Mat); // the special prior for vec(B_spec)
  
  Sum_logdet = 0;
  // log-likelihood
  for (s in 1 : S) {
    Part_s = nu * Sigma + R_s[s] + quad_form(Q_s_inv[s], B - E_s[s]);
    Sum_logdet = Sum_logdet + log_determinant(Part_s);
  }
  
  target += S * (lmgamma(p, 0.5 * (T + nu)) - lmgamma(p, 0.5 * nu))
            + 0.5 * S * nu * log_determinant(nu * Sigma)
            - 0.5 * (nu + T) * Sum_logdet; // add log-likelihood term to the posterior
}
