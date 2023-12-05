data {
  int < lower = 0 > N; // num of respondents
  int < lower = 0 > K; // num of subpopulations
  int y[N,K]; // ARD
  vector[N] offset; // municipality population of respondent i
  int < lower = 1 > J; // num groups
  array[N] int < lower = 1, upper = J > jj;  // group for individual i
}

parameters {
  array[J] row_vector[K] rho_j;
  vector[K] mu_rho;
  vector < lower = 0 > [K] sigma_rho;

  vector[N] delta;
  real < lower = 0 > sigma_delta;

  cholesky_factor_corr[K] L;
  vector[K] tau_n;

  matrix[N,K] epsilon;

  vector[K] reciprocal_w;
}

transformed parameters {
  matrix[N,K] lambda;
  vector < lower = 0 > [K] w;

  {
    for(k in 1:K){
      w[k] = 1.0 / sqrt(abs(reciprocal_w[k]));
    }

    vector[K] mu;
    vector[K] tau;
    
    mu = log(1.0 ./ sqrt(1.0 + square(tau_n)));
    tau = sqrt(log(1.0 + square(tau_n)));

    matrix[N,K] rho;
    for(n in 1:N){
      rho[n] = rho_j[jj[n]];
    }
    
    lambda = exp(rho + rep_matrix(sigma_delta * delta, K) + 
                  rep_matrix(mu, N)' + (diag_pre_multiply(tau, L) * epsilon')') .* rep_matrix(offset, K);
  }
}

model {
  mu_rho ~ normal(0, 10);
  sigma_rho ~ cauchy(0, 2.5);

  for(j in 1:J){
    rho_j[j] ~ normal(mu_rho, sigma_rho);
  } 

  sigma_delta ~ cauchy(0, 2.5);
  delta ~ std_normal();
  
  tau_n ~ cauchy(0, 2.5);

  to_vector(epsilon) ~ std_normal();

  reciprocal_w ~ std_normal();

  L ~ lkj_corr_cholesky(2);
  
  for(k in 1:K){
    y[,k] ~ neg_binomial_2(lambda[,k], w[k]);
  }
}
