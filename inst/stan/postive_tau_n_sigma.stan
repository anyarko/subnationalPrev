data {
  int < lower = 0 > N; // num of respondents
  int < lower = 0 > K; // num of subpopulations
  int y[N,K]; // ARD
  vector[N] offset; // municipality population of respondent i
  int < lower = 1 > J; // num groups
  array[N] int < lower = 1, upper = J > jj;  // group for individual i
}

parameters {
  array[J] row_vector < lower = 0 > [K] rho_j;
  vector < lower = 0 > [K] mu_rho;
  vector < lower = 0 > [K] sigma_rho;
  
  vector < lower = 0 > [N] delta;

  cholesky_factor_corr[K] L;

  vector < lower = 0 > [K] tau_n;

  matrix[N,K] epsilon;

  vector < lower = 0 > [K] w;
}

transformed parameters {
  matrix[N,K] lambda;
  vector < lower = 0 > [K] shape_rho;
  vector < lower = 0 > [K] rate_rho;

  {
    vector[K] mu;
    vector[K] tau;

    shape_rho = square(mu_rho) ./ square(sigma_rho);
    rate_rho = mu_rho ./ square(sigma_rho);
    
    mu = log(1.0 ./ sqrt(1.0 + square(tau_n)));
    tau = sqrt(log(1.0 + square(tau_n)));

    matrix[N,K] rho;
    for(n in 1:N){
      rho[n] = rho_j[jj[n]] * -1;
    }

    lambda = exp(rho - rep_matrix(delta, K) + 
                  rep_matrix(mu, N)' + 
                  (diag_pre_multiply(tau, L) * epsilon')') .* rep_matrix(offset, K);                  
  }
}

model {
  mu_rho ~ cauchy(0, 2.5);
  sigma_rho ~ cauchy(0, 2.5);

  for(j in 1:J){
    rho_j[j] ~ gamma(shape_rho, rate_rho);
  } 

  delta ~ weibull(35, 10);
  
  tau_n ~ cauchy(0, 2.5);

  to_vector(epsilon) ~ std_normal();

  w ~ cauchy(0, 2.5);

  L ~ lkj_corr_cholesky(2);
  
  for(k in 1:K){
    y[,k] ~ neg_binomial_2(lambda[,k], w[k]);
  }
}
