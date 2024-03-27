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
  vector < lower = 0 > [K] shape_rho;
  vector < lower = 0 > [K] scale_rho;
  
  vector < lower = 0 > [N] delta;

  cholesky_factor_corr[K] L;

  vector < lower = 0 > [K] tau_n;

  matrix[N,K] epsilon;

  vector < lower = 0 > [K] w;
}

transformed parameters {
  matrix[N,K] lambda;

  {
    vector[K] mu;
    vector[K] tau;
    
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
  // shape_rho ~ weibull(8.5, 5.5);
  // rate_rho ~ weibull(3.5, 1);

  shape_rho ~ cauchy(0, 2.5);
  scale_rho ~ cauchy(0, 2.5);

  for(j in 1:J){
    rho_j[j] ~ weibull(shape_rho, scale_rho);
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
