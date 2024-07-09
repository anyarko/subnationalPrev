data {
  int < lower = 0 > N; // num of respondents
  int < lower = 0 > K; // num of subpopulations
  int y[N,K]; // ARD
  int < lower = 1 > J; // num groups
  array [N] int < lower = 1, upper = J > jj;  // group for individual i
}

parameters {
  row_vector < lower = 0 > [K] rho_j [J];
  row_vector < lower = 0 > [K] mu_rho;
  row_vector < lower = 0 > [K] sigma_rho;
  
  real < lower = 0, upper = 1.5 > sigma_delta;
  vector < lower = 0 > [N] delta;

  row_vector < lower = 0 > [K] w;
}

transformed parameters {
  row_vector < lower = 0 > [K] shape_rho;
  row_vector < lower = 0 > [K] rate_rho;
  row_vector < lower = 0 > [K] lambda [N];
  {

    shape_rho = square(mu_rho) ./ square(sigma_rho);
    rate_rho = mu_rho ./ square(sigma_rho);

    for(n in 1:N){
      lambda[n] = exp(delta[n] - rho_j[jj[n]]);
    }
  }
}

model {
  mu_rho ~ normal(0, 10);
  sigma_rho ~ std_normal();

  for(j in 1:J){
    rho_j[j] ~ gamma(shape_rho, rate_rho);
  } 

  delta ~ normal(5.5, sigma_delta);

  w ~ normal(0, 10);

  for(n in 1:N){
    y[n] ~ neg_binomial_2(lambda[n], w);
  }
}
