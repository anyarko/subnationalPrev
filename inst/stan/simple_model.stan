data {
  int < lower = 0 > N; // num of respondents
  int < lower = 0 > K; // num of subpopulations
  int y[N,K]; // ARD
  int < lower = 1 > J; // num groups
  array [N] int < lower = 1, upper = J > jj;  // group for individual i
}

parameters {
  row_vector < lower = 0 > [K] rho_j [N];
  vector < lower = 0 > [K] mu_rho;
  vector [K] sigma_rho;

  vector < lower = 0 > [N] delta;

  row_vector < lower = 0 > [K] w;
}

transformed parameters {
  vector < lower = 0 > [K] shape_rho;
  vector < lower = 0 > [K] rate_rho;
  row_vector [K] lambda [N];

  {
    shape_rho = square(mu_rho) ./ square(sigma_rho);
    rate_rho = mu_rho ./ square(sigma_rho);
    
    for(n in 1:N){
      lambda[n] = exp(rho_j[jj[n]] .* -1 + rep_vector(delta[n], K)');
    }
  }
}

model {
  mu_rho ~ cauchy(0, 2.5);
  sigma_rho ~ cauchy(0, 2.5);

  for(j in 1:J){
    rho_j[j] ~ gamma(shape_rho, rate_rho);
  } 

  delta ~ weibull(16, 4);

  for(n in 1:N){
    y[n,] ~ neg_binomial_2(lambda[n], w);
  }
}
