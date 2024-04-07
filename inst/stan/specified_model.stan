data {
  int < lower = 0 > N; // num of respondents
  int < lower = 0 > K; // num of subpopulations
  int y[N,K]; // ARD
  int < lower = 1 > J; // num groups
  array [N] int < lower = 1, upper = J > jj;  // group for individual i
  vector < lower = 0 > [K] shape_rho;
  vector < lower = 0 > [K] rate_rho;
  row_vector < lower = 0 > [K] bias [N];

}

parameters {
  row_vector < lower = 0 > [K] rho_j [J];
  row_vector < lower = 0 > [K] w;
}

transformed parameters {
  row_vector [K] lambda [N];
  {
    for(n in 1:N){
      lambda[n] = exp(rho_j[jj[n]] .* bias[n]);
    }
  }
}

model {
  for(j in 1:J){
    rho_j[j] ~ gamma(shape_rho, rate_rho);
  } 

  w ~ cauchy(0, 2.5);

  for(n in 1:N){
    y[n] ~ neg_binomial_2(lambda[n], w);
  }
}
