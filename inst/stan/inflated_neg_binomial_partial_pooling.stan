functions {
  int num_zeros(array[] int y) {
    int sum = 0;
    for (n in 1:size(y)) {
      sum += (y[n] == 0);
    }
    return sum;
  }
}

data {
  int < lower = 0 > N; // num of respondents
  int < lower = 0 > K; // num of subpopulations
  int y[N,K]; // ARD
  int < lower = 1 > J; // num groups
  array [N] int < lower = 1, upper = J > jj;  // group for individual i
}

transformed data {
  int < lower = 0 > N_zero[K];
  int < lower = 0 > N_nonzero[K];
  for(k in 1:K){
    N_zero[k] = num_zeros(y[,k]);
    N_nonzero[k] = N - N_zero[k];
  }

  int total_nonzero = sum(N_nonzero);
  int total_zero = (N*K) - total_nonzero;

  int y_nonzero[total_nonzero];
  int index_nz[total_nonzero];
  int index_z[total_zero];
  {
    int index1 = 1;
    int index2 = 1;  
    for(k in 1:K){
      for(n in 1:N){
        if (y[n, k] == 0){
          index_z[index1] = n*k;
          index1 += 1;
        }else{
          y_nonzero[index2] = y[n, k];
          index_nz[index2] = n*k;
          index2 += 1;
        }
      }
    }
  }
}

parameters {
  row_vector < lower = 0 > [K] rho_j [J];
  row_vector < lower = 0 > [K] mu_rho;
  row_vector < lower = 0 > [K] sigma_rho;
  
  real < lower = 0, upper = 1.5 > sigma_delta;
  vector < lower = 0 > [N] delta;

  row_vector < lower = 0 > [K] w;
  array[K] real < lower = 0 , upper = 1 > theta;
}

transformed parameters {
  vector[sum(N_zero)] lambda_zero;
  vector[sum(N_nonzero)] lambda_nonzero;

  row_vector < lower = 0 > [K] shape_rho;
  row_vector < lower = 0 > [K] rate_rho;
  matrix[N,K] lambda;
  {

    shape_rho = square(mu_rho) ./ square(sigma_rho);
    rate_rho = mu_rho ./ square(sigma_rho);

    for(n in 1:N){
      lambda[n] = exp(delta[n] - rho_j[jj[n]]);
    }

  vector[N*K] flattened_lambda = to_vector(lambda);
  lambda_nonzero = flattened_lambda[index_nz];
  lambda_zero = flattened_lambda[index_z];  
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

  {
    int index1 = 1;
    int index2 = 1;

    for(k in 1:K){
      target += N_zero[k] 
                * log_sum_exp(log(theta[k]), 
                              log1m(theta[k]) 
                              + neg_binomial_2_lpmf( 0 | segment(lambda_zero, index1, N_zero[k]), w[k]));
      target += N_nonzero[k] * log1m(theta[k]);
      target += neg_binomial_2_lpmf( segment(y_nonzero, index2, N_nonzero[k])  
                                    | segment(lambda_nonzero, index2, N_nonzero[k]), w[k]);
      index1 += N_zero[k];
      index2 += N_nonzero[k];
    }
  }
}
