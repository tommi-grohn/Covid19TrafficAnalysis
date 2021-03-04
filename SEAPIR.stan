functions {
  
  # t: current time, y: current state value, theta: c_1 and c_2, x_r: traffic, times, constants, x_i: int constants
  real[] seapir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
    
      int N = x_i[1];
      int n_training = x_i[2];

      real tc[n_training] = x_r[(n_training + 1):(2*n_training)];

      real D_e = x_r[2*n_training + 1];
      real D_p = x_r[2*n_training + 2];
      real D_i = x_r[2*n_training + 3];
      real r   = x_r[2*n_training + 4];
      real r_a = x_r[2*n_training + 5];
      real r_p = x_r[2*n_training + 6];

      real S = y[1];
      real E = y[2];
      real A = y[3];
      real P = y[4];
      real I = y[5];
      real R = y[6];
      
      real dS_dt;
      real dE_dt;
      real dA_dt;
      real dP_dt;
      real dI_dt;
      real dR_dt;

      # Beta consists of c_1 term and c_2 * tau
      real beta = theta[1];
      
      # Step function
      int i = 1;
      while (tc[i] < t)
        i = i + 1; 

      beta = beta + x_r[i] * theta[2];
      
      dS_dt = - (r_a*beta*S*A/N + r_p*beta*S*P/N + beta*S*I/N);
      dE_dt =  r_a*beta*S*A/N + r_p*beta*S*P/N + beta*S*I/N - E/D_e;
      dA_dt =  r*E/D_e - A/D_i;
      dP_dt =  (1-r)*E/D_e - P/D_p;
      dI_dt =  P/D_p - I/D_i;
      dR_dt =  A/D_i + I/D_i;
      
      return {dS_dt, dE_dt, dA_dt, dP_dt, dI_dt, dR_dt};
  }
  
}

data {

  real prior_means[2];
  real prior_stds[2];
    
  int<lower=1>  n_training;
  int<lower=0>  n_test;
  real<lower=0> y0[6];  // 6 stages
  real<lower=0> t0;
  real<lower=0> t_training[n_training];
  real<lower=0> t_test[n_test];
  int<lower=1>  N;
  
  real<lower=0> D_e;   // average exposure time
  real<lower=0> D_p;   // average presymptomatic time
  real<lower=0> D_i;   // average infected time
  real<lower=0> r;     // proportion of asymptotic cases out of all infected
  real<lower=0> r_a;   // coefficient of asymptomatics
  real<lower=0> r_p;   // coefficient of presymptomatics
  
  int<lower=0> infections[n_training];
  int<lower=0> infections_pred[n_test];
  real traffic[n_training];  # this matrix should be changed to vector
  real traffic_pred[n_test];  # this matrix should be changed to vector

}

transformed data {
  
  int n_sum = n_training + n_test;
  int x_i[2] = { N, n_training }; 
  
  real x_r[2*n_training+6];
  real x_r_test[2*n_sum+6]; # time interval, length: n_training + n_test
  int x_i_test[2] = { N, n_sum }; 


  # These contain traffic component, time interval, 6 constants
  # x_r is used only for training
  for (i in 1:n_training) {
      x_r[i] = traffic[i];
  }
  x_r[(n_training + 1):(2*n_training)] = t_training;
  x_r[(2*n_training + 1):(2*n_training + 6)] = {D_e, D_p, D_i, r, r_a, r_p};
  
  # This structure is created for lfo because n_training is used for different values
  for (i in 1:n_sum) {
    if (i < n_training + 1) {
      x_r_test[i] = traffic[i];
    }
    else {
      x_r_test[i] = traffic_pred[i-n_training];
      }
  }
  x_r_test[(n_sum + 1):(2*n_sum)] = append_array(t_training, t_test);
  x_r_test[(2*n_sum + 1):(2*n_sum + 6)] = {D_e, D_p, D_i, r, r_a, r_p};

}

parameters {

  real <lower=0> constant;
  real traffic_slope;
  
  real<lower=0> phi_inv; 
  
}

transformed parameters{

  real phi = 1. / phi_inv; 

  real traffic_coeff[2] = append_array({constant}, {traffic_slope});
  
  real<lower=0> y[n_training, 6];
  vector<lower=0> [n_training] lambda;

  y = integrate_ode_rk45(seapir, y0, t0, t_training, traffic_coeff, x_r, x_i);
  lambda = to_vector(y[,4]) / D_p;

}

model {
  
  // priors
  constant ~ normal(prior_means[1], prior_stds[1]); 
  traffic_slope ~ normal(prior_means[2], prior_stds[2]); 
  
  phi_inv ~ exponential(5);
  
  // sampling distribution
  infections ~  neg_binomial_2(lambda, phi); 
}

generated quantities {

  real y_hat[n_sum, 6];
  vector [n_sum] lambda_hat;
  int <lower=0> infections_hat[n_sum];
  
  real log_lik_training;
  real log_lik_pred;

  y_hat = integrate_ode_rk45(seapir, y0, t0, append_array(t_training, t_test), traffic_coeff, x_r_test, x_i_test);
  for (i in 1:n_sum) {
    lambda_hat[i] = max([1e-9, y_hat[i,4] / D_p]); 
  }
  infections_hat = neg_binomial_2_rng(lambda_hat, phi);
  
  log_lik_training = neg_binomial_2_lpmf(infections| lambda, phi);
  log_lik_pred = neg_binomial_2_lpmf(infections_pred| lambda_hat[n_training + 1:], phi);

}
