functions {

  // Temperature function for ODE forcing.
  real temp_func(real t, real temp_ref, real temp_rise) {
    return temp_ref + (temp_rise * t) / (80 * 24 * 365) + 10 * sin((2 * pi() / 24) * t) + 10 * sin((2 * pi() / (24 * 365)) * t);
  }

  vector temp_func_vec(vector ts, real temp_ref, real temp_rise) {
    return temp_ref + (temp_rise * ts) / (80 * 24 * 365) + 10 * sin((2 * pi() / 24) * ts) + 10 * sin((2 * pi() / (24 * 365)) * ts);
  }

  // Exogenous SOC input function.
  real i_s_func(real t) {
    return 0.001 + 0.0005 * sin((2 * pi() / (24 * 365)) * t);
  }

  // Exogenous DOC input function.
  real i_d_func(real t) {
    return 0.0001 + 0.00005 * sin((2 * pi() / (24 * 365)) * t);
  }

  // Function for enforcing Arrhenius temperature dependency of ODE parameter.
  real arrhenius_temp(real input, real temp, real Ea, real temp_ref) {
    return input * exp(-Ea / 0.008314 * (1 / temp - 1 / temp_ref));
  }

  vector arrhenius_temp_vec(real input, vector temps, real Ea, real temp_ref) {
    return input * exp(-Ea / 0.008314 * (1 ./ temps - 1 ./ temp_ref));
  }

  // Function for enforcing linear temperature dependency of ODE parameter.
  real linear_temp(real input, real temp, real Q, real temp_ref) {
    return input - Q * (temp - temp_ref);
  }

  vector linear_temp_vec(real input, vector temps, real Q, real temp_ref) {
    return input - Q * (temps - temp_ref);
  }

  /*
  AWB-ECA (equilibrium chemistry approximation) variant of AWB model.
  C[1] is soil organic carbon (SOC) density.
  C[2] is dissolved organic carbon (DOC) density.
  C[3] is microbial biomass carbon (MBC) density.
  C[4] is extracellular enzyme carbon (EEC) density.
  */

  vector AWB_ECA_ODE(real t, vector C,
                     real u_Q_ref, // Reference carbon use efficiency.
                     real Q, // Carbon use efficiency linear dependence negative slope.
                     real a_MSA, // AWB MBC-to-SOC transfer fraction.
                     real K_DE, // SOC decomposition K_m.
                     real K_UE, // DOC uptake K_m.
                     real V_DE_ref, // Reference SOC decomposition V_max.
                     real V_UE_ref, // Reference DOC uptake V_max.
                     real Ea_V_DE, // SOC V_max activation energy.
                     real Ea_V_UE, // DOC V_max activation energy.
                     real r_M, // MBC turnover rate.
                     real r_E, // Enzyme production rate.
                     real r_L, // Enzyme loss rate.
                     real temp_ref,
                     real temp_rise) {

    // Initiate exogenous input and forcing variables for future assignment.
    real temp;
    real i_s;
    real i_d;
    real u_Q; // Forced u_Q_ref.
    real V_DE; // Forced V_DE.
    real V_UE; // Forced V_UE.

    // Initiate dependent expressions.
    real F_S; // SOC decomposition
    real F_D; // DOC uptake

    // Assign input and forcing variables to appropriate value at time t.
    temp = temp_func(t, temp_ref, temp_rise); // x_r[1] is temp_ref 283.
    i_s = i_s_func(t);
    i_d = i_d_func(t);

    // Force temperature dependent parameters.
    u_Q = linear_temp(u_Q_ref, temp, Q, temp_ref);
    V_DE = arrhenius_temp(V_DE_ref, temp, Ea_V_DE, temp_ref);
    V_UE = arrhenius_temp(V_UE_ref, temp, Ea_V_UE, temp_ref);

    // Assign dependent expressions.
    F_S = V_DE * C[4] * C[1] / (K_DE + C[4] + C[1]);
    F_D = V_UE * C[3] * C[2] / (K_UE + C[3] + C[2]);

    // Initiate vector for storing derivatives.
    vector[4] dCdt;

    // Compute derivatives.
    dCdt[1] = i_s + a_MSA * r_M * C[3] - F_S;
    dCdt[2] = i_d + (1 - a_MSA) * r_M * C[3] + F_S + r_L * C[4] - F_D;
    dCdt[3] = u_Q * F_D * - (r_M + r_E) * C[3];
    dCdt[4] = r_E * C[3] - r_L * C[4];
    return dCdt;
  }

  // Calculate model output CO2 observations from states x.
  vector calc_AWB_ECA_CO2(real[] ts,
                          array[] vector C,
                          real u_Q_ref, // Reference carbon use efficiency.
                          real Q, // Carbon use efficiency linear dependence negative slope.
                          real K_UE, // DOC uptake K_m.
                          real V_UE_ref, // Reference DOC uptake V_max.
                          real Ea_V_UE, // DOC V_max activation energy.
                          real temp_ref,
                          real temp_rise) {
    
    vector[size(ts)] ts_vec;
    vector[size(ts)] u_Q;
    vector[size(ts)] V_UE;
    vector[size(ts)] temp_vec;
    vector[size(ts)] CO2;
    
    ts_vec = to_vector(ts);
    //print("ts_vec", ts_vec);
    temp_vec = temp_func_vec(ts_vec, temp_ref, temp_rise);
    //print("temp_vec", temp_vec);
    u_Q = linear_temp_vec(u_Q_ref, temp_vec, Q, temp_ref);
    //print("u_Q", u_Q);
    V_UE = arrhenius_temp_vec(V_UE_ref, temp_vec, Ea_V_UE, temp_ref);
    //print("V_UE", V_UE);
    CO2 = (1 - u_Q) .* (V_UE .* C[3,] .* C[2,]) ./ (K_UE + C[3,] + C[2,]);
  
    return CO2;    
  }

  // From https://discourse.mc-stan.org/t/rng-for-truncated-distributions/3122/12.
  real normal_lb_ub_rng(real mu, real sigma, real lb, real ub) {
      real p1 = normal_cdf(lb, mu, sigma);  // cdf with lower bound
      real p2 = normal_cdf(ub, mu, sigma);  // cdf with upper bound
      real u = uniform_rng(p1, p2);
      return (sigma * inv_Phi(u)) + mu;  // inverse cdf 
  }
}

data {
  int<lower=1> state_dim; // Number of state dimensions (4 for AWB).
  int<lower=1> N_t; // Number of observations.
  array[N_t] real<lower=0> ts; // Univariate array of observation time steps.
  array[state_dim+1] vector<lower=0>[N_t] y; // Multidimensional array of state observations and CO2 bounded at 0. y in [state_dim, N_t] shape to facilitate likelihood sampling.
  real<lower=0> temp_ref; // Reference temperature for temperature forcing.
  real<lower=0> temp_rise; // Assumed increase in temperature (°C/K) over next 80 years.
  real<lower=0> prior_scale_factor; // Factor multiplying parameter means to obtain prior standard deviations.
  real<lower=0> obs_error_scale; // Observation noise factor multiplying observations of model output x_hat.
  vector<lower=0>[state_dim] x_hat0; // Initial ODE conditions.
  // [1] is prior mean, [2] is prior lower bound, [3] is prior upper bound.
  array[3] real<lower=0> u_Q_ref_prior_dist_params;
  array[3] real<lower=0> Q_prior_dist_params;
  array[3] real<lower=0> a_MSA_prior_dist_params;
  array[3] real<lower=0> K_DE_prior_dist_params;
  array[3] real<lower=0> K_UE_prior_dist_params;
  array[3] real<lower=0> V_DE_ref_prior_dist_params;
  array[3] real<lower=0> V_UE_ref_prior_dist_params;
  array[3] real<lower=0> Ea_V_DE_prior_dist_params;
  array[3] real<lower=0> Ea_V_UE_prior_dist_params;
  array[3] real<lower=0> r_M_prior_dist_params;
  array[3] real<lower=0> r_E_prior_dist_params;
  array[3] real<lower=0> r_L_prior_dist_params;
}

transformed data {
  real t0 = 0; // Initial time.
}

parameters {
  real<lower = u_Q_ref_prior_dist_params[2], upper = u_Q_ref_prior_dist_params[3]> u_Q_ref; // Reference carbon use efficiency.
  real<lower = Q_prior_dist_params[2], upper = Q_prior_dist_params[3]> Q; // Carbon use efficiency linear dependence negative slope.
  real<lower = a_MSA_prior_dist_params[2], upper = a_MSA_prior_dist_params[3]> a_MSA; // AWB MBC-to-SOC transfer fraction.
  real<lower = K_DE_prior_dist_params[2], upper = K_DE_prior_dist_params[3]> K_DE; // SOC decomposition K_m.
  real<lower = K_UE_prior_dist_params[2], upper = K_UE_prior_dist_params[3]> K_UE; // DOC uptake K_m.
  real<lower = V_DE_ref_prior_dist_params[2], upper = V_DE_ref_prior_dist_params[3]> V_DE_ref; // Reference SOC decomposition V_max.
  real<lower = V_UE_ref_prior_dist_params[2], upper = V_UE_ref_prior_dist_params[3]> V_UE_ref; // Reference DOC uptake V_max.
  real<lower = Ea_V_DE_prior_dist_params[2], upper = Ea_V_DE_prior_dist_params[3]> Ea_V_DE; // SOC V_max activation energy.
  real<lower = Ea_V_UE_prior_dist_params[2], upper = Ea_V_UE_prior_dist_params[3]> Ea_V_UE; // DOC V_max activation energy.
  real<lower = r_M_prior_dist_params[2], upper = r_M_prior_dist_params[3]> r_M; // MBC turnover rate.
  real<lower = r_E_prior_dist_params[2], upper = r_E_prior_dist_params[3]> r_E; // Enzyme production rate.
  real<lower = r_L_prior_dist_params[2], upper = r_L_prior_dist_params[3]> r_L; // Enzyme loss rate.
}

transformed parameters {
  vector<lower=0>[N_t] x_hat_CO2;
  array[state_dim+1] vector<lower=0>[N_t] x_hat_add_CO2;

  // Solve ODE.
  array[N_t] vector<lower=0>[state_dim] x_hat_intmd = ode_rk45(AWB_ECA_ODE, x_hat0, t0, ts, u_Q_ref, Q, a_MSA, K_DE, K_UE, V_DE_ref, V_UE_ref, Ea_V_DE, Ea_V_UE, r_M, r_E, r_L, temp_ref, temp_rise);

  // Transform model output to match observations y in shape, [state_dim, N_t].
  array[state_dim] vector<lower=0>[N_t] x_hat;
  for (i in 1:N_t) {
    for (j in 1:state_dim) {
      x_hat[j, i] = x_hat_intmd[i, j];
    }
  }

  // Compute CO2.
  x_hat_CO2 = calc_AWB_ECA_CO2(ts, x_hat, u_Q_ref, Q, K_UE, V_UE_ref, Ea_V_UE, temp_ref, temp_rise);   

  // Add CO2 vector to x_hat.
  x_hat_add_CO2[1:4,] = x_hat;
  x_hat_add_CO2[5,] = x_hat_CO2;
  
  //print("Leapfrog x: ", x_hat);
  //print("Leapfrog CO2: ", x_hat_CO2);
  //print("Leapfrog x add CO2: ", x_hat_add_CO2);  
}

model {
  u_Q_ref ~ normal(u_Q_ref_prior_dist_params[1], u_Q_ref_prior_dist_params[1] * prior_scale_factor) T[u_Q_ref_prior_dist_params[2], u_Q_ref_prior_dist_params[3]];
  Q ~ normal(Q_prior_dist_params[1], Q_prior_dist_params[1] * prior_scale_factor) T[Q_prior_dist_params[2], Q_prior_dist_params[3]];
  a_MSA ~ normal(a_MSA_prior_dist_params[1], a_MSA_prior_dist_params[1] * prior_scale_factor) T[a_MSA_prior_dist_params[2], a_MSA_prior_dist_params[3]];
  K_DE ~ normal(K_DE_prior_dist_params[1], K_DE_prior_dist_params[1] * prior_scale_factor) T[K_DE_prior_dist_params[2], K_DE_prior_dist_params[3]];
  K_UE ~ normal(K_UE_prior_dist_params[1], K_UE_prior_dist_params[1] * prior_scale_factor) T[K_UE_prior_dist_params[2], K_UE_prior_dist_params[3]];
  V_DE_ref ~ normal(V_DE_ref_prior_dist_params[1], V_DE_ref_prior_dist_params[1] * prior_scale_factor) T[V_DE_ref_prior_dist_params[2], V_DE_ref_prior_dist_params[3]];
  V_UE_ref ~ normal(V_UE_ref_prior_dist_params[1], V_UE_ref_prior_dist_params[1] * prior_scale_factor) T[V_UE_ref_prior_dist_params[2], V_UE_ref_prior_dist_params[3]];
  Ea_V_DE ~ normal(Ea_V_DE_prior_dist_params[1], Ea_V_DE_prior_dist_params[1] * prior_scale_factor) T[Ea_V_DE_prior_dist_params[2], Ea_V_DE_prior_dist_params[3]];
  Ea_V_UE ~ normal(Ea_V_UE_prior_dist_params[1], Ea_V_UE_prior_dist_params[1] * prior_scale_factor) T[Ea_V_UE_prior_dist_params[2], Ea_V_UE_prior_dist_params[3]];
  r_M ~ normal(r_M_prior_dist_params[1], r_M_prior_dist_params[1] * prior_scale_factor) T[r_M_prior_dist_params[2], r_M_prior_dist_params[3]];
  r_E ~ normal(r_E_prior_dist_params[1], r_E_prior_dist_params[1] * prior_scale_factor) T[r_E_prior_dist_params[2], r_E_prior_dist_params[3]];
  r_L ~ normal(r_L_prior_dist_params[1], r_L_prior_dist_params[1] * prior_scale_factor) T[r_L_prior_dist_params[2], r_L_prior_dist_params[3]];
  //print("Leapfrog theta: ", "u_Q_ref = ", u_Q_ref, ", Q = ", Q, ", a_MSA = ", a_MSA, ", K_DE = ", K_DE, ", K_UE = ", K_UE, ", V_DE_ref = ", V_DE_ref, ", V_UE_ref = ", V_UE_ref, ", Ea_V_DE = ", Ea_V_DE, ", Ea_V_UE = ", Ea_V_UE, ", r_M = ", r_M, ", r_E = ", r_E, ", r_L = ", r_L);

  // Likelihood evaluation.
  for (i in 1:state_dim+1) {
    y[i,] ~ normal(x_hat_add_CO2[i,], obs_error_scale * x_hat_add_CO2[i,]);
  }
}

generated quantities {
  array[N_t] vector<lower=0>[state_dim] x_hat_post_pred_intmd;
  array[state_dim] vector<lower=0>[N_t] x_hat_post_pred;
  vector<lower=0>[N_t] x_hat_post_pred_CO2;
  array[state_dim+1] vector<lower=0>[N_t] x_hat_post_pred_add_CO2;
  array[state_dim+1, N_t] real<lower=0> y_hat_post_pred;

  print("Iteration theta: ", "u_Q_ref = ", u_Q_ref, ", Q = ", Q, ", a_MSA = ", a_MSA, ", K_DE = ", K_DE, ", K_UE = ", K_UE, ", V_DE_ref = ", V_DE_ref, ", V_UE_ref = ", V_UE_ref, ", Ea_V_DE = ", Ea_V_DE, ", Ea_V_UE = ", Ea_V_UE, ", r_M = ", r_M, ", r_E = ", r_E, ", r_L = ", r_L);

  x_hat_post_pred_intmd = ode_rk45(AWB_ECA_ODE, x_hat0, t0, ts, u_Q_ref, Q, a_MSA, K_DE, K_UE, V_DE_ref, V_UE_ref, Ea_V_DE, Ea_V_UE, r_M, r_E, r_L, temp_ref, temp_rise);
  // Transform posterior predictive model output to match observations y in dimensions, [state_dim, N_t].
  for (i in 1:N_t) {
    for (j in 1:state_dim) {
      x_hat_post_pred[j, i] = x_hat_post_pred_intmd[i, j];
    }
  }

  // Compute posterior predictive CO2.
  x_hat_post_pred_CO2 = calc_AWB_ECA_CO2(ts, x_hat_post_pred, u_Q_ref, Q, K_UE, V_UE_ref, Ea_V_UE, temp_ref, temp_rise);   

  // Append CO2 vector to posterior predictive x_hat.
  x_hat_post_pred_add_CO2[1:4,] = x_hat_post_pred;
  x_hat_post_pred_add_CO2[5,] = x_hat_post_pred_CO2;

  // Add observation noise to posterior predictive model output to obtain posterior predictive samples.
  for (i in 1:state_dim+1) {
    y_hat_post_pred[i,] = normal_rng(x_hat_post_pred_add_CO2[i,], obs_error_scale * x_hat_post_pred_add_CO2[i,]);
  }
  print("Iteration posterior predictive y observation: ", y_hat_post_pred);

  // Obtain prior predictive samples. 
  real u_Q_ref_prior_pred = normal_lb_ub_rng(u_Q_ref_prior_dist_params[1], u_Q_ref_prior_dist_params[1] * prior_scale_factor, u_Q_ref_prior_dist_params[2], u_Q_ref_prior_dist_params[3]);
  real Q_prior_pred = normal_lb_ub_rng(Q_prior_dist_params[1], Q_prior_dist_params[1] * prior_scale_factor, Q_prior_dist_params[2], Q_prior_dist_params[3]);
  real a_MSA_prior_pred = normal_lb_ub_rng(a_MSA_prior_dist_params[1], a_MSA_prior_dist_params[1] * prior_scale_factor, a_MSA_prior_dist_params[2], a_MSA_prior_dist_params[3]);
  real K_DE_prior_pred = normal_lb_ub_rng(K_DE_prior_dist_params[1], K_DE_prior_dist_params[1] * prior_scale_factor, K_DE_prior_dist_params[2], K_DE_prior_dist_params[3]);
  real K_UE_prior_pred = normal_lb_ub_rng(K_UE_prior_dist_params[1], K_UE_prior_dist_params[1] * prior_scale_factor, K_UE_prior_dist_params[2], K_UE_prior_dist_params[3]);
  real V_DE_ref_prior_pred = normal_lb_ub_rng(V_DE_ref_prior_dist_params[1], V_DE_ref_prior_dist_params[1] * prior_scale_factor, V_DE_ref_prior_dist_params[2], V_DE_ref_prior_dist_params[3]);
  real V_UE_ref_prior_pred = normal_lb_ub_rng(V_UE_ref_prior_dist_params[1], V_UE_ref_prior_dist_params[1] * prior_scale_factor, V_UE_ref_prior_dist_params[2], V_UE_ref_prior_dist_params[3]);
  real Ea_V_DE_prior_pred = normal_lb_ub_rng(Ea_V_DE_prior_dist_params[1], Ea_V_DE_prior_dist_params[1] * prior_scale_factor, Ea_V_DE_prior_dist_params[2], Ea_V_DE_prior_dist_params[3]);
  real Ea_V_UE_prior_pred = normal_lb_ub_rng(Ea_V_UE_prior_dist_params[1], Ea_V_UE_prior_dist_params[1] * prior_scale_factor, Ea_V_UE_prior_dist_params[2], Ea_V_UE_prior_dist_params[3]);
  real r_M_prior_pred = normal_lb_ub_rng(r_M_prior_dist_params[1], r_M_prior_dist_params[1] * prior_scale_factor, r_M_prior_dist_params[2], r_M_prior_dist_params[3]);
  real r_E_prior_pred = normal_lb_ub_rng(r_E_prior_dist_params[1], r_E_prior_dist_params[1] * prior_scale_factor, r_E_prior_dist_params[2], r_E_prior_dist_params[3]);
  real r_L_prior_pred = normal_lb_ub_rng(r_L_prior_dist_params[1], r_L_prior_dist_params[1] * prior_scale_factor, r_L_prior_dist_params[2], r_L_prior_dist_params[3]);
}
