functions {

  // Temperature function for ODE forcing.
  real temp_func(real t, real temp_ref, real temp_rise) {
    return temp_ref + (temp_rise * t) / (80 * 24 * 365) + 10 * sin((2 * pi() / 24) * t) + 10 * sin((2 * pi() / (24 * 365)) * t);
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

  // Function for enforcing linear temperature dependency of ODE parameter.
  real linear_temp(real input, real temp, real Q, real temp_ref) {
    return input - Q * (temp - temp_ref);
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

    // Compute derivatives.
    vector[4] dCdt; // Initiate vector for storing derivatives.
    dCdt[1] = i_s + a_MSA * r_M * C[3] - F_S;
    dCdt[2] = i_d + (1 - a_MSA) * r_M * C[3] + F_S + r_L * C[4] - F_D;
    dCdt[3] = u_Q * F_D * - (r_M + r_E) * C[3];
    dCdt[4] = r_E * C[3] - r_L * C[4];

    return dCdt;
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
  array[N_t] vector<lower=0>[state_dim] y; // Multidimensional array of state observations bounded at 0.
  real<lower=0> temp_ref; // Reference temperature for temperature forcing.
  real<lower=0> temp_rise; // Assumed increase in temperature (°C/K) over next 80 years.
  real<lower=0> prior_scale_factor; // Factor multiplying parameter means to obtain prior standard deviations.
  real<lower=0> obs_error_scale; // Observation noise factor multiplying observations of model output x_hat.
  vector<lower=0>[state_dim] x_hat0; // Initial ODE conditions. [60.35767864, 5.16483124, 2.0068896, 0.99331202] for AWB-ECA instance of data.
  real<lower=0> u_Q_ref_mean;
  real<lower=0> Q_mean;
  real<lower=0> a_MSA_mean;
  real<lower=0> K_DE_mean;
  real<lower=0> K_UE_mean;
  real<lower=0> V_DE_ref_mean;
  real<lower=0> V_UE_ref_mean;
  real<lower=0> Ea_V_DE_mean;
  real<lower=0> Ea_V_UE_mean;
  real<lower=0> r_M_mean;
  real<lower=0> r_E_mean;
  real<lower=0> r_L_mean;
}

transformed data {
  real t0 = 0; // Initial time.
}

parameters {
  real<lower=0> u_Q_ref; // Reference carbon use efficiency.
  real<lower=0> Q; // Carbon use efficiency linear dependence negative slope.
  real<lower=0> a_MSA; // AWB MBC-to-SOC transfer fraction.
  real<lower=0> K_DE; // SOC decomposition K_m.
  real<lower=0> K_UE; // DOC uptake K_m.
  real<lower=0> V_DE_ref; // Reference SOC decomposition V_max.
  real<lower=0> V_UE_ref; // Reference DOC uptake V_max.
  real<lower=0> Ea_V_DE; // SOC V_max activation energy.
  real<lower=0> Ea_V_UE; // DOC V_max activation energy.
  real<lower=0> r_M; // MBC turnover rate.
  real<lower=0> r_E; // Enzyme production rate.
  real<lower=0> r_L; // Enzyme loss rate.
}

transformed parameters {
  array[N_t] vector<lower=0>[state_dim] x_hat = ode_ckrk(AWB_ECA_ODE, x_hat0, t0, ts, u_Q_ref, Q, a_MSA, K_DE, K_UE, V_DE_ref, V_UE_ref, Ea_V_DE, Ea_V_UE, r_M, r_E, r_L, temp_ref, temp_rise); 
}

model {
  u_Q_ref ~ normal(u_Q_ref_mean, u_Q_ref_mean * prior_scale_factor) T[0, 1];
  Q ~ normal(Q_mean, Q_mean * prior_scale_factor) T[0, 0.1];
  a_MSA ~ normal(a_MSA_mean, a_MSA_mean * prior_scale_factor) T[0, 1];
  K_DE ~ normal(K_DE_mean, K_DE_mean * prior_scale_factor) T[0, 5000];
  K_UE ~ normal(K_UE_mean, K_UE_mean * prior_scale_factor) T[0, 1];
  V_DE_ref ~ normal(V_DE_ref_mean, V_DE_ref_mean * prior_scale_factor) T[0, 1];
  V_UE_ref ~ normal(V_UE_ref_mean, V_UE_ref_mean * prior_scale_factor) T[0, 0.1];
  Ea_V_DE ~ normal(Ea_V_DE_mean, Ea_V_DE_mean * prior_scale_factor) T[5, 80];
  Ea_V_UE ~ normal(Ea_V_UE_mean, Ea_V_UE_mean * prior_scale_factor) T[5, 80];
  r_M ~ normal(r_M_mean, r_M_mean * prior_scale_factor) T[0, 0.1];
  r_E ~ normal(r_E_mean, r_E_mean * prior_scale_factor) T[0, 0.1];
  r_L ~ normal(r_L_mean, r_L_mean * prior_scale_factor) T[0, 0.1];

  // Likelihood evaluation.
  // y ~ normal(x_hat, obs_error_scale * x_hat);
  y ~ normal(x_hat, 1);
}

generated quantities {
  array[N_t] vector<lower=0>[state_dim] x_hat_post_pred;
  array[N_t] vector<lower=0>[state_dim] y_hat_post_pred;
  x_hat_post_pred = ode_ckrk(AWB_ECA_ODE, x_hat0, t0, ts, u_Q_ref, Q, a_MSA, K_DE, K_UE, V_DE_ref, V_UE_ref, Ea_V_DE, Ea_V_UE, r_M, r_E, r_L);
  y_hat_post_pred = normal_rng(x_hat_post_pred, obs_error_scale * x_hat_post_pred);
  real u_Q_ref_post = normal_lb_ub_rng(u_Q_ref_mean, u_Q_ref_mean * prior_scale_factor, 0, 1);
  real Q_post = normal_lb_ub_rng(Q_mean, Q_mean * prior_scale_factor, 0, 0.1);
  real a_MSA_post = normal_lb_ub_rng(a_MSA_mean, a_MSA_mean * prior_scale_factor, 0, 1);
  real K_DE_post = normal_lb_ub_rng(K_DE_mean, K_DE_mean * prior_scale_factor, 0, 5000);
  real K_UE_post = normal_lb_ub_rng(K_UE_mean, K_UE_mean * prior_scale_factor, 0, 1);
  real V_DE_ref_post = normal_lb_ub_rng(V_DE_ref_mean, V_DE_ref_mean * prior_scale_factor, 0, 1);
  real V_UE_ref_post = normal_lb_ub_rng(V_UE_ref_mean, V_UE_ref_mean * prior_scale_factor, 0, 0.1);
  real Ea_V_DE_post = normal_lb_ub_rng(Ea_V_DE_mean, Ea_V_DE_mean * prior_scale_factor, 5, 80);
  real Ea_V_UE_post = normal_lb_ub_rng(Ea_V_UE_mean, Ea_V_UE_mean * prior_scale_factor, 5, 80);
  real r_M_post = normal_lb_ub_rng(r_M_mean, r_M_mean * prior_scale_factor, 0, 0.1);
  real r_E_post = normal_lb_ub_rng(r_E_mean, r_E_mean * prior_scale_factor, 0, 0.1);
  real r_L_post = normal_lb_ub_rng(r_L_mean, r_L_mean * prior_scale_factor, 0, 0.1);
}
