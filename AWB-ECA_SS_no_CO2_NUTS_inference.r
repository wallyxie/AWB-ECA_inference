library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

color_scheme_set("viridisA")

#Data to be passed to Stan.
state_dim <- 4
temp_ref <- 283
temp_rise <- 5 #High estimate of 5 celsius temperature rise by 2100.
prior_scale_factor <- 0.25
obs_error_scale <- 0.1
obs_every <- 10 #Observations every 10 hours.
t <- 3000 #Total time span of ODE simulation.
x_hat0 <- c(99.69113179, 1.94470102, 1.87967951, 2.01791824) #Originally sampled values used for Euler-Maruyama solution.
y_full <- read_csv('generated_data/SAWB-ECA-SS_no_CO2_trunc_short_2021_12_23_22_21_sample_y_t_3000_dt_0-01_sd_scale_0-25.csv')
y <- y_full %>% filter(hour <= t) %>% tail(-1)
ts <- y$hour
N_t <- length(ts)
y <- y %>% select(-hour)
#y <- split(y, 1:nrow(y)) #Convert data observations to list of rows to correspond to Stan's array of vectors type.
y <- as.list(y) #Convert data observations to list of columns to correspond to Stan's array of vectors type.

#Parameter prior means
u_Q_ref_prior_dist_params <- c(0.22, 1e-2, 1)
Q_prior_dist_params <- c(0.001, 0, 0.1)
a_MSA_prior_dist_params <- c(0.5, 0, 1)
K_DE_prior_dist_params <- c(1000, 100, 5000)
K_UE_prior_dist_params <- c(0.1, 1e-2, 1)
V_DE_ref_prior_dist_params <- c(0.04, 1e-3, 1)
V_UE_ref_prior_dist_params <- c(0.005, 1e-4, 0.1)
Ea_V_DE_prior_dist_params <- c(40, 5, 80)
Ea_V_UE_prior_dist_params <- c(30, 5, 80)
r_M_prior_dist_params <- c(0.00016667, 1e-5, 0.1)
r_E_prior_dist_params <- c(0.0002, 1e-5, 0.1)
r_L_prior_dist_params <- c(0.0004, 1e-5, 0.1)

#Create list of lists to pass prior means as initial theta values in Stan corresponding to four chains.
init_theta_single = list(
                  u_Q_ref = u_Q_ref_prior_dist_params[1],
                  Q = Q_prior_dist_params[1],
                  a_MSA = a_MSA_prior_dist_params[1],
                  K_DE = K_DE_prior_dist_params[1],
                  K_UE = K_UE_prior_dist_params[1],
                  V_DE_ref = V_DE_ref_prior_dist_params[1],
                  V_UE_ref = V_UE_ref_prior_dist_params[1],
                  Ea_V_DE = Ea_V_DE_prior_dist_params[1],
                  Ea_V_UE = Ea_V_UE_prior_dist_params[1],
                  r_M = r_M_prior_dist_params[1],
                  r_E = r_E_prior_dist_params[1],
                  r_L = r_L_prior_dist_params[1]
                  )
init_theta = list(init_theta_single)[rep(1, 4)]

data_list = list(
    state_dim = state_dim,
    N_t = N_t,
    ts = ts,
    y = y,
    temp_ref = temp_ref,
    temp_rise = temp_rise,
    prior_scale_factor = prior_scale_factor,
    obs_error_scale = obs_error_scale,
    x_hat0 = x_hat0,
    u_Q_ref_prior_dist_params = u_Q_ref_prior_dist_params,
    Q_prior_dist_params = Q_prior_dist_params,
    a_MSA_prior_dist_params = a_MSA_prior_dist_params,
    K_DE_prior_dist_params = K_DE_prior_dist_params,
    K_UE_prior_dist_params = K_UE_prior_dist_params,
    V_DE_ref_prior_dist_params = V_DE_ref_prior_dist_params,
    V_UE_ref_prior_dist_params = V_UE_ref_prior_dist_params,
    Ea_V_DE_prior_dist_params = Ea_V_DE_prior_dist_params,
    Ea_V_UE_prior_dist_params = Ea_V_UE_prior_dist_params,
    r_M_prior_dist_params = r_M_prior_dist_params,
    r_E_prior_dist_params = r_E_prior_dist_params,
    r_L_prior_dist_params = r_L_prior_dist_params
    )

file_path <- 'AWB-ECA_SS_no_CO2_cont_time.stan' #Read in Stan model code.
lines <- readLines(file_path, encoding = "ASCII")
for (n in 1:length(lines)) cat(lines[n],'\n')
model <- cmdstan_model(file_path)

AWB_ECA_stan_fit <- model$sample(data = data_list, seed = 1234, refresh = 20, init = init_theta, iter_sampling = 7500, iter_warmup = 1100, chains = 4, parallel_chains = 4, adapt_delta = 0.95)
