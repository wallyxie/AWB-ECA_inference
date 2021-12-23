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
t <- 2000 #Total time span of ODE simulation.
dt <- 0.1 #Euler time step size corresponding to data generation. 
x_hat0 <- c(78.81405796, 5.66343552, 2.2209765, 1.19214226) #Originally sampled values used for Euler-Maruyama solution.
y_full <- read_csv('generated_data/SAWB-ECA-SS_no_CO2_trunc_short_2021_12_20_18_20_sample_y_t_2000_dt_0-01_sd_scale_0-25.csv')
y <- y_full %>% filter(hour <= t) %>% tail(-1)
ts <- y$hour
N_t <- length(ts)
y <- y %>% select(-hour)
#y <- split(y, 1:nrow(y)) #Convert data observations to list of rows to correspond to Stan's array of vectors type.
y <- as.list(y) #Convert data observations to list of columns to correspond to Stan's array of vectors type.

#Parameter prior means
u_Q_ref_prior_dist_params <- c(0.25, 0, 1)
Q_prior_dist_params <- c(0.001, 0, 0.1)
a_MSA_prior_dist_params <- c(0.5, 0, 1)
K_DE_prior_dist_params <- c(1000, 0, 5000)
K_UE_prior_dist_params <- c(0.1, 0, 1)
V_DE_ref_prior_dist_params <- c(0.11, 0, 1)
V_UE_ref_prior_dist_params <- c(0.0044, 0, 0.1)
Ea_V_DE_prior_dist_params <- c(40, 5, 80)
Ea_V_UE_prior_dist_params <- c(30, 5, 80)
r_M_prior_dist_params <- c(0.0001667, 0, 0.1)
r_E_prior_dist_params <- c(0.0002, 0, 0.1)
r_L_prior_dist_params <- c(0.0004, 0, 0.1)

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

AWB_ECA_stan_fit <- model$sample(data = data_list, iter_sampling = 20, iter_warmup = 10, refresh = 1, chains = 3, parallel_chains = 3, seed = 123, adapt_delta = 0.95)
