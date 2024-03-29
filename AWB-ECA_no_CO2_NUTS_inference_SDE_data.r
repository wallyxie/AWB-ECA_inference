#Stan model coded for Stan version 2.29.1, not compatible with Stan version >= 2.32.0.
library(cmdstanr)
library(posterior)
library(tidyverse)
library(bayesplot)

options(mc.cores = parallel::detectCores())

num_chains <- 4

#Data to be passed to Stan.
state_dim <- 4
temp_ref <- 283
temp_rise <- 5 #High estimate of 5 celsius temperature rise by 2100.
prior_scale_factor <- 0.25
obs_error_scale <- 0.1
obs_every <- 10 #Observations every 10 hours.
t <- 2000 #Total time span of ODE simulation.
x_hat0 <- c(114.572059600027, 0.881891334730684, 1.51538971282248, 2.14394308159024) #Plugging observed y0 values for ODE solver x0.
y_full <- read_csv('generated_data/SAWB-ECA-SS_no_CO2_trunc_short_2021_12_24_12_46_sample_y_t_2000_dt_0-01_sd_scale_0-25.csv')
y <- y_full %>% filter(hour <= t) %>% tail(-1)
ts <- y$hour
N_t <- length(ts)
y <- y %>% select(-hour)
#y <- split(y, 1:nrow(y)) #Convert data observations to list of rows to correspond to Stan's array of vectors type.
y <- as.list(y) #Convert data observations to list of columns to correspond to Stan's array of vectors type.

#Parameter prior means
u_Q_ref_prior_dist_params <- c(0.22, 0, 1)
Q_prior_dist_params <- c(0.001, 0, 0.1)
a_MSA_prior_dist_params <- c(0.5, 0, 1)
K_DE_prior_dist_params <- c(1000, 0, 5000)
K_UE_prior_dist_params <- c(0.1, 0, 1)
V_DE_ref_prior_dist_params <- c(0.04, 0, 1)
V_UE_ref_prior_dist_params <- c(0.005, 0, 0.1)
Ea_V_DE_prior_dist_params <- c(40, 5, 80)
Ea_V_UE_prior_dist_params <- c(30, 5, 80)
r_M_prior_dist_params <- c(0.00016667, 0, 0.1)
r_E_prior_dist_params <- c(0.0002, 0, 0.1)
r_L_prior_dist_params <- c(0.0004, 0, 0.1)

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
init_theta = list(init_theta_single)[rep(1, num_chains)] #num_chains copies of initial theta proposals for each of the HMC chains to be used.

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

file_path <- 'AWB-ECA_no_CO2_cont_time.stan' #Read in Stan model code.
lines <- readLines(file_path, encoding = "ASCII")
for (n in 1:length(lines)) cat(lines[n],'\n')
model <- cmdstan_model(file_path)

start = Sys.time()
AWB_ECA_stan_fit_no_CO2 <- model$sample(data = data_list, seed = 1234, refresh = 10, init = init_theta, iter_sampling = 1250, iter_warmup = 500, chains = num_chains, parallel_chains = num_chains, adapt_delta = 0.95)
stop = Sys.time()
elapsed = stop - start
cat(elapsed, file = 'NUTS_results/AWB-ECA_no_CO2_NUTS_inference_SDE_data_elapsed.txt')

#Save Stan fit object and NUTS inference results.
AWB_ECA_stan_fit_no_CO2$save_object(file = "NUTS_results/AWB-ECA_no_CO2_NUTS_inference_SDE_data.rds")
AWB_ECA_stan_fit_no_CO2_theta_post <- as_tibble(AWB_ECA_stan_fit_no_CO2$draws(c("u_Q_ref", "Q", "a_MSA", "K_DE", "K_UE", "V_DE_ref", "V_UE_ref", "Ea_V_DE", "Ea_V_UE", "r_M", "r_E", "r_L"), format = "draws_df"))
write_csv(AWB_ECA_stan_fit_no_CO2_theta_post, "NUTS_results/AWB-ECA_no_CO2_NUTS_inference_SDE_data_theta_post.csv")
AWB_ECA_stan_fit_no_CO2_x_post <- as_tibble(AWB_ECA_stan_fit_no_CO2$draws(c("x_hat"), format = "draws_df"))
write_csv(AWB_ECA_stan_fit_no_CO2_x_post, "NUTS_results/AWB-ECA_no_CO2_NUTS_inference_SDE_data_x_post.csv")
AWB_ECA_stan_fit_no_CO2_post_summary <- as_tibble(AWB_ECA_stan_fit_no_CO2$summary(c("u_Q_ref", "Q", "a_MSA", "K_DE", "K_UE", "V_DE_ref", "V_UE_ref", "Ea_V_DE", "Ea_V_UE", "r_M", "r_E", "r_L", "x_hat")))
write_csv(AWB_ECA_stan_fit_no_CO2_post_summary, "NUTS_results/AWB-ECA_no_CO2_NUTS_inference_SDE_data_post_summary.csv")
