library(cmdstanr)
library(posterior)
library(bayesplot)

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
y <- y_full %>% filter(hour <= t)
ts <- y$hour
N_t <- length(ts)
y <- y %>% select(-hour)
y <- split(y, 1:nrow(y))

#Parameter prior means
u_Q_ref_mean <- 0.25
Q_mean <- 0.001
a_MSA_mean <- 0.5
K_DE_mean <- 1000
K_UE_mean <- 0.1
V_DE_ref_mean <- 0.11
V_UE_ref_mean <- 0.0044
Ea_V_DE_mean <- 40
Ea_V_UE_mean <- 30
r_M_mean <- 0.0001667
r_E_mean <- 0.0002
r_L_mean <- 0.0004

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
    u_Q_ref_mean = u_Q_ref_mean,
    Q_mean = Q_mean,
    a_MSA_mean = a_MSA_mean,
    K_DE_mean = K_DE_mean,
    K_UE_mean = K_UE_mean,
    V_DE_ref_mean = V_DE_ref_mean,
    V_UE_ref_mean = V_UE_ref_mean,
    Ea_V_DE_mean = Ea_V_DE_mean,
    Ea_V_UE_mean = Ea_V_UE_mean,
    r_M_mean = r_M_mean,
    r_E_mean = r_E_mean,
    r_L_mean = r_L_mean
    )

file_path <- 'AWB-ECA_SS_no_CO2_cont_time.stan' #Read in Stan model code.
lines <- readLines(file_path, encoding = "ASCII")
for (n in 1:length(lines)) cat(lines[n],'\n')
model <- cmdstan_model(file_path)

AWB_ECA_stan_fit <- model$sample(data = data_list, iter = 1000, warmup = 500, refresh = 10, chains = 3, seed = 123)
