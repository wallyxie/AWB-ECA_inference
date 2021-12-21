functions {

// Temperature function for ODE forcing.
real temp_func(real t, real temp_ref) {
return temp_ref + (5 * t) / (80 * 24 * 365) + 10 * sin((2 * pi() / 24) * t) + 10 * sin((2 * pi() / (24 * 365)) * t);
}

// Exogenous SOC input function.
real i_s_func(real t) {
return 0.001 + 0.0005 * sin((2 * pi() / (24 * 365)) * t)
}

// Exogenous DOC input function.
real i_d_func(real t) {
return 0.0001 + 0.00005 * sin((2 * pi() / (24 * 365)) * t)
}

// Function for enforcing Arrhenius temperature dependency of ODE parameter.
real arrhenius_temp(real input, real Ea, real temp, real temp_ref) {
return input * exp(-Ea / 0.008314 * (1 / temp - 1 / temp_ref));
}

// Function for enforcing linear temperature dependency of ODE parameter.
real linear_temp(real input, real Q, real temp, real temp_ref) {
return input - Q * (temp - temp_ref);
}

/*
AWB-ECA (equilibrium chemistry approximation) variant of AWB model.
C[1] is soil organic carbon (SOC) density.
C[2] is dissolved organic carbon (DOC) density.
C[3] is microbial biomass carbon (MBC) density.
C[4] is extracellular enzyme carbon (EEC) density.
*/

real[] AWB_ECA_ODE(real t, real[] C, real[] theta, real[] x_r, int[] x_i) {

real dC_dt[4];

// Initiate theta variables for future assignment.

real u_Q_ref; // Reference carbon use efficiency.
real Q; // Carbon use efficiency linear dependence negative slope.
real a_MSA; // AWB MBC-to-SOC transfer fraction. 
real V_DE_ref; // Reference SOC decomposition V_max 
real V_UE_ref; // Reference DOC uptake V_max.
real K_DE; // SOC decomposition K_m.
real K_UE; // DOC uptake K_m.
real Ea_V; // SOC V_max activation energy
real Ea_VU; // DOC V_max activation energy
real Ea_K; // SOC K_m activation energy
real Ea_KU; // DOC K_m activation energy 
real r_M; // MBC turnover rate
real r_E; // Enzyme production rate
real r_L; // Enzyme loss rate
real a_MS; // Fraction of dead MBC transferred to SOC
real E_C_f; // Reference carbon use efficiency (CUE)
real m_t; // CUE slope


// Dependent expressions

real F_S; // SOC decomposition
real F_D; // DOC uptake
real F_M; // MBC death
real F_E; // Enzyme production
real F_L; // Enzyme loss

/*
*Match parameters to thetas
*/

V_DE_ref = theta[1]; // SOC reference V_max
V_UE_ref = theta[2]; // DOC reference V_max
K_DE_ref = theta[3]; // SOC reference K_m
K_UE_ref = theta[4]; // DOC reference K_m
Ea_V = theta[5]; // SOC V_max activation energy
Ea_VU = theta[6]; // DOC V_max activation energy
Ea_K = theta[7]; // SOC K_m activation energy
Ea_KU = theta[8]; // DOC K_m activation energy 
r_M = theta[9]; // MBC turnover rate
a_MS = theta[10]; // Fraction of dead MBC transferred to SOC
E_C_f = theta[11]; // Carbon use efficiency

r_L = x_r[FILL IN HERE];
r_E = x_r[FILL IN HERE];

F_S = V_f * C[4] * C[1] / (K_f + C[1]);
F_D = V_U_f * C[3] * C[2] / (K_U_f + C[2]);
F_M = r_M * C[3];
F_E = r_E * C[3];
F_L = r_L * C[4];

dC_dt[1] = x_r[3] + a_MS * F_M - F_S;
dC_dt[2] = x_r[4] + (1 - a_MS) * F_M + F_S + F_L - F_D;
dC_dt[3] = F_D * E_C_f - F_M - F_E;
dC_dt[4] = F_E - F_L;

return dC_dt;}
