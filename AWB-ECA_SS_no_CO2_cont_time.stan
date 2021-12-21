functions {

  // Temperature function for ODE forcing.
  real temp_func(real t, real temp_ref) {
    return temp_ref + (5 * t) / (80 * 24 * 365) + 10 * sin((2 * pi() / 24) * t) + 10 * sin((2 * pi() / (24 * 365)) * t);
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

  real[] AWB_ECA_ODE(real t, real[] C, real[] theta, real[] x_r, int[] x_i) {

    // Initiate theta variables for future assignment.
    real u_Q_ref; // Reference carbon use efficiency.
    real Q; // Carbon use efficiency linear dependence negative slope.
    real a_MSA; // AWB MBC-to-SOC transfer fraction. 
    real V_DE_ref; // Reference SOC decomposition V_max 
    real V_UE_ref; // Reference DOC uptake V_max.
    real K_DE; // SOC decomposition K_m.
    real K_UE; // DOC uptake K_m.
    real Ea_V_DE; // SOC V_max activation energy.
    real Ea_V_UE; // DOC V_max activation energy.
    real r_M; // MBC turnover rate.
    real r_E; // Enzyme production rate.
    real r_L; // Enzyme loss rate.

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

    // Assign variables to thetas.
    u_Q_ref = theta[1];
    Q = theta[2];
    a_MSA = theta[3];
    V_DE_ref = theta[1];
    V_UE_ref = theta[2];
    K_DE = theta[3];
    K_UE = theta[4];
    Ea_V_DE = theta[5];
    Ea_V_UE = theta[6];
    r_M = theta[7];
    r_E = theta[8];
    r_L = theta[9];

    // Assign input and forcing variables to appropriate value at time t.
    temp = temp_func(t, x_r[1]); // x_r[1] is temp_ref 283.
    i_s = i_s_func(t);
    i_d = i_d_func(t);

    // Force temperature dependent parameters.
    u_Q = linear_temp(u_Q_ref, temp, Q, x_r[1]);
    V_DE = arrhenius_temp(V_DE_ref, temp, Ea_V_DE, x_r[1]);
    V_UE = arrhenius_temp(V_UE_ref, temp, Ea_V_UE, x_r[1]);

    // Assign dependent expressions.
    F_S = V_DE * C[4] * C[1] / (K_DE + C[4] + C[1]);
    F_D = V_UE * C[3] * C[2] / (K_UE + C[3] + C[2]);

    // Compute derivatives.
    vector[4] dCdt; // Initiate vector for storing derivatives.
    dCdt[1] = i_s + a_MSA * r_M * C[3] - F_S;
    dCdt[2] = i_d + (1 - a_MSA) * r_M * C[3] + F_S + r_L * C[4] - F_D;
    dCdt[3] = u_Q * F_D * - (r_M + r_E) * C[3];
    dCdt[4] = r_E * C[3] - r_L * C[4];

    return to_array_1d(dCdt);
  }
}

data {
}

transformed data {
}

parameters {
}

model {
}

generated quantities {
}
