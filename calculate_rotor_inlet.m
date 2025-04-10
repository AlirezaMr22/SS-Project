function [rotor_inlet_state] = calculate_rotor_inlet(P01, T01, alpha1, beta1, betap1, rh, rt, RPM)

    h01 = thermodynamic_calculator('H', 'P', P01, 'T', T01); % J/kg
    s1 = thermodynamic_calculator('S', 'P', P01, 'T', T01); % J/kg-K
    rho01 = thermodynamic_calculator('D', 'P', P01, 'T', T01); % kg/m^3


%%  rm for first stage only

rm_in = 0.5 * (rh(1,1) + rt(1,1));
omega_rad_s = RPM * (2 * pi / 60); % Rotational speed in rad/s
U1 = omega_rad_s * rm_in;

beta1_rad = deg2rad(beta1);
alpha1_rad = deg2rad(alpha1);
tan_alpha1 = tan(alpha1_rad);
tan_beta1 = tan(beta1_rad);

Cm1 = U1 / (-tan_alpha1 + tan_beta1);

C1 = Cm1 / cos(alpha1_rad);

h1 = h01 - (0.5 * (C1^2)); % J/kg

    T1 = thermodynamic_calculator('T', 'H', h1, 'S', s1); % K
    P1 = thermodynamic_calculator('P', 'H', h1, 'S', s1); % Pa
    rho1 = thermodynamic_calculator('D', 'H', h1, 'S', s1); % kg/m^3
    
W1 = Cm1 / cos(beta1_rad);    
h01_rel = h1 +(0.5 * (W1^2)); % J/kg

    T01_rel = thermodynamic_calculator('T', 'H', h01_rel, 'S', s1); % K
    P01_rel = thermodynamic_calculator('P', 'H', h01_rel, 'S', s1); % Pa
    
H_in = rt(1,1) - rh(1,1);

A1 = 2 * pi * rm_in * H_in ;

m_dot = rho1 * A1 * Cm1;
    
% --- Store outputs in a structure ---
rotor_inlet_state = struct();
rotor_inlet_state.P01 = P01;
rotor_inlet_state.T01 = T01;
rotor_inlet_state.h01 = h01;
rotor_inlet_state.s1 = s1;
rotor_inlet_state.rho01 = rho01;
rotor_inlet_state.alpha1 = alpha1;
rotor_inlet_state.beta1 = beta1;
rotor_inlet_state.betap1 = betap1; 
rotor_inlet_state.U1 = U1;
rotor_inlet_state.Cm1 = Cm1;
rotor_inlet_state.C1 = C1;
rotor_inlet_state.P1 = P1;
rotor_inlet_state.T1 = T1;
rotor_inlet_state.h1 = h1;
rotor_inlet_state.s1 = s1;
rotor_inlet_state.rho1 = rho1;
rotor_inlet_state.W1 = W1;
rotor_inlet_state.P01_rel = P01_rel;
rotor_inlet_state.T01_rel = T01_rel;
rotor_inlet_state.h01_rel = h01_rel;
rotor_inlet_state.m_dot = m_dot;
