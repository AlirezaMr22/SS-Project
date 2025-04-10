function [w_total, w_profile, w_annulus, w_sec] = calculate_loss_coefficients(beta1, beta2, Cm1, sigma, c, H, incidence, istar, istar_low_mach, ic, is, rm1, rm2)

beta1_rad = deg2rad(beta1);
beta2_rad = deg2rad(beta2);

Wm1 = Cm1;
Wm2 = Cm1; 

Wteta1 = Wm1*tan(beta1_rad);
Wteta2 = Wm2*tan(beta2_rad);
% Wteta1 = Cm1 * tan(beta1_rad);
% Wteta2 = Cm1 * tan(beta2_rad);
K = 0.0117;

Wmax_W1 =real( 1.12 + (0.61 * (cos(beta1_rad)^2 / sigma) * (((rm1*Wteta1) - (rm2*Wteta2)) / (rm1*Wm1))) + K * (incidence - istar)^1.43);

%%  Profile Loss (omega_profile)

W1 = Wm1 / cos(beta1_rad);
W2 = Wm2 / cos(beta2_rad);

Deq = Wmax_W1 / (W1/W2);

% Calculate Profile Loss Coefficient (omega_profile)
term1_profile = 2 * sigma / cos(beta2_rad);
term2_profile = (W2 / W1)^2; % (W2*/W1*)^2 term
term3_profile = 0.004 * (1 + (3.1 * (Deq - 1)^2) + (0.4 * (Deq - 1)^8)); 

w_profile = (term1_profile * term2_profile * term3_profile);

%% Annulus Loss (omega_annulus)

% Calculate Blade Spacing (s)
s = c / sigma;

% Calculate Annulus Drag Coefficient
CDannulus = 0.02 * (s / H);

tan_beta_m = 0.5 * (tan(beta1_rad) + tan(beta2_rad));
beta_m_rad = atan(tan_beta_m);

% Calculate Annulus Loss Coefficient (omega_annulus)
cos_beta1_squared = cos(beta1_rad)^2; % More efficient
cos_beta_m_cubed = cos(beta_m_rad)^3;
w_annulus = (sigma * CDannulus * cos_beta1_squared) / (cos_beta_m_cubed); 

%% Secondary Loss (omega_sec)

% Calculate Lift Coefficient (CL)

CDprofile_annulus = (w_profile + w_annulus) * (cos_beta_m_cubed / (sigma * cos_beta1_squared));
CL_term1 = (2 / sigma) * (tan(beta1_rad) - tan(beta2_rad)) * cos(beta_m_rad);
CL_term2 = CDprofile_annulus * tan(beta_m_rad); % 
CL = CL_term1 - CL_term2; % 

% Calculate Secondary Drag Coefficient (CDsec)
CDsec = 0.018 * CL^2;

% Calculate Secondary Loss Coefficient (omega_sec)
w_sec = (sigma * CDsec * cos_beta1_squared) / (cos_beta_m_cubed); % Same factor as annulus

%% Total Loss 
w_total_star = w_profile + w_annulus + w_sec; 

%% off design loss

% if incidence >= istar
%     landa = (incidence - istar ) / (is - istar);
%     w_total = w_total_star * (1+ (landa)^2);
% 
% else
%     landa = (incidence - istar ) / (istar - ic);
%     w_total = w_total_star * (1+ (landa)^2);
% 
% end % End of function   

%% high mach loss
Rs = is - istar;
w_total_star = w_total_star * (1+ ((istar - istar_low_mach)^2)/(Rs^2));

if incidence >= istar
    w_total = w_total_star + w_total_star * ((incidence - istar)/(is - istar))^2;
else
    w_total = w_total_star + w_total_star * ((incidence - istar)/(ic - istar))^2;
end
