% MATLAB code defining variables and arrays

% --- General Compressor Parameters ---
N = 1;                     % Number of stages 
P01 = 101325;         % Inlet stagnation pressure (Pa)
T01 = 288;            % Inlet stagnation temperature (K)
alpha1 = 15;      % Absolute inlet flow angle (degrees) - Assumed Axial (No IGV)
NumSpeeds = 5;             % Number of speed lines 
RPM = [0.8*16000,0.9*16000,16000,1.1*16000,1.2*16000];      % Array of RPM values 

% Radii at 3 sections (Rotor Inlet, Rotor Outlet, Stator Outlet) 
% Rows = stages, Columns = sections 


rh = [0.1244, 0.1471, 0.1602];% 0.1624, 0.1734, 0.1811; 0.1822, 0.1877, 0.1917]; %  Hub Radius (m)
rt = [0.2540, 0.2473, 0.2458];% 0.2455, 0.2407, 0.2396; 0.2393, 0.2349, 0.2327]; % Tip Radius (m)



% --- Stage-Wise Blade Parameters (N=3 stages) ---

% Blade Angles (at mean radius, degrees) 
% betap1 = [54, 53, 53];    %  - Rotor inlet blade angle
betap1 = 54;
% betap2 = [39, 37, 37];    %  - Rotor outlet blade angle
betap2 = 39;
% alphap2 = [40, 39, 40];    % - Stator inlet blade angle
alphap2 = 40;
% alphap3 = [-10, -9, -11];   %  - Stator outlet blade angle
alphap3 = -10;

% Solidity 
% sigma_R = [1.16, 1.41, 1.32];       %  Rotor solidity
sigma_R = 1.16;
% sigma_S = [1.45, 1.365, 1.332];     %  Stator solidity
sigma_S = 1.45;

% Chord Length (meters) 
% Chord_R_cm = [7.57, 5.955, 4.62];
Chord_R_cm = 7.57;
% Chord_S_cm = [5.70, 4.019, 3.369];
Chord_S_cm = 5.70;

Chord_R = Chord_R_cm / 100;       % Rotor chord length (m)
Chord_S = Chord_S_cm / 100;       % Stator chord length (m)

% Thickness (Max thickness-to-chord ratio, t/c) 
% tb_R_cm = [0.368, 0.362, 0.311]; % Max thickness in cm
tb_R_cm = 0.368;
% tb_S_cm = [0.329, 0.274, 0.265]; % Max thickness in cm
tb_S_cm = 0.329;

tb_c_R = tb_R_cm ./ Chord_R_cm;   % Rotor max thickness/chord ratio (Dimensionless)
tb_c_S = tb_S_cm ./ Chord_S_cm;   % Stator max thickness/chord ratio (Dimensionless)

% Stagger Angle (degrees) - NOT EXPLICITLY GIVEN in PDF
% Needs to be calculated or assumed based on profile definition/other data.
% Initialize as NaN or Zero, add comment.
% Stagger_R = nan(1, N); % ?r - Rotor stagger angle 
% Stagger_S = nan(1, N); % ?s - Stator stagger angle 

% Camber Angle (degrees) - 
theta_R =abs( betap1 - betap2);   %  - Rotor camber angle 
theta_S = abs(alphap2 - alphap3);   %  - Stator camber angle 

% TipClearance_R = ones(1, N) * 0.0005; % Rotor tip clearance (Taw_R) (Assumed 0.5 mm)
% TipClearance_S = zeros(1, N);         % Stator tip/hub clearance (Taw_S) (Assumed 0 for shrouded)

% --- Display Example Variable Values (for verification) ---
% disp('--- General Parameters ---');
% disp(['N = ', num2str(N)]);
% disp(['P0_inlet (Pa) = ', num2str(P01)]);
% disp(['T0_inlet (K) = ', num2str(T01)]);
% disp(['Alpha1_inlet (deg) = ', num2str(alpha1)]);
% disp(['RPM = ', num2str(RPM)]); % Displaying the design RPM
% 
% % disp('--- Stage 1 Data ---');
% % disp(['Blade Type R1: ', Blade_type_R{1}, ', S1: ', Blade_type_S{1}]);
% disp(['rh (m): RotorIn=', num2str(rh(1,1)), ' RotorOut=', num2str(rh(1,2)), ' StatorOut=', num2str(rh(1,3))]);
% disp(['rt (m): RotorIn=', num2str(rt_(1,1)), ' RotorOut=', num2str(rt_(1,2)), ' StatorOut=', num2str(rt_(1,3))]);
% disp(['Betap1 (deg): ', num2str(Betap1_deg(1)), ', Betap2 (deg): ', num2str(Betap2_deg(1))]);
% disp(['Alfap2 (deg): ', num2str(Alfap2_deg(1)), ', Alfap3 (deg): ', num2str(Alfap3_deg(1))]);
% disp(['Sigma R: ', num2str(Sigma_R(1)), ', Sigma S: ', num2str(Sigma_S(1))]);
% disp(['Chord R (m): ', num2str(Chord_R(1)), ', Chord S (m): ', num2str(Chord_S(1))]);
% disp(['Thickness R (t/c): ', num2str(Thickness_R(1)), ', Thickness S (t/c): ', num2str(Thickness_S(1))]);
% disp(['Theta R (deg): ', num2str(Theta_R_deg(1)), ', Theta S (deg): ', num2str(Theta_S_deg(1))]);
% % disp(['Tip Clearance R (m): ', num2str(TipClearance_R(1)), ', Tip Clearance S (m): ', num2str(TipClearance_S(1))]);
% % Add similar display blocks for Stage 2 and Stage 3 if needed