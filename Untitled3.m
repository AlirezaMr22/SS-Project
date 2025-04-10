clc
clear 
close all

% Main script for first stage 

% Step 1: Input Data Acquisition
data_input_design_point; % Load all parameters from data_input.m

% Step 2: Select Machine Speed and Inlet Angle
NumSpeeds = length(RPM);

% Step 3: Initialize Stage Counter 
N = 1;

% Number of incidence points to generate
num_inc = 10;  
    
% --- Calculate Incidence Angle (istar) ---
% Parameters for incidence calculation (example parameters are already provided)
tb_c = tb_c_R;       % Thickness-to-chord ratio for rotor
sigma  = sigma_R;    % Rotor solidity
theta  = theta_R;    % Camber angle

% Preallocate storage for results
% Dimensions: rows for each RPM, columns for each incidence condition
m_dot_store   = zeros(NumSpeeds, num_inc);   
PR            = zeros(NumSpeeds, num_inc);       
etha_is_stage = zeros(NumSpeeds, num_inc);      

% Main loop over RPM values
for rpm_idx = 1:NumSpeeds
    current_RPM = RPM(rpm_idx);
    fprintf('\n====== Processing RPM = %d (Index %d/%d) ======\n', current_RPM, rpm_idx, NumSpeeds);
    
    % Calculate the design incidence properties (istar, ic, is) at the current RPM
    [istar, ic, is, istar_low_mach] = incidence_calc(tb_c, sigma, alpha1, rh, rt, current_RPM, T01, betap1, theta);
    
    incidence_range = linspace(ic+1, is, num_inc);  % Equally spaced incidence angles
    
    % Loop through incidence angles for current RPM
    for inc_idx = 1:num_inc
        incidence = incidence_range(inc_idx);  % Current incidence angle
        fprintf('\n  -- Processing Incidence = %.2f degrees (Index %d/%d) --\n', incidence, inc_idx, num_inc);

        % --- Calculate Rotor Design Deviation Angle (delta_star) ---
        % Adjust inlet flow angle beta1 based on the current incidence
        beta1 = betap1 + incidence; % Inlet flow angle
        % Use deviation function to get the design deviation
        delta_star = deviation(tb_c, sigma, beta1, betap2, theta);
        % --- Calculate Rotor Inlet State ---
        [rotor_inlet_state] = calculate_rotor_inlet(P01, T01, alpha1, beta1, betap1, rh, rt, current_RPM);
        
        % Extract necessary variables from rotor inlet state
        Cm1 = rotor_inlet_state.Cm1;
        beta1_rad = deg2rad(rotor_inlet_state.beta1);
        rm1 = 0.5*(rh(1) + rt(1));
        rm2 = 0.5*(rh(2) + rt(2)); % Mean radius at rotor outlet
        
        % --- Calculate Rotor off-Design Deviation Angle ---
        [delta] = calculate_off_design_deviation(delta_star, istar, incidence, Cm1, sigma, beta1);
        beta2 = betap2 + delta;   % Rotor outlet blade angle

        % --- Calculate Rotor Loss Coefficients ---
        c = Chord_R;         % Rotor chord length (from data_input)
        H = (rt(2) - rh(2));   % Rotor blade height
        [w_total_rotor, ~, ~, ~] = calculate_loss_coefficients(beta1, beta2, Cm1, sigma, c, H, incidence, istar, istar_low_mach, ic, is, rm1, rm2);
        
        % Extract additional required variables from rotor inlet state for later use
        U1 = rotor_inlet_state.U1;
        h1 = rotor_inlet_state.h1;
        P01_rel = rotor_inlet_state.P01_rel;
        P1 = rotor_inlet_state.P1;
        h01 = rotor_inlet_state.h01;
        % Store current mass flow rate in a temporary variable
        m_dot_current = rotor_inlet_state.m_dot;
        
        % --- Calculate Rotor Outlet ---
        [rotor_exit_props] = calculate_rotor_outlet(...
            h01, P1, h1, P01_rel, m_dot_current, alpha1, beta1, betap2, delta, Cm1, U1, rh, rt, rm2, current_RPM, beta2, w_total_rotor);
        
        %% --- Calculate Stator Exit ---
        c = Chord_S;                      % Stator chord length (from data_input)
        H_stator = rt(3) - rh(3);           % Stator blade height
        H = H_stator;
        rm1 = 0.5*(rh(2) + rt(2));
        rm2 = 0.5*(rh(3) + rt(3));
        Cm2 = rotor_exit_props.Cm2;         % Axial velocity at stator inlet
        alpha2 = rotor_exit_props.alpha2;     % Stator inlet angle (rotor outlet absolute angle)
        beta1 = alpha2;
        betap2 = alphap3;
        theta = theta_S;                    % Stator camber angle

        % Calculate deviation for the stator stage
        tb_c = tb_c_S;     % Stator thickness-to-chord ratio
        sigma = sigma_S;   % Stator solidity
        
        delta_star = deviation(tb_c, sigma, beta1, betap2, theta);
        [delta] = calculate_off_design_deviation(delta_star, istar, incidence, Cm2, sigma, beta1);
        alpha3 = alphap3 + delta; % Stator outlet blade angle
        beta2 = alpha3;
        
        % --- Calculate Stator Loss Coefficients ---
        [w_total_stator, ~, ~, ~] = calculate_loss_coefficients(beta1, beta2, Cm1, sigma, c, H, incidence, istar, istar_low_mach, ic, is, rm1, rm2);
        
        % --- Calculate Stator Exit State ---
        % Extract necessary values from rotor exit properties for stator calculations
        P02 = rotor_exit_props.P02;
        P2 = rotor_exit_props.P2;
        h02 = rotor_exit_props.h02;
        [stator_exit_props] = calculate_stator_exit(...
            P02, h02, P2, alpha2, alphap2, Cm2, m_dot_current, rh, rt, alphap3, delta, w_total_stator);
        
        % --- Calculate Performance Metrics ---
        h2 = rotor_exit_props.h2;
        s1 = rotor_inlet_state.s1;
        P03 = stator_exit_props.P03;
        h03 = stator_exit_props.h03;
        h3 = stator_exit_props.h3;
        % Obtain isentropic head at stator exit (using a thermodynamic calculator function)
        h03ss = thermodynamic_calculator('H', 'P', P03, 'S', s1); % J/kg
        
        % Calculate current performance metrics
        PR_current = P03 / P01;
        etha_is_current = (h03ss - h01) / (h03 - h01);
        
        % Store the performance metrics for the current RPM and incidence point
        PR(rpm_idx, inc_idx) = PR_current;
        etha_is_stage(rpm_idx, inc_idx) = etha_is_current;
        m_dot_store(rpm_idx, inc_idx) = m_dot_current;
        
        % Print current calculation results summary
        fprintf('    Results: Mass Flow = %.4f kg/s, PR = %.4f, Efficiency = %.4f\n', m_dot_current, PR_current, etha_is_current);
    end % end incidence loop
    
    fprintf('\n  Completed RPM = %d calculations\n', current_RPM);
end % end RPM loop

fprintf('\n====== All calculations complete, generating plots ======\n');

%% Plot Results for each RPM

% Pressure Ratio vs Mass Flow Rate
figure;
hold on;
for i = 1:NumSpeeds
    plot(m_dot_store(i, :), PR(i, :), 'o-', 'LineWidth', 1.5);
end
xlabel('Mass Flow Rate (kg/s)');
ylabel('Pressure Ratio (P03/P01)');
title('Compressor Performance Map: Pressure Ratio vs Mass Flow Rate');
grid on;
legend(arrayfun(@(x) sprintf('RPM = %d', x), RPM, 'UniformOutput', false));
hold off;

% Isentropic Efficiency vs Mass Flow Rate
figure;
hold on;
for i = 1:NumSpeeds
    plot(m_dot_store(i, :), etha_is_stage(i, :), 'o-', 'LineWidth', 1.5);
end
xlabel('Mass Flow Rate (kg/s)');
ylabel('Isentropic Efficiency');
title('Compressor Performance Map: Isentropic Efficiency vs Mass Flow Rate');
grid on;
legend(arrayfun(@(x) sprintf('RPM = %d', x), RPM, 'UniformOutput', false));
hold off;