function stator_exit_props = calculate_stator_exit(...
    P02, h02, P2, alpha2, alphap2, Cm2, m_dot, rh, rt, alphap3, delta, w_total_stator)

%

fprintf('\n--- Calculating Stator Exit Properties ---\n');

% --- checking stator incidence angle---
is_stator = alpha2 - alphap2;
fprintf('Stator Incidence Angle (is_stator): %.2f deg\n', is_stator);

% --- Constants for iteration ---
tol = 0.001;     % Convergence tolerance
max_iter = 100;  % Maximum iterations
% --- Initialize iteration ---
Cm3 = Cm2; % Initial guess
converged = false;
iter = 0;

while ~converged && iter < max_iter
    iter = iter + 1;
    
    % --- Calculate properties with current Cm3 ---
    alpha3 = alphap3 + delta;
    alpha3_rad = deg2rad(alpha3);
    
    C3 = Cm3 / cos(alpha3_rad);
    
    P03 = P02 - w_total_stator * (P02 - P2);
    h03 = h02;
    s3 = thermodynamic_calculator('S','H',h03,'P',P03);
    
    h3 = h03 - (C3^2 / 2);
    
    T3 = thermodynamic_calculator('T','H',h3,'S',s3);
    P3 = thermodynamic_calculator('P','T',T3,'S',s3);
    rho3 = thermodynamic_calculator('D','T',T3,'P',P3);
    
    % --- Calculate new Cm3 ---
    rm3 = 0.5 * (rh(3) + rt(3));
    H_out = rt(3) - rh(3);
    A3 = 2 * pi * rm3 * H_out;
    Cm3_new = m_dot / (rho3 * A3);
    
    % --- Check convergence and update Cm3 ---
    diff = abs(Cm3_new - Cm3);
    
    if diff < tol
        converged = true;
        fprintf('Converged after %d iterations.\n', iter);
    else
        Cm3 = Cm3_new; % Direct update for next iteration
    end
end

if ~converged
    warning('Did not converge after %d iterations. Final difference: %.6f', max_iter, abs(Cm3_new - Cm3));
end

% --- Final calculations with converged Cm3 ---
alpha3 = alphap3 + delta;
alpha3_rad = deg2rad(alpha3);
C3 = Cm3 / cos(alpha3_rad);
h3 = h03 - (C3^2 / 2);

% --- Store results ---
stator_exit_props = struct(...
    'Cm3', Cm3, ...
    'alpha3_deg', alpha3, ...
    'C3', C3, ...
    'rm3', rm3, ...
    'P3', P3, ...
    'T3', T3, ...
    'rho3', rho3, ...
    'h3', h3, ...
    's3', s3, ...
    'P03', P03, ...
    'h03', h03, ...
    'iterations', iter, ...
    'converged', converged ...
);

% --- Display results ---
fprintf('\n--- Final Stator Exit Properties ---\n');
fprintf('Outlet Axial Velocity (Cm3): %.4f m/s\n', Cm3);
fprintf('Outlet Absolute Flow Angle (alpha3): %.2f deg\n', alpha3);
fprintf('Outlet Absolute Velocity (C3): %.2f m/s\n', C3);
fprintf('Outlet Static Pressure (P3): %.2f Pa\n', P3);
fprintf('Outlet Static Temperature (T3): %.2f K\n', T3);
fprintf('Outlet Density (rho3): %.4f kg/m^3\n', rho3);
fprintf('Number of iterations: %d\n', iter);
fprintf('Converged: %s\n', string(converged));