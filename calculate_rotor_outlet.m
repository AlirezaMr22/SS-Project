function [rotor_exit_props] = calculate_rotor_outlet(...
     h01, P1, h1, P01_rel, m_dot, alpha1, beta1, betap2, delta, Cm1, U1, rh, rt, rm2, RPM, beta2, w_total_rotor)

tol = 0.001; % Convergence tolerance
max_iter = 200; % Maximum iterations
converged = false;
iter = 0;
Cm2 = 0.8*Cm1; % Initial guess
Wm1 = Cm1;

while ~converged && iter < max_iter
    iter = iter + 1;
    
%     % Ensure Cm2 is real before using it
%     if ~isreal(Cm2)
%         Cm2 = real(Cm2);
%         if Cm2 <= 0
%             Cm2 = Cm1; % Reset to initial value if problematic
%         end
%     end
    
    Wm2 = Cm2;
    beta2 = betap2 + delta;  
    beta2_rad = deg2rad(beta2);
    alpha1_rad = deg2rad(alpha1);
    beta1_rad = deg2rad(beta1);
    
    rm2 = 0.5 * (rh(2) + rt(2));
    U2 = (2 * pi * rm2 * RPM) / 60;
    
    % Add a limiter to prevent extreme values
    % U2_Cm2_ratio = min(max(U2 / Cm2, -100), 100);  % Limit the ratio
    U2_Cm2_ratio = U2/Cm2;
    tan_alpha2 = U2_Cm2_ratio + tan(beta2_rad);

    alpha2_rad = atan(tan_alpha2);
    alpha2 = rad2deg(alpha2_rad);
    
    Ctheta1 = Cm1 * tan(alpha1_rad);
    Ctheta2 = Cm2 * tan(alpha2_rad);
    
    delta_h0 = (U2 * Ctheta2) - (U1 * Ctheta1);
    
    W1 = Wm1 / cos(beta1_rad);
    W2 = Wm2 / cos(beta2_rad);
    
    h2 = h1 + (W1.^2 / 2) - (U1^2 / 2) - (W2.^2 / 2) + (U2^2 / 2);
%     if ~isreal(h2), h2 = real(h2); end
%     
    h02_rel = h2 + (W2.^2 / 2);
%     if ~isreal(h02_rel), h02_rel = real(h02_rel); end
%     
    P02_rel = P01_rel - w_total_rotor * (P01_rel - P1);
%     if ~isreal(P02_rel), P02_rel = real(P02_rel); end
%     
    s2 = thermodynamic_calculator('S','H',h02_rel,'P',P02_rel);
%     if ~isreal(s2), s2 = real(s2); end
%     
    h02 = h01 + delta_h0;
%     if ~isreal(h02), h02 = real(h02); end
%     
    T02 = thermodynamic_calculator('T','H',h02,'S',s2);
%     if ~isreal(T02), T02 = real(T02); end
%     
    P02 = thermodynamic_calculator('P','T',T02,'S',s2);
%     if ~isreal(P02), P02 = real(P02); end
%     
    T2 = thermodynamic_calculator('T','H',h2,'S',s2);
%     if ~isreal(T2), T2 = real(T2); end
%     
    P2 = thermodynamic_calculator('P','T',T2,'S',s2);
%     if ~isreal(P2), P2 = real(P2); end
%     
%     % Ensure P2 is physically realistic
%     if P2 <= 0
%         fprintf('WARNING: Negative P2 detected, resetting to a small positive value\n');
%         P2 = 1000; % Set to a small positive pressure
%     end
%     
    rho2 = thermodynamic_calculator('D','T',T2,'P',P2);
%     if ~isreal(rho2), rho2 = real(rho2); end
%     
%     % Ensure rho2 is physically realistic
%     if rho2 <= 0
%         fprintf('WARNING: Negative rho2 detected, resetting to a small positive value\n');
%         rho2 = 0.01; % Set to a small positive density
%     end
    
    % Check for near-choked condition
    W_sound = sqrt(1.4 * 287 * T2); % Local speed of sound
    if W2 > 0.95 * W_sound
        warning('Flow approaching choked condition: W2/a = %f', W2/W_sound);
    end
    
    % --- Calculate new Cm2 ---
    H_out = rt(2) - rh(2);
    A2 = 2 * pi * rm2 * H_out;
    Cm2_new = m_dot / (rho2 * A2);
    
    % Ensure Cm2_new is real
    if ~isreal(Cm2_new)
        fprintf('WARNING: Complex Cm2_new value detected: real=%.6f, imag=%.6f\n', real(Cm2_new), imag(Cm2_new));
        Cm2_new = real(Cm2_new);
    end
    
%     % Ensure Cm2_new is physically realistic
%     if Cm2_new <= 0
%         fprintf('WARNING: Non-positive Cm2_new detected, using absolute value\n');
%         Cm2_new = abs(Cm2_new) + 0.1;
%     end
    
    % --- Check convergence and update Cm2 ---
    diff = abs(Cm2_new - Cm2);
    
%     % Output diagnostics
%     fprintf('\n--- Iteration %d ---\n', iter);
%     fprintf('Cm2 = %.6f m/s\n', Cm2);
%     fprintf('Cm2_new = %.6f m/s\n', Cm2_new);
%     fprintf('diff = %.6e\n', diff);
    
    if diff < tol
        converged = true;
        fprintf('Converged after %d iterations.\n', iter);
    else
        % Adaptive relaxation - use smaller values for difficult convergence
        if iter < 10
            relax = 0.5;  % More aggressive initially
        else
            relax = 0.2;  % More conservative later
        end
        
        Cm2 = (1-relax)*Cm2 + relax*Cm2_new;
    end
end

if ~converged
    warning('Did not converge after %d iterations. Final difference: %.6f', max_iter, diff);
end



% --- Final calculations with converged Cm2 ---
C2 = sqrt(Cm2^2 + (Cm2*tan(alpha2_rad))^2); % Absolute velocity at outlet

% --- Store results ---
rotor_exit_props = struct(...
    'Cm2', Cm2, ...
    'alpha2', alpha2, ...
    'beta2', beta2, ...
    'U2', U2, ...
    'rm', rm2, ...
    'P2', P2, ...
    'T2', T2, ...
    'rho2', rho2, ...
    'h2', h2, ...
    'P02', P02, ...
    'T02', T02, ...
    'h02', h02, ...
    'P02_rel', P02_rel, ...
    'h02_rel', h02_rel, ...
    's2', s2, ...
    'C2', C2, ...
    'iterations', iter, ...
    'converged', converged ...
);

% --- Display results ---
fprintf('\n--- Final Rotor Exit Properties ---\n');
fprintf('Outlet Axial Velocity (Cm2): %.4f m/s\n', Cm2);
fprintf('Outlet Absolute Flow Angle (alpha2): %.2f deg\n', alpha2);
fprintf('Outlet Static Pressure (P2): %.2f Pa\n', P2);
fprintf('Outlet Static Temperature (T2): %.2f K\n', T2);
fprintf('Outlet Density (rho2): %.4f kg/m^3\n', rho2);
fprintf('Outlet Absolute Velocity (C2): %.2f m/s\n', C2);
fprintf('Number of iterations: %d\n', iter);
fprintf('Converged: %s\n', string(converged));