function [istar , ic , is, istar_low_mach] = incidence_calc (tb_c , sigma , alpha1, rh, rt, RPM, T01, betap1,theta)


Ksh = 1;  % Shape factor Ksh

% Initialize variables
beta1_guess = betap1;      % Initial guess for Beta1
tol = 0.001;           % Convergence tolerance
max_iter = 100;        % Maximum iterations
converged = false;
iter = 0;

while ~converged && iter < max_iter
    iter = iter + 1;
    beta1_prev = beta1_guess;
    
    q = 0.28 / (0.1 + (tb_c^0.3));

    Kt_i = (10 * tb_c)^q;

    p = 0.914 + (sigma^3) / 160;

    i_0_10_term1 = (beta1_guess^p) / (5 + 46 * exp(-2.3 * sigma));
    i_0_10_term2 = 0.1 * (sigma^3) * (exp((beta1_guess - 70) / 4));
    i_0_10 = i_0_10_term1 - i_0_10_term2;

    n_numerator = (beta1_guess / 90)^(1 + 1.2 * sigma);
    n_denominator = 1.5 + 0.43 * sigma;
    n_term = n_numerator / n_denominator;
    n = 0.025 * sigma - 0.06 - n_term;

    istar_low_mach = Ksh * Kt_i * i_0_10 + n * theta;

    % Update beta1 guess
    beta1_guess = istar_low_mach + betap1;
    
    % Check convergence
    if abs(beta1_guess - beta1_prev) < tol
        converged = true;
    end
end



% Display results
if converged
    fprintf('istar_minimumloss Converged in %d iterations:\n', iter);
else
    fprintf('Stopped after %d iterations (did not fully converge):\n', iter);
end

fprintf('istar_minimumloss = %.2f\n',istar_low_mach);
%%
% % --- Define Input Variables ---
% beta1c_guess = istar_low_mach - 8;   % Initial guess for inlet relative flow angle at choke
% beta1s_guess = istar_low_mach + 8;   % Initial guess for inlet relative flow angle at stall
% 
% % Convergence tolerance
% tol = 0.001;
% max_iter = 100;
% converged = false;
% iter = 0;
% 
% while ~converged && iter < max_iter
%     iter = iter + 1;
%     
%     % Store previous values for convergence check
%     beta1c_prev = beta1c_guess;
%     beta1s_prev = beta1s_guess;
%     
%     term1_ic =  1 - real(30 ./ beta1c_guess)^0.48;
%     term2_ic = theta / 4.176;
%     ic = istar_low_mach - 9 + term1_ic * term2_ic;
%     
%     term1_is =  2.92 - real(beta1s_guess / 15.6);
%     term2_is = theta / 8.2;
%     is = istar_low_mach + 10.3 + term1_is * term2_is;
%     
%     % Update beta1 guesses
%     beta1c_guess = real(ic) + betap1;
%     beta1s_guess = is + betap1;
%     
%     % Check convergence
%     if abs(beta1c_guess - beta1c_prev) < tol && abs(beta1s_guess - beta1s_prev) < tol
%         converged = true;
%     end
% end
% 
% if ~converged
%     warning('Failed to converge within %d iterations', max_iter);
% end
% 
% % Final values after convergence
% beta1c = real(beta1c_guess);
% beta1s = real(beta1s_guess);
% 
% ic = istar_low_mach - 9 + (1 - (30 ./ beta1c)^0.48) * (theta ./ 4.176);
% ic = real(ic);
% 
% is = istar_low_mach + 10.3 + (2.92 - (beta1s / 15.6)) * (theta / 8.2);
% 
% 
% fprintf('ic = %.2f\n',ic);
% fprintf('is = %.2f\n',is);
%% choke and stall incidence

term1_Rc = 10 - theta * ( (betap1 - 40) / 450 );
term2_Rc = 0.5 + 5 * tb_c;
Rc = term1_Rc * term2_Rc;
ic = istar_low_mach - Rc;

term1_Rs = 10 + theta * ( (55 - betap1) / 150 );
term2_Rs = 0.5 + 5 * tb_c;
Rs = 1.5 * term1_Rs * term2_Rs; 
is = istar_low_mach + Rs;

% --- Display the Results ---
fprintf('Calculated Choke Incidence Angle (ic1): %f degrees\n', ic);
fprintf('Calculated Stall Incidence Angle (is1): %f degrees\n', is);

%%   Mach Corrected Incidence Angles

beta1 = betap1 + istar_low_mach;
beta1_rad = deg2rad(beta1);
alpha1_rad = deg2rad(alpha1);

rm = 0.5 * (rh(1) + rt(1));
omega_rad_s = RPM * (2 * pi / 60); % Rotational speed in rad/s
U1 = omega_rad_s * rm;
Cm1 = U1 / (-tan(alpha1_rad) + tan(beta1_rad));
C1= Cm1 / cos (alpha1_rad);
T1 = T01 - ((0.5/1005)*C1^2);
W1 = Cm1 / cos (beta1_rad);
MW1= W1 / (sqrt(1.4 * 287 * T1 ));

%

Rc = istar_low_mach - ic;
Rs= is - istar_low_mach;
% Equation for Mach-corrected Choke Incidence (ic)
ic = istar_low_mach - (Rc ./ (1 + 0.5 * (MW1.^3)));

% Equation for Mach-corrected Stall Incidence (is)
is = istar_low_mach + (Rs ./ (1 + 0.5 * (Ksh .* MW1).^3));

% Equation for Mach-corrected 'Mid-range' Incidence (im)
istar = ic + ((is - ic) .* Rc) ./ (Rc + Rs);



fprintf('\nMach-corrected Choke Incidence (ic): %.2f deg\n', ic);
fprintf('Mach-corrected Stall Incidence (is): %.2f deg\n', is);
fprintf('Mach-corrected Mid-range Incidence (im): %.2f deg\n', istar);

fprintf('--- Mach Incidence Correction  Complete ---\n');

