function delta_star = deviation( tb_c, sigma, beta1, betap2, theta)

% --- Define Input Variables ---
% You MUST replace these example values with your actual values for the
% specific stage and profile you are analyzing.
% profile_type = 'NACA65'; % Options: 'NACA65', 'C4', 'DCA' (Case-insensitive)
% tb_c = 0.058;           % Example: Max thickness-to-chord ratio (t_b/c)
% sigma = 1.45;           % Example: Solidity (c/s)
% beta1 = 27.83;        % Example: Inlet flow angle  in DEGREES
% theta = 30;        % Example: Camber angle  in DEGREES
%betap2

%'NACA65'

Ksh = 1.0;
    
% --- Calculate Kt,delta (Thickness correction factor) ---
Kt_delta = (6.25 * tb_c) + (37.5 * (tb_c^2));


term1_delta010 = 0.01 * sigma * beta1;
term2_delta010_base = 0.74 * sigma^1.9 + 3 * sigma;
term3_delta010_power = (beta1 / 90)^(1.67 + 1.09 * sigma); 
delta_0_10 = term1_delta010 + term2_delta010_base * term3_delta010_power;
% 
% x = beta1 / 100;
% 
% 
% m10 = 0.17 - 0.0333*x + 0.333*x^2;
% 
% 
% b = 0.9625 - 0.17*x - 0.85*x^3;
% 
% m = m10 / (sigma^b);
% 
% 
% 
% delta_star = Ksh * Kt_delta * delta_0_10 + (m * theta);
%%
a_c = 0.4;
delta_star = ((0.92*((a_c)^2) + 0.002*betap2)*theta) / (sqrt(sigma)-0.002*theta)...
    +(Ksh*Kt_delta -1)*delta_0_10 ;



% --- Display the Results ---
fprintf('\n--- Final Result ---\n');
fprintf('delta* (Design Deviation Angle): %f degrees\n', delta_star);
