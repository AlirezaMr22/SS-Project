function result = thermodynamic_calculator(output_prop, prop1, val1, prop2, val2)
    % thermodynamic_calculator - Calculate thermodynamic properties for air
    % Inputs:
    %   output_prop - Desired output property ('P', 'T', 'D', 'H', 'S', 'U', 'V')
    %   prop1       - First known property ('P', 'T', 'H', 'S', 'D')
    %   val1        - Value of the first known property
    %   prop2       - Second known property ('P', 'T', 'H', 'S', 'D')
    %   val2        - Value of the second known property
    % Output:
    %   result      - Calculated value of the desired property
    % Units:
    %   P: Pa, T: K, D: kg/m³, H: J/kg, S: J/kg·K, U: J/kg, V: m³/kg
    % Note: This function uses the Optimization Toolbox for fsolve.

    % Define gas properties for air mixture (hardcoded)
    gas(1).name = 'Oxygen';
    gas(1).coeffs_below_1000K = [0.362598E+01, -0.187854E-02, 0.705540E-05, -0.676351E-08, 0.215559E-11, -0.10475226E+04, 0.43052778E+01];
    gas(1).coeffs_above_1000K = [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00];
    gas(1).MW = 31.9988;
    gas(1).y = 0.2314212;  % Mass fraction of oxygen

    gas(2).name = 'Nitrogen';
    gas(2).coeffs_below_1000K = [0.36748261E+01, -0.12081500E-02, 0.23240102E-05, -0.63217559E-09, -0.22577253E-11, -0.10611588E+04, 0.23580424E+01];
    gas(2).coeffs_above_1000K = [0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04, 0.05980528E+02];
    gas(2).MW = 28.0134;
    gas(2).y = 0.7552179;  % Mass fraction of nitrogen

    gas(3).name = 'Argon';
    gas(3).coeffs_below_1000K = [0.02500000E+02, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, -0.07453750E+04, 0.04366000E+02];
    gas(3).coeffs_above_1000K = [0.02500000E+02, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, -0.07453750E+04, 0.04366000E+02];
    gas(3).MW = 39.948;
    gas(3).y = 0.01288;    % Mass fraction of argon

    gas(4).name = 'Carbon Dioxide';
    gas(4).coeffs_below_1000K = [2.35677352E+00, 8.98459677E-03, -7.12356269E-06, 2.45919022E-09, -1.43699548E-13, -4.83719697E+04, 9.90105222E+00];
    gas(4).coeffs_above_1000K = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00];
    gas(4).MW = 44.009;
    gas(4).y = 0.0004809;  % Mass fraction of carbon dioxide

    % Constants
    R = 8.314;       % Universal gas constant, J/(mol·K)
    P_ref = 100000;  % Reference pressure, Pa (1 bar)

    % Calculate mixture molar mass (g/mol)
    mass_fractions = [gas.y];
    molar_masses = [gas.MW];
    MW_mix = 1 / sum(mass_fractions ./ molar_masses); % g/mol

    % Conversion factor for molar to mass basis (mol/kg)
    conv = 1000 / MW_mix;

    % Calculate mole fractions
    for i = 1:length(gas)
        gas(i).x = (gas(i).y / gas(i).MW) / sum(mass_fractions ./ molar_masses);
    end

    % Parse input properties
    props = {prop1, prop2};
    vals = [val1, val2];

    % Find indices of provided properties
    idx_P = find(strcmp(props, 'P'));
    idx_T = find(strcmp(props, 'T'));
    idx_H = find(strcmp(props, 'H'));
    idx_S = find(strcmp(props, 'S'));
    idx_D = find(strcmp(props, 'D'));

    % Validate exactly two properties are provided
    if length(props) ~= 2
        error('Exactly two input properties must be provided.');
    end

    % Convert mass-based inputs to molar-based if necessary
    if ~isempty(idx_H)
        vals(idx_H) = vals(idx_H) / conv;  % Convert H from J/kg to J/mol
    end
    if ~isempty(idx_S)
        vals(idx_S) = vals(idx_S) / conv;  % Convert S from J/kg·K to J/mol·K
    end
    if ~isempty(idx_D)
        v_mass = 1 / vals(idx_D);  % m^3/kg
        v_molar = v_mass / conv;   % m^3/mol
    end

    % Calculate intermediate properties based on input pair
    if ~isempty(idx_P) && ~isempty(idx_T)
        P = vals(idx_P);
        T = vals(idx_T);
        [v_mass, h_mass, s_mass, u_mass, rho] = calc_PT(P, T, gas, R, P_ref, conv);
    elseif ~isempty(idx_P) && ~isempty(idx_H)
        P = vals(idx_P);
        h_molar = vals(idx_H);
        [T, v_mass, s_mass, u_mass, rho] = calc_PH(P, h_molar, gas, R, P_ref, conv);
    elseif ~isempty(idx_P) && ~isempty(idx_S)
        P = vals(idx_P);
        s_molar = vals(idx_S);
        [T, v_mass, h_mass, u_mass, rho] = calc_PS(P, s_molar, gas, R, P_ref, conv);
    elseif ~isempty(idx_T) && ~isempty(idx_S)
        T = vals(idx_T);
        s_molar = vals(idx_S);
        [P, v_mass, h_mass, u_mass, rho] = calc_TS(T, s_molar, gas, R, P_ref, conv);
    elseif ~isempty(idx_H) && ~isempty(idx_S)
        h_molar = vals(idx_H);
        s_molar = vals(idx_S);
        [T, P, v_mass, u_mass, rho] = calc_HS(h_molar, s_molar, gas, R, P_ref, conv);
    else
        error('Unsupported property combination.');
    end

    % Return the desired output property
    switch upper(output_prop)
        case 'P'
            result = P;
        case 'T'
            result = T;
        case 'D'
            result = rho;
        case 'H'
            if exist('h_mass', 'var')
                result = h_mass;
            else
                h_molar = calc_enthalpy(T, gas);
                result = h_molar * conv;  % J/kg
            end
        case 'S'
            if exist('s_mass', 'var')
                result = s_mass;
            else
                s_molar = calc_entropy(T, P, gas, R, P_ref);
                result = s_molar * conv;  % J/kg·K
            end
        case 'U'
            if exist('u_mass', 'var')
                result = u_mass;
            else
                v_molar = R * T / P;
                h_molar = calc_enthalpy(T, gas);
                u_molar = h_molar - P * v_molar;
                result = u_molar * conv;  % J/kg
            end
        case 'V'
            if exist('v_mass', 'var')
                result = v_mass;
            else
                v_molar = R * T / P;
                result = v_molar * conv;  % m^3/kg
            end
        default
            error('Unsupported output property: %s', output_prop);
    end
end

% Sub-function: Calculate properties from P and T
function [v_mass, h_mass, s_mass, u_mass, rho] = calc_PT(P, T, gas, R, P_ref, conv)
    v_molar = R * T / P;  % m^3/mol
    [h_molar, s_molar] = calc_thermodynamic_properties(T, P, gas, R, P_ref);
    u_molar = h_molar - P * v_molar;
    v_mass = v_molar * conv;  % m^3/kg
    h_mass = h_molar * conv;  % J/kg
    s_mass = s_molar * conv;  % J/kg·K
    u_mass = u_molar * conv;  % J/kg
    rho = 1 / v_mass;         % kg/m^3
end

% Sub-function: Calculate properties from P and H
function [T, v_mass, s_mass, u_mass, rho] = calc_PH(P, h_molar, gas, R, P_ref, conv)
    T_guess = 300;  % Initial guess, K
    options = optimset('Display', 'off');
    f = @(T) calc_enthalpy(T, gas) - h_molar;
    T = fsolve(f, T_guess, options);
    v_molar = R * T / P;
    [~, s_molar] = calc_thermodynamic_properties(T, P, gas, R, P_ref);
    u_molar = h_molar - P * v_molar;
    v_mass = v_molar * conv;
    s_mass = s_molar * conv;
    u_mass = u_molar * conv;
    rho = 1 / v_mass;
end

% Sub-function: Calculate properties from P and S
function [T, v_mass, h_mass, u_mass, rho] = calc_PS(P, s_molar, gas, R, P_ref, conv)
    T_guess = 300;
    options = optimset('Display', 'off');
    f = @(T) calc_entropy(T, P, gas, R, P_ref) - s_molar;
    T = fsolve(f, T_guess, options);
    v_molar = R * T / P;
    h_molar = calc_enthalpy(T, gas);
    u_molar = h_molar - P * v_molar;
    v_mass = v_molar * conv;
    h_mass = h_molar * conv;
    u_mass = u_molar * conv;
    rho = 1 / v_mass;
end

% Sub-function: Calculate properties from T and S
function [P, v_mass, h_mass, u_mass, rho] = calc_TS(T, s_molar, gas, R, P_ref, conv)
    s_0_molar = calc_entropy(T, P_ref, gas, R, P_ref);
    P = P_ref * exp((s_0_molar - s_molar) / R);
    v_molar = R * T / P;
    h_molar = calc_enthalpy(T, gas);
    u_molar = h_molar - P * v_molar;
    v_mass = v_molar * conv;
    h_mass = h_molar * conv;
    u_mass = u_molar * conv;
    rho = 1 / v_mass;
end

% Sub-function: Calculate properties from H and S
function [T, P, v_mass, u_mass, rho] = calc_HS(h_molar, s_molar, gas, R, P_ref, conv)
    T_guess = 300;
    options = optimset('Display', 'off');
    f_T = @(T) calc_enthalpy(T, gas) - h_molar;
    T = fsolve(f_T, T_guess, options);
    s_0_molar = calc_entropy(T, P_ref, gas, R, P_ref);
    P = P_ref * exp((s_0_molar - s_molar) / R);
    v_molar = R * T / P;
    u_molar = h_molar - P * v_molar;
    v_mass = v_molar * conv;
    u_mass = u_molar * conv;
    rho = 1 / v_mass;
end

% Helper function: Calculate enthalpy for the air mixture
function h_molar = calc_enthalpy(T, gas)
    h_mix = 0;
    for i = 1:length(gas)
        if T < 1000
            coeffs = gas(i).coeffs_below_1000K;
        else
            coeffs = gas(i).coeffs_above_1000K;
        end
        a = coeffs;
        % NASA polynomial for enthalpy: h / R = a1*T + a2*T^2/2 + ... + a6
        h_i = 8.314 * (a(1)*T + a(2)*T^2/2 + a(3)*T^3/3 + a(4)*T^4/4 + a(5)*T^5/5 + a(6));
        h_mix = h_mix + gas(i).x * h_i;
    end
    h_molar = h_mix;
end

% Helper function: Calculate entropy for the air mixture
function s_molar = calc_entropy(T, P, gas, R, P_ref)
    s_mix = 0;
    for i = 1:length(gas)
        if T < 1000
            coeffs = gas(i).coeffs_below_1000K;
        else
            coeffs = gas(i).coeffs_above_1000K;
        end
        a = coeffs;
        % NASA polynomial for entropy: s / R = a1*ln(T) + a2*T + ... + a7
        s0_i = 8.314 * (a(1)*log(T) + a(2)*T + a(3)*T^2/2 + a(4)*T^3/3 + a(5)*T^4/4 + a(7));
        s_i = s0_i - R * log(gas(i).x * P / P_ref);
        s_mix = s_mix + gas(i).x * s_i;
    end
    s_molar = s_mix;
end

% Helper function: Calculate both enthalpy and entropy
function [h_molar, s_molar] = calc_thermodynamic_properties(T, P, gas, R, P_ref)
    h_molar = calc_enthalpy(T, gas);
    s_molar = calc_entropy(T, P, gas, R, P_ref);
end