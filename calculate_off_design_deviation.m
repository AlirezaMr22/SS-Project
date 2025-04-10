function [delta] = calculate_off_design_deviation(delta_star, istar, incidence, Cm1, sigma, beta1)

    numerator_slope = 1 + ((sigma + 0.25 * sigma^4) * ((beta1 / 53)^2.5));
    denominator_slope = exp(3.1 * sigma);

    slope_d_delta_d_i = numerator_slope / denominator_slope;

    incidence_diff = incidence - istar;
    
    Wm1=Cm1;
    Wm2=Cm1;
    velocity_correction = 10 * (1 - (Wm2 / Wm1));

    delta = delta_star + (slope_d_delta_d_i * incidence_diff) + velocity_correction;
    fprintf('  Calculated off-design deviation  = %.2f degrees\n', delta);

