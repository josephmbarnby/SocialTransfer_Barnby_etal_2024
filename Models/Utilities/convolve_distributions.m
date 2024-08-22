function [convolve_self_v, convolve_self_m] = convolve_distributions(self, self_v, self_v_ref, sigma_par, mu_par)

    % Precompute frequently used terms
    self_v_dash          = self_v^2;
    self_v_ref_dash      = self_v_ref^2;
    sigma_par_dash       = sigma_par^2;
    
    % Compute convolve_self_v
    self_var             = self_v_dash^-1 + (2 * self_v_ref_dash + sigma_par_dash)^-1;
    convolve_self_v_dash = self_var^-1;
    convolve_self_v      = sqrt(self_var^-1);
    
    % Compute convolve_self_m
    convolve_self_m = convolve_self_v_dash * (self_v_dash^-1 * self...
        + (2*self_v_ref_dash + sigma_par_dash)^-1 * mu_par);

end