function [sigma_par] = calculate_std(grid, marg, mu)
    marg = marg/sum(marg(:));
    var_par = sum((grid - mu).^2 .* marg);
    sigma_par = sqrt(var_par);
end