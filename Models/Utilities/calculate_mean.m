function mu = calculate_mean(grid, marg)    
    % Calculate the mean
    mu = sum(grid .* marg);
end