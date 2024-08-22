function [val1, val2] = calculate_utilities(alpha_grid, beta_grid, s1, o1, s2, o2)
% Check if alpha_grid is provided as a non-empty argument
    if nargin < 6 || isempty(alpha_grid)
        % If alpha_grid is not provided, set it to a default value (e.g., 0)
        alpha_grid = 1;
    end

    % Calculate utilities
    val1 = (alpha_grid * s1) + (beta_grid * max(s1 - o1, 0));
    val2 = (alpha_grid * s2) + (beta_grid * max(s2 - o2, 0));
    
end