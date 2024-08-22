function pabg_par = update_beliefs(val1, val2, actual_choice, pabg_par)
    if (actual_choice == 1)
        pchoose2 = (1 ./ (1 + exp(-(val1 - val2)))); % probability of 1
    else
        pchoose2 = (1 ./ (1 + exp(-(val2 - val1)))); % probability of 2
    end
    pabg_par = pchoose2 .* pabg_par; % Apply Bayes' rule
    pabg_par = pabg_par ./ sum(pabg_par(:)); % Normalize
end