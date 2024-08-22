function subject_netp = calculate_subject_netp(val1, val2, belief)
    subject_estimate_pchoose1 = (1 ./ (1 + exp(-(val1 - val2))));
    tmp_par = subject_estimate_pchoose1 .* belief;
    subject_netp = sum(tmp_par(:));
end