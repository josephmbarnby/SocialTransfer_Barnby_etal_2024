function simA_t = simulate_action(subject_netp)
    simA_t = randsample([1, 2], 1, true, [subject_netp, 1 - subject_netp]);
end