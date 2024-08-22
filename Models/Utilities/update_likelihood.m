function [newlik, prob_t] = update_likelihood(oldlik, subject_estimate, subject_netp)
    if (subject_estimate == 1)
        newlik = oldlik + log(subject_netp); % log likelihood
        prob_t = subject_netp;
    else
        newlik = oldlik + log(1 - subject_netp);
        prob_t = 1 - subject_netp;
    end
end