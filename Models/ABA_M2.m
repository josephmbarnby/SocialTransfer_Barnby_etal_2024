% CC JMBarnby (2024)
%
% josephmbarnby@gmail.com

% FULL BAYESIAN GENERATIVE MODEL
% Feed the model the desired parameters when simulating (e.g. alpha = 2,
% beta = -5)
% Feed the model normally distributed parameters when fitting
% e.g. m = 0, v = 7.5

%%%%%% DATA STRUCTURE REQUIRED %%%%%%%
% colnames:
% 1 ID                    | 1....n
% 2 Trial                 | 1....n (max trials of T1 + T2 + T3)
% 3 O1-Self               | 8, 10, 10...n
% 4 O1-Other              | 8, 6,  5 ...n
% 5 O2-Self               | 6, 9,  8 ...n
% 6 O2-Other              | 2, 9,  4 ...n
% 7 PPT choice/prediction | 2, 2,  1 ...n % if simulating this column is irrelevant
% 8 Partner action        | 1, 2,  1 ...n
% 9 Phase                 | 1, 1,  1 ...n

%% Model

function[F, results] = ABA_insert_noconv_Gen(parms, data, sim)

% Check if the sim argument is provided
    if nargin < 3
        sim = 0;  % Set default value for sim
    end

% Validate input parameters
   if length(parms) ~= 6
       error('Parameter vector should have 6 elements.');
   end
   if size(data, 2) ~= 9
       error('Data should have 9 columns.');
   end

   % Initialise
   res = 30;       % resolution of belief grid
   inc = 0.25;     % res of inc grid
   T1  = 36;       % length(find(data(:,9) == 1));  % trials for phase 1
   T2  = 54 + T1;  % length(find(data(:,9) == 2)) + T1;  % trials for phase 2
   T3  = T1 + T2;  % trials for phase 3

   %phase 1 parms
   alpha_raw       = parms(1); % subjects alpha for phase 1
   alpha           = res*(1./(1+exp(-alpha_raw))); % restrict between 0 and res
   beta            = parms(2); % subjects beta for phase 1
   alpha_v          = exp(parms(3)); % restrict the variance to above 0
   beta_v           = exp(parms(4));
   
       %phase 2 parms
   alpha_v_ref      = exp(parms(5)); % restrict the variance to above 0
   beta_v_ref       = exp(parms(6));

   %parameters space of alpha and beta
   [alpha_grid,beta_grid]=meshgrid(0:(inc/2):res,-res:inc:res);
   
    %generate standardised grid to form priors for preferences
   pabg=...
       (normpdf(alpha_grid,alpha,alpha_v).*...
       normpdf(beta_grid,beta,beta_v))+eps;
   
   pabg=pabg ./ sum(pabg(:));
   
   alpha_marg1 = squeeze(sum(pabg,[1 3]));
   beta_marg1  = squeeze(sum(pabg,[2 3]))';

    % grid for a subjects beliefs over their partner in phase 2

   pabg_par=...
       (normpdf(alpha_grid,alpha,alpha_v_ref).*...
       normpdf(beta_grid,beta,beta_v_ref))+eps;
   
   pabg_par=pabg_par ./ sum(pabg_par(:));
   
   alpha_marg2a = squeeze(sum(pabg_par,[1 3]));
   beta_marg2a  = squeeze(sum(pabg_par,[2 3]))';
   
       % initialised dummy values
   
   lik1        = 0;   % likelihood for choices in phase 1
   prob1       = zeros(T1, 1);
   lik2        = 0;   % likelihood for guesses in phase 2
   prob2       = zeros(T2-T1, 1);
   lik3        = 0;
   prob3       = zeros(T1, 1);
   
   simA        = zeros(T3,1);
   simAFix     = zeros(T3,1);
   rew         = zeros(T3,4);
   
    % Phase 1 choices of the participant

    for t=1:T1

    %Extracts vals
    [s1, o1, s2, o2] = extract_normalized_data(data, t);

    %Save Rewards
    rew(t,:) = [s1, o1, s2, o2]; %save rewards

    %Calculate Utilities
    [val1, val2] = calculate_utilities(alpha_grid, beta_grid, s1, o1, s2, o2);

    %Sim P(choice==1)
    subject_netp1 = calculate_subject_netp(val1, val2, pabg);

    %Sim Action
    simA(t) = simulate_action(subject_netp1);

    %Participant estimate
    subject_estimate = data(t, 7); % subject choice

        if sim==1 
          subject_estimate = simA(t); 
        end    
    
    %Likelihood of choice
    [lik1,prob1(t)] = update_likelihood(lik1, subject_estimate, subject_netp1);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % Phase2

    alpha_cont = zeros(length(alpha_grid), length((T1+1):T2));
    beta_cont  = zeros(length(beta_grid), length((T1+1):T2));

    for t=(T1+1):T2

    alpha_cont(:,t) = squeeze(sum(pabg_par,[1 3]));
    beta_cont(:,t)  = squeeze(sum(pabg_par,[2 3]))';

    %Extracts vals
    [s1, o1, s2, o2] = extract_normalized_data(data, t);

    %Save Rewards
    rew(t,:) = [s1, o1, s2, o2]; %save rewards

    %Calculate Utilities
    [val1, val2] = calculate_utilities(alpha_grid, beta_grid, s1, o1, s2, o2);

    %Sim P(choice==1)
    subject_netp2     = calculate_subject_netp(val1, val2, pabg_par);
    subject_netp2_fix = calculate_subject_netp(val1, val2, pabg);

    %Sim Action
    simA(t) = simulate_action(subject_netp2);
    simAFix(t) = simulate_action(subject_netp2_fix);

    %Participant estimate
    subject_estimate = data(t, 7); % subject choice

        if sim==1 
          subject_estimate = simA(t); 
        end    
    
    %Likelihood of choice
    [lik2,prob2(t)] = update_likelihood(lik2, subject_estimate, subject_netp2);

    %Partner decision
    actual_choice = data(t, 8); % what did our partner 'choose'

    %Belief update
    pabg_par = update_beliefs(val1, val2, actual_choice, pabg_par);

    end

    alpha_marg2b = squeeze(sum(pabg_par,[1 3]));
    beta_marg2b  = squeeze(sum(pabg_par,[2 3]))';

    % Phase 3 choices of the participant

    for t=(T2+1):T3

    %Extracts vals
    [s1, o1, s2, o2] = extract_normalized_data(data, t);

    %Save Rewards
    rew(t,:) = [s1, o1, s2, o2]; %save rewards

    %Calculate Utilities
    [val1, val2] = calculate_utilities(alpha_grid, beta_grid, s1, o1, s2, o2);

    %Sim P(choice==1)
    subject_netp3 = calculate_subject_netp(val1, val2, pabg);

    %Sim Action
    simA(t) = simulate_action(subject_netp3);

    %Participant estimate
    subject_estimate = data(t, 7); % subject choice

        if sim==1 
          subject_estimate = simA(t); 
        end    
    
    %Likelihood of choice
    [lik3,prob3(t)] = update_likelihood(lik3, subject_estimate, subject_netp3);

    end

    alpha_marg3 = squeeze(sum(pabg,[1 3]));
    beta_marg3  = squeeze(sum(pabg,[2 3]))';

    F  = lik1 + lik2 + lik3 + eps;

    results = struct;
    results.lik1 = lik1;
    results.lik2 = lik2;
    results.lik3 = lik3;
    results.simA = simA;
    results.simAFix = simAFix;
    results.prob1 = prob1;
    results.prob2 = prob2;
    results.prob3 = prob3;
    results.rew = rew;
    results.val1 = val1;
    results.val2 = val2;
    results.alpha_marg1  = alpha_marg1;
    results.beta_marg1   = beta_marg1;
    results.alpha_marg2a = alpha_marg2a;
    results.beta_marg2a  = beta_marg2a;
    results.alpha_cont   = alpha_cont;
    results.beta_cont    = beta_cont;    
    results.alpha_marg2b = alpha_marg2b;
    results.beta_marg2b  = beta_marg2b;
    results.alpha_marg3  = alpha_marg3 ;
    results.beta_marg3   = beta_marg3 ;

end

