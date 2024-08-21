#Utility script for:
#Transfer of social information reduces uncertainty about the self and others
#CC Joe Barnby 2024

# Full Generative Model ---------------------------------------------------

# Feed the model the desired parameters when simulating (e.g. alpha = 2,
                                                         # beta = -5)
# Feed the model normally distributed parameters when fitting
# e.g. m = 0, v = 7.5

# Estimates:

#1  alpha_mu; alpha_sigma and beta_mu; beta_sigma of
#   participant over trials 1:n in phase 1

#2  The variance of the alpha and beta priors a participant holds about
#   their partner

#3  A participant's inferred beliefs about the partner given their
#   partner's decisions over trials 1:36.

#4  A participant's alpha_mu;alpha_sigma and beta_mu; beta_sigma given
#   their partner's decisions over another 18 trials

###### DATA STRUCTURE REQUIRED ###
# colnames:
  # 1 ID                  | 1....n
# 2 Trial                 | 1....n (max trials of T1 + T2 + T3)
# 3 O1-Self               | 8, 10, 10...n
# 4 O1-Other              | 8, 6,  5 ...n
# 5 O2-Self               | 6, 9,  8 ...n
# 6 O2-Other              | 2, 9,  4 ...n
# 7 PPT choice/prediction | 2, 2,  1 ...n # if simulating this column is irrelevant
# 8 Partner action        | 1, 2,  1 ...n

## Model

ABA_shift_Gen <- function(parms, df, sim = 0, plot = 0, phase = 3, conv = 1){

  #Arguments

  #Necessary
  # parms = parameter vector (either 4 or 6 parameters; see conv)
  # df    = dataframe of structure as listed above

  #Optional
  # sim   = do you want the function to simlate data (1=yes, 0=no)
  # plot  = do you want a plot at the end (1=yes, 0=no)
  # phase = how many phases of the game (can only be 1 2 or 3)
  # conv  = is the third phase a convolution of the first and second phase? If so
  #         the function requires 2 extra parameters to define the variance of the
  #         participants beliefs in the first phase to convolve with the partners

  # Initialise
  res = 30; # bounds of belief grid
  inc = 0.125; # resolution of grid
  eps = 0.02/(length(seq(-res, res, inc))^2) #noise floor for

  if (phase == 1) {
    data <- df[df$Phase == 1, ]
  } else if (phase == 2) {
    data <- df[df$Phase %in% c(1, 2), ]
  } else if (phase == 3) {
    data <- df[df$Phase %in% c(1, 2, 3), ]
  }

  T1 <- sum(data$Phase == 1)
  T2 <- if (phase != 1) sum(data$Phase == 2) else 0
  T3 <- if (phase == 3) sum(data$Phase == 3) else 0

  if(phase > 3){stop('This function currently only supports (up to) three-phase estimation')}

  #Phase 1
  alpha <- ifelse(sim, parms[1], (1 / (1 + exp(-parms[1]))) * res)
  beta  <- parms[2]

  if (!sim) {
    parmsP <- c(alpha, beta)
    if (conv) {
      alpha_v <- exp(parms[3])
      beta_v  <- exp(parms[4])
      parmsP  <- c(parmsP, alpha_v, beta_v)
    }
  }
  if (sim) {
    parmsP <- c(alpha, beta)
    if (conv) {
      alpha_v <- parms[3]
      beta_v  <- parms[4]
      parmsP  <- c(parmsP, alpha_v, beta_v)
    }
  }

  if(phase != 1){
  #Phase 2
    if(!sim){
      alpha_v_ref      = ifelse(conv, exp(parms[5]), exp(parms[3]))
      beta_v_ref       = ifelse(conv, exp(parms[6]), exp(parms[4]))
      parmsP = c(parmsP, alpha_v_ref, beta_v_ref)
    }
    if(sim){
      alpha_v_ref      = ifelse(conv, parms[5], parms[3])
      beta_v_ref       = ifelse(conv, parms[6], parms[4])
      parmsP = c(parmsP, alpha_v_ref, beta_v_ref)
    }
  }

  #output list
  simulatedDF <- create_output_file(conv=conv, res=res, data=data, inc=inc, T1, T2, T3, parmsP)

  simulatedDF$Priors[1:length(parms)] <- parmsP[1:length(parms)]

  # grid for a subjects beliefs over their preferences in phase 1

  #parameters space of alpha and beta
  grid <- meshgrid(seq(res-res, res, inc),seq(-res, res, inc*2));
  alpha_grid <- grid$x
  beta_grid  <- grid$y

  if(conv){
  #generate standardised grid to form priors for preferences
  pabg <- (dnorm(alpha_grid,alpha,alpha_v)+eps)*(dnorm(beta_grid,beta,beta_v)+eps);
  pabg <- pabg / sum(as.vector(pabg));
  simulatedDF$marginals$Phase1_alpha = colSums(pabg);
  simulatedDF$marginals$Phase1_beta  = rowSums(pabg);
  }

  # initialised dummy values

  lik1        <- 0;   # likelihood for choices in phase 1
  prob1       <- rep(NA, T1)
  lik2        <- 0;   # likelihood for guesses in phase 2
  prob2       <- rep(NA, T2)
  lik3        <- 0;
  prob3       <- rep(NA, T3)
  LL          <- 0

  simA        <- rep(NA, (T1+T2+T3))
  simAFix     <- simA
  correct     <- simA
  rew         <- matrix(NA, nrow = (T1+T2+T3), ncol = 4)

  # Phase 1 choices of the participant

    for (t in 1:T1){

    s1 = as.numeric(data[t, 3]/10);
    o1 = as.numeric(data[t, 4]/10);
    s2 = as.numeric(data[t, 5]/10);
    o2 = as.numeric(data[t, 6]/10);

    rew[t,1] = s1
    rew[t,2] = o1
    rew[t,3] = s2
    rew[t,4] = o2

    if(conv){
    val1 = alpha_grid*s1 + beta_grid*max(s1-o1,0) ;
    val2 = alpha_grid*s2 + beta_grid*max(s2-o2,0) ;
    } else {
    val1 = alpha*s1 + beta*max(s1-o1,0) ;
    val2 = alpha*s2 + beta*max(s2-o2,0) ;
    }

    pchoose1=(1/(1+exp(-(val1 - val2)))); # probability of 1

    tmp_ppt = ifelse(conv, (pchoose1*pabg), pchoose1);

    subject_netp1 = sum(as.vector(tmp_ppt));
    simA[t] = sample(c(1,2),1, prob = c(subject_netp1, 1-subject_netp1));

      if (sim){
        actual_choice = simA[t];
      } else {
        actual_choice = data[t, 7];
      }

      if (actual_choice==1){
      lik1 = lik1 + log(subject_netp1); # log likelihood of 1
      prob1[t] = subject_netp1;
      } else {
      lik1 = lik1 + log(1-subject_netp1);
      prob1[t] = 1-subject_netp1;
      }

    }

####

  # how the experimenter learns how the subject learns online about the partner in phase 2

  #generate standardised grid to form priors for beliefs
  pabg_par <- (dnorm(alpha_grid,alpha,alpha_v_ref)+eps)*(dnorm(beta_grid,beta,beta_v_ref)+eps);
  pabg_par <- pabg_par/sum(as.vector(pabg_par));

  simulatedDF$marginals$Phase2a_alpha = colSums(pabg_par);
  simulatedDF$marginals$Phase2a_beta  = rowSums(pabg_par);

    # Phase2

    if(phase != 1){

    for (t in (T1+1):(T2+T1)){

    s1 = as.numeric(data[t, 3]/10);
    o1 = as.numeric(data[t, 4]/10);
    s2 = as.numeric(data[t, 5]/10);
    o2 = as.numeric(data[t, 6]/10);

    rew[t,1] = s1
    rew[t,2] = o1
    rew[t,3] = s2
    rew[t,4] = o2

    val1 = alpha_grid*s1 + beta_grid*max(s1-o1,0) ;
    val2 = alpha_grid*s2 + beta_grid*max(s2-o2,0) ;

    subject_estimate_pchoose1 = (1/(1+exp(-(val1 - val2))));

    tmp_par = subject_estimate_pchoose1 * pabg_par;
    subject_netp2 = sum(as.vector(tmp_par));

    if(conv){
      tmp_fix = subject_estimate_pchoose1 * pabg;
      subject_netp2_fix = sum(as.vector(tmp_fix));
    } else {
      val1 = alpha*s1 + beta*max(s1-o1,0) ;
      val2 = alpha*s2 + beta*max(s2-o2,0) ;
      subject_netp2_fix = (1/(1+exp(-(val1 - val2))))
    }

    simA[t]    = sample(c(1,2),1,prob = c(subject_netp2, 1-subject_netp2));
    simAFix[t] = sample(c(1,2),1,prob = c(subject_netp2_fix, 1-subject_netp2_fix));

      if (sim){
        subject_estimate = simA[t];
      } else {
        subject_estimate = data[t, 7]; # subject choice
      }

      if (subject_estimate == 1){
        lik2        = lik2 + log(subject_netp2); # log likelihood
        prob2[t-T1] = subject_netp2;
      } else {
        lik2        = lik2 + log(1-subject_netp2);
        prob2[t-T1] = 1-subject_netp2;
      }

    actual_choice = data[t, 8]; # what did our partner 'choose'
    if(actual_choice == subject_estimate){correct[t] <- 1} else {correct[t] <- 0}

      if (actual_choice==1){
        pchoose2=(1./(1+exp(-(val1 - val2)))); # probability of 1
      } else {
        pchoose2=(1./(1+exp(-(val2 - val1)))); # probability of 2
      }

    pabg_par = pchoose2*pabg_par; # Bayes rule
    pabg_par = pabg_par / sum(as.vector(pabg_par)); #normalised distribution

    }

    simulatedDF$marginals$Phase2_alpha = colSums(pabg_par);
    simulatedDF$marginals$Phase2_beta  = rowSums(pabg_par);

    #find the indices of maximum likelihood after seeing the partner
    mu_alphapar   = sum(alpha_grid[1, ] * simulatedDF$marginals$Phase2_alpha)
    mu_betapar    = sum(beta_grid[, 1] * simulatedDF$marginals$Phase2_beta)
    sddevs  = std_dev_prob(alpha_grid, beta_grid,
                           mu_alphapar, mu_betapar,
                           simulatedDF$marginals$Phase2_alpha,
                           simulatedDF$marginals$Phase2_beta)
    sigma_alphapar= sddevs[1]
    sigma_betapar = sddevs[2]

    simulatedDF$PartnerParms[1:4] <- c(mu_alphapar, mu_betapar, sigma_alphapar, sigma_betapar)

    #convolve two probability distributions for third phase parameters
    #Equations adapted from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004965
    #(Moutoussis et al., 2016)

    # Phase 3 choices of the participant
    if(phase == 3){
      if(conv){

      conv_a <- convolve_distributions(alpha, alpha_v, alpha_v_ref, sigma_alphapar, mu_alphapar)
      conv_b <- convolve_distributions(beta, beta_v, beta_v_ref, sigma_betapar, mu_betapar)

      #create and normalise the third phase distribution
      pabg_shift=(dnorm(alpha_grid,conv_a$convolve_self_m,conv_a$convolve_self_v)+eps)*
        (dnorm(beta_grid,conv_b$convolve_self_m,conv_b$convolve_self_v)+eps)
      pabg_shift=pabg_shift / sum(as.vector(pabg_shift))

      }

    for (t in ((T1+T2)+1):(T1+T2+T3)){

    s1 = as.numeric(data[t, 3]/10);
    o1 = as.numeric(data[t, 4]/10);
    s2 = as.numeric(data[t, 5]/10);
    o2 = as.numeric(data[t, 6]/10);

    rew[t,1] = s1
    rew[t,2] = o1
    rew[t,3] = s2
    rew[t,4] = o2

    if(conv){
    val1 = (alpha_grid*s1) + (beta_grid*max(s1-o1,0)) ;
    val2 = (alpha_grid*s2) + (beta_grid*max(s2-o2,0)) ;
    } else {
    val1 = alpha*s1 + beta*max(s1-o1,0) ;
    val2 = alpha*s2 + beta*max(s2-o2,0) ;
    }

    pchoose3 = (1/(1+exp(-(val1 - val2))));
    tmp3 = ifelse(conv, pchoose3*pabg_shift, pchoose3);
    subject_netp3 = sum(as.vector(tmp3));

    simA[t] = sample(c(1,2),1,prob = c(subject_netp3, 1-subject_netp3));

    if (sim){
      subject_estimate = simA[t];
    } else {
      subject_estimate = data[t, 7]; # subject choice
    }

    if (subject_estimate == 1){
        lik3 = lik3 + log(subject_netp3); # log likelihood
        prob3[t-(T1+T2)] = subject_netp3;
    } else {
        lik3 = lik3 + log(1-subject_netp3);
        prob3[t-(T1+T2)] = 1-subject_netp3;
    }

  }

      if(conv){
        simulatedDF$marginals$Phase3_alpha = colSums(pabg_shift);
        simulatedDF$marginals$Phase3_beta  = rowSums(pabg_shift);

      #find the indices of maximum likelihood after seing the partner
      mu_alphashift   = sum(alpha_grid[1,] * simulatedDF$marginals$Phase3_alpha)
      mu_betashift    = sum(beta_grid[,1] * simulatedDF$marginals$Phase3_beta)

      sddevs2  = std_dev_prob(alpha_grid, beta_grid,
                             mu_alphashift, mu_betashift,
                             simulatedDF$marginals$Phase3_alpha,
                             simulatedDF$marginals$Phase3_beta)

      sigma_alphashift= sddevs2[1]
      sigma_betashift = sddevs2[2]
      simulatedDF$Posteriors[1:4] <- c(mu_alphashift, mu_betashift, sigma_alphashift, sigma_betashift)
     }
    }

  }

  simulatedDF$Participant$PPTActions  <- simA
  simulatedDF$Participant$Correct     <- correct
  simulatedDF$Participant$FixActions  <- simAFix
  simulatedDF$Participant$ProbAction  <- c(prob1, prob2, prob3)

  if(phase == 3 & conv){
  simulatedDF$Shift      <- c(mu_alphashift - alpha,
                              mu_betashift - beta,
                              sigma_alphashift - alpha_v,
                              sigma_betashift - beta_v)
  }

  simulatedDF$LL  = lik1 + lik2 + lik3
  LL = lik1 + lik2 + lik3

  if(plot & conv){
    alphaPlot <- simulatedDF$marginals %>%
      as.data.frame() %>%
      dplyr::select(1, 3, 7, 9) %>%
      pivot_longer(2:4, names_to = 'Phase', values_to = 'Probability') %>%
      mutate(Phase = ifelse(Phase == 'Phase1_alpha', '1', ifelse(Phase == 'Phase2_alpha', '2', '3'))) %>%
      ggplot(aes(alpha_grid, Probability, colour = Phase))+
      geom_line()+
      labs(x = 'Alpha')+
      scale_color_manual(values = c('#095256', '#D33F49', '#087F8C'), labels = c(expression(paste(alpha[ppt])),
                                                                                 expression(paste(alpha[par])),
                                                                                 expression(paste(hat(alpha)[ppt]))))+
      theme_tq()+
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            legend.position.inside = c(0.9, 0.7),
            legend.background = element_rect(colour = 'black'),
            legend.title = element_blank(),
            legend.text = element_text(size = 14))
    betaPlot <- simulatedDF$marginals %>%
      as.data.frame() %>%
      dplyr::select(2, 4, 8, 10) %>%
      pivot_longer(2:4, names_to = 'Phase', values_to = 'Probability') %>%
      ggplot(aes(beta_grid, Probability, colour = Phase))+
      geom_line()+
      labs(x = 'Beta')+
      scale_color_manual(values = c('#095256', '#D33F49', '#087F8C'), labels = c(expression(paste(beta[ppt])),
                                                                                 expression(paste(beta[par])),
                                                                                 expression(paste(hat(beta)[ppt]))))+
      theme_tq()+
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            legend.position.inside = c(0.9, 0.7),
            legend.title = element_blank(),
            legend.background = element_rect(colour = 'black'),
            legend.text = element_text(size = 14))

    if(phase != 1){
    logistic <- simulatedDF$Participant %>%
      as.data.frame() %>%
      filter(Phase == 2) %>%
      dplyr::select(Trial, Correct) %>%
      na.omit() %>%
      mutate(Trial = 1:T2)
    save <- glm(Correct ~ Trial, data = logistic, family = binomial)
    newrange <- seq(1, T2, 1)
    guessdf <- data.frame(Trial = newrange)
    rownames(guessdf) <- guessdf$Trial
    log_curve    <- predict(save, newdata = guessdf, type = 'response') %>% as.data.frame() %>% mutate(Trial = newrange) %>% rename(Prob = 1)
    learningPlot <- ggplot( logistic,
                      aes(Trial, Correct)) +
      geom_point()+
      geom_line(data = log_curve, aes(x = Trial, y = Prob), linewidth = 1)+
      #geom_text(aes(x = 45, y = 0.3, label = paste('Correct = ', round((sum(simulatedDF$Participant$Correct[(T1+1):T2])/(T2-T1))*100,1), '%')), check_overlap = T,  show.legend = F, color ='black')+
      ggpubr::stat_cor(label.y.npc = 0.5, label.x.npc = 0.6)+
      scale_y_continuous(breaks = c(0, 1), labels = c(0,1))+
      labs(x = 'Trial', y = expression(paste('p(Correct | ', theta[ppt], ';', D[par], ')')))+
      scale_color_brewer(palette = 'Set2')+
      theme_tq()+
      theme(legend.position = 'none',
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14))
    } else {
      learningPlot <- ggplot() + theme_void()
    }

    showResults <- ((alphaPlot/betaPlot) | learningPlot) + plot_annotation(tag_levels = 'A')
    simulatedDF$Plot <- showResults
  }

  if(sim){
  return(simulatedDF)
  } else {
  return(LL)
  }

}

ABA_shift_Gen_integratedprior <- function(parms, df, sim = 0, phase = 3, conv=1){

  #Arguments

  #Necessary
  # parms = parameter vector (either 4 or 6 parameters; see conv)
  # df    = dataframe of structure as listed above

  #Optional
  # sim   = do you want the function to simlate data (1=yes, 0=no)
  # plot  = do you want a plot at the end (1=yes, 0=no)
  # phase = how many phases of the game (can only be 1 2 or 3)
  # conv  = is the third phase a convolution of the first and second phase? If so
  #         the function requires 2 extra parameters to define the variance of the
  #         participants beliefs in the first phase to convolve with the partners

  # Initialise
  res = 30; # bounds of belief grid
  inc = 0.125; # resolution of grid
  eps = 0.02/(length(seq(-res, res, inc))^2) #noise floor for

  if (phase == 1) {
    data <- df[df$Phase == 1, ]
  } else if (phase == 2) {
    data <- df[df$Phase %in% c(1, 2), ]
  } else if (phase == 3) {
    data <- df[df$Phase %in% c(1, 2, 3), ]
  }

  T1 <- sum(data$Phase == 1)
  T2 <- if (phase != 1) sum(data$Phase == 2) else 0
  T3 <- if (phase == 3) sum(data$Phase == 3) else 0

  if(phase > 3){stop('This function currently only supports (up to) three-phase estimation')}

  #Phase 1
  alpha <- ifelse(sim, parms[1], (1 / (1 + exp(-parms[1]))) * res)
  beta  <- parms[2]

  if (!sim) {
    parmsP <- c(alpha, beta)
    alpha_v <- exp(parms[3])
    beta_v  <- exp(parms[4])
    parmsP  <- c(parmsP, alpha_v, beta_v)
  }
  if (sim) {
    parmsP <- c(alpha, beta)
    alpha_v <- parms[3]
    beta_v  <- parms[4]
    parmsP  <- c(parmsP, alpha_v, beta_v)
  }

  if(phase != 1){
  #Phase 2
    if(!sim){
      alpha_m_par      = 1/(1+exp(-parms[5]))*30
      beta_m_par       = parms[6]
      alpha_v_ref      = exp(parms[7])
      beta_v_ref       = exp(parms[8])
      w                = 1/(1+exp(-parms[9]))
      parmsP = c(parmsP, alpha_m_par, beta_m_par, alpha_v_ref, beta_v_ref, w)
    }
    if(sim){
      alpha_m_par      = parms[5]
      beta_m_par       = parms[6]
      alpha_v_ref      = parms[7]
      beta_v_ref       = parms[8]
      w                = parms[9]
      parmsP = c(parmsP, alpha_m_par, beta_m_par, alpha_v_ref, beta_v_ref, w)
    }
  }

  #output list
  simulatedDF <- create_output_file(conv=conv, res=res, data=data, inc=inc, T1, T2, T3, parmsP)

  simulatedDF$Priors[1:length(parms)] <- parmsP[1:length(parms)]

  # grid for a subjects beliefs over their preferences in phase 1

  #parameters space of alpha and beta
  grid <- meshgrid(seq(res-res, res, inc),seq(-res, res, inc*2));
  alpha_grid <- grid$x
  beta_grid  <- grid$y

  if(conv){
  #generate standardised grid to form priors for preferences
  pabg <- (dnorm(alpha_grid,alpha,alpha_v)+eps)*(dnorm(beta_grid,beta,beta_v)+eps);
  pabg <- pabg / sum(as.vector(pabg));
  simulatedDF$marginals$Phase1_alpha = colSums(pabg);
  simulatedDF$marginals$Phase1_beta  = rowSums(pabg);
  }

  # initialised dummy values

  lik1        <- 0;   # likelihood for choices in phase 1
  prob1       <- rep(NA, T1)
  lik2        <- 0;   # likelihood for guesses in phase 2
  prob2       <- rep(NA, T2)
  lik3        <- 0;
  prob3       <- rep(NA, T3)
  LL          <- 0

  simA        <- rep(NA, (T1+T2+T3))
  simAFix     <- simA
  correct     <- simA
  rew         <- matrix(NA, nrow = (T1+T2+T3), ncol = 4)

  # Phase 1 choices of the participant

    for (t in 1:T1){

    s1 = as.numeric(data[t, 3]/10);
    o1 = as.numeric(data[t, 4]/10);
    s2 = as.numeric(data[t, 5]/10);
    o2 = as.numeric(data[t, 6]/10);

    rew[t,1] = s1
    rew[t,2] = o1
    rew[t,3] = s2
    rew[t,4] = o2

    if(conv){
    val1 = alpha_grid*s1 + beta_grid*max(s1-o1,0) ;
    val2 = alpha_grid*s2 + beta_grid*max(s2-o2,0) ;
    } else {
    val1 = alpha*s1 + beta*max(s1-o1,0) ;
    val2 = alpha*s2 + beta*max(s2-o2,0) ;
    }

    pchoose1=(1/(1+exp(-(val1 - val2)))); # probability of 1

    tmp_ppt = ifelse(conv, (pchoose1*pabg), pchoose1);

    subject_netp1 = sum(as.vector(tmp_ppt));
    simA[t] = sample(c(1,2),1, prob = c(subject_netp1, 1-subject_netp1));

      if (sim){
        actual_choice = simA[t];
      } else {
        actual_choice = data[t, 7];
      }

      if (actual_choice==1){
      lik1 = lik1 + log(subject_netp1); # log likelihood of 1
      prob1[t] = subject_netp1;
      } else {
      lik1 = lik1 + log(1-subject_netp1);
      prob1[t] = 1-subject_netp1;
      }

    }

####

  # how the experimenter learns how the subject learns online about the partner in phase 2

  #generate standardised grid to form priors for beliefs

  alpha_m_ref <- (alpha * w) + (alpha_m_par * (1-w))
  beta_m_ref  <- (beta * w) + (beta_m_par * (1-w))
  pabg_par    <- (dnorm(alpha_grid,alpha_m_ref,alpha_v_ref)+eps)*(dnorm(beta_grid,beta_m_ref,beta_v_ref)+eps);
  pabg_par    <- pabg_par/sum(as.vector(pabg_par));

  simulatedDF$marginals$Phase2a_alpha = colSums(pabg_par);
  simulatedDF$marginals$Phase2a_beta  = rowSums(pabg_par);

    # Phase2

    if(phase != 1){

    for (t in (T1+1):(T2+T1)){

    s1 = as.numeric(data[t, 3]/10);
    o1 = as.numeric(data[t, 4]/10);
    s2 = as.numeric(data[t, 5]/10);
    o2 = as.numeric(data[t, 6]/10);

    rew[t,1] = s1
    rew[t,2] = o1
    rew[t,3] = s2
    rew[t,4] = o2

    val1 = alpha_grid*s1 + beta_grid*max(s1-o1,0) ;
    val2 = alpha_grid*s2 + beta_grid*max(s2-o2,0) ;

    subject_estimate_pchoose1 = (1/(1+exp(-(val1 - val2))));

    tmp_par = subject_estimate_pchoose1 * pabg_par;
    subject_netp2 = sum(as.vector(tmp_par));

    if(conv){
      tmp_fix = subject_estimate_pchoose1 * pabg;
      subject_netp2_fix = sum(as.vector(tmp_fix));
    } else {
      val1 = alpha*s1 + beta*max(s1-o1,0) ;
      val2 = alpha*s2 + beta*max(s2-o2,0) ;
      subject_netp2_fix = (1/(1+exp(-(val1 - val2))))
    }

    simA[t]    = sample(c(1,2),1,prob = c(subject_netp2, 1-subject_netp2));
    simAFix[t] = sample(c(1,2),1,prob = c(subject_netp2_fix, 1-subject_netp2_fix));

      if (sim){
        subject_estimate = simA[t];
      } else {
        subject_estimate = data[t, 7]; # subject choice
      }

      if (subject_estimate == 1){
        lik2        = lik2 + log(subject_netp2); # log likelihood
        prob2[t-T1] = subject_netp2;
      } else {
        lik2        = lik2 + log(1-subject_netp2);
        prob2[t-T1] = 1-subject_netp2;
      }

    actual_choice = data[t, 8]; # what did our partner 'choose'
    if(actual_choice == subject_estimate){correct[t] <- 1} else {correct[t] <- 0}

      if (actual_choice==1){
        pchoose2=(1./(1+exp(-(val1 - val2)))); # probability of 1
      } else {
        pchoose2=(1./(1+exp(-(val2 - val1)))); # probability of 2
      }

    pabg_par = pchoose2*pabg_par; # Bayes rule
    pabg_par = pabg_par / sum(as.vector(pabg_par)); #normalised distribution

    }

    simulatedDF$marginals$Phase2_alpha = colSums(pabg_par);
    simulatedDF$marginals$Phase2_beta  = rowSums(pabg_par);

    #find the indices of maximum likelihood after seeing the partner
    mu_alphapar   = sum(alpha_grid[1, ] * simulatedDF$marginals$Phase2_alpha)
    mu_betapar    = sum(beta_grid[, 1]  * simulatedDF$marginals$Phase2_beta)
    sddevs  = std_dev_prob(alpha_grid, beta_grid,
                           mu_alphapar, mu_betapar,
                           simulatedDF$marginals$Phase2_alpha,
                           simulatedDF$marginals$Phase2_beta)
    sigma_alphapar= sddevs[1]
    sigma_betapar = sddevs[2]

    simulatedDF$PartnerParms[1:4] <- c(mu_alphapar, mu_betapar, sigma_alphapar, sigma_betapar)

    #convolve two probability distributions for third phase parameters
    #Equations adapted from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004965
    #(Moutoussis et al., 2016)

    # Phase 3 choices of the participant
    if(phase == 3){
      if(conv){

      conv_a <- convolve_distributions(alpha, alpha_v, alpha_v_ref, sigma_alphapar, mu_alphapar)
      conv_b <- convolve_distributions(beta, beta_v, beta_v_ref, sigma_betapar, mu_betapar)

      #create and normalise the third phase distribution
      pabg_shift=(dnorm(alpha_grid,conv_a$convolve_self_m,conv_a$convolve_self_v)+eps)*
        (dnorm(beta_grid,conv_b$convolve_self_m,conv_b$convolve_self_v)+eps)
      pabg_shift=pabg_shift / sum(as.vector(pabg_shift))

      }

    for (t in ((T1+T2)+1):(T1+T2+T3)){

    s1 = as.numeric(data[t, 3]/10);
    o1 = as.numeric(data[t, 4]/10);
    s2 = as.numeric(data[t, 5]/10);
    o2 = as.numeric(data[t, 6]/10);

    rew[t,1] = s1
    rew[t,2] = o1
    rew[t,3] = s2
    rew[t,4] = o2

    if(conv){
    val1 = (alpha_grid*s1) + (beta_grid*max(s1-o1,0)) ;
    val2 = (alpha_grid*s2) + (beta_grid*max(s2-o2,0)) ;
    } else {
    val1 = alpha*s1 + beta*max(s1-o1,0) ;
    val2 = alpha*s2 + beta*max(s2-o2,0) ;
    }

    pchoose3 = (1/(1+exp(-(val1 - val2))));
    tmp3 = ifelse(conv, pchoose3*pabg_shift, pchoose3);
    subject_netp3 = sum(as.vector(tmp3));

    simA[t] = sample(c(1,2),1,prob = c(subject_netp3, 1-subject_netp3));

    if (sim){
      subject_estimate = simA[t];
    } else {
      subject_estimate = data[t, 7]; # subject choice
    }

    if (subject_estimate == 1){
        lik3 = lik3 + log(subject_netp3); # log likelihood
        prob3[t-(T1+T2)] = subject_netp3;
    } else {
        lik3 = lik3 + log(1-subject_netp3);
        prob3[t-(T1+T2)] = 1-subject_netp3;
    }

  }

      if(conv){
        simulatedDF$marginals$Phase3_alpha = colSums(pabg_shift);
        simulatedDF$marginals$Phase3_beta  = rowSums(pabg_shift);

      #find the indices of maximum likelihood after seing the partner
      mu_alphashift   = sum(alpha_grid[1,] * simulatedDF$marginals$Phase3_alpha)
      mu_betashift    = sum(beta_grid[,1] * simulatedDF$marginals$Phase3_beta)

      sddevs2  = std_dev_prob(alpha_grid, beta_grid,
                             mu_alphashift, mu_betashift,
                             simulatedDF$marginals$Phase3_alpha,
                             simulatedDF$marginals$Phase3_beta)

      sigma_alphashift= sddevs2[1]
      sigma_betashift = sddevs2[2]
      simulatedDF$Posteriors[1:4] <- c(mu_alphashift, mu_betashift, sigma_alphashift, sigma_betashift)
     }
    }

  }

  simulatedDF$Participant$PPTActions  <- simA
  simulatedDF$Participant$Correct     <- correct
  simulatedDF$Participant$FixActions  <- simAFix
  simulatedDF$Participant$ProbAction  <- c(prob1, prob2, prob3)

  if(phase == 3 & conv){
  simulatedDF$Shift      <- c(mu_alphashift - alpha,
                              mu_betashift - beta,
                              sigma_alphashift - alpha_v,
                              sigma_betashift - beta_v)
  }

  simulatedDF$LL  = lik1 + lik2 + lik3
  LL = lik1 + lik2 + lik3

  if(sim){
  return(simulatedDF)
  } else {
  return(LL)
  }

}

clean_gorilla_csv <- function(x){

  y <- x %>%
  rename(ID = `Participant.Public.ID`,
         S1 = `playerPoints_option1`,
         O1 = `partnerPoints_option1`,
         S2 = `playerPoints_option2`,
         O2 = `partnerPoints_option2`,
         choice = `selectedOption_player`,
         real_answer = `realAnswer`,
         correct = `correctGuess`,
         attention = `attention_check_correct`,
         SI = `inferenceResponse_Selfish`,
         HI = `inferenceResponse_Harm`,
         agency = `agencyResponse`,
         RT = `trialDuration`
         ) %>%
  dplyr::select(ID, trial, S1, O1, S2, O2,
                choice, real_answer, correct, SI, HI, agency, classification, display, attention,
                server_alpha_ppt, server_beta_ppt, server_alpha_par, server_beta_par, RT) %>%
  filter(display %in% c('playerChoice', 'playerGuess', 'playerChoice2',
                        'inference', 'classification', 'agency', 'matching') |
         !is.na(attention)) %>%
  group_by(ID) %>%
  mutate(real_answer = ifelse(real_answer == 'Option 2', 2, 1),
         Phase = ifelse(display == 'playerChoice', 1,
                        ifelse(display == 'playerChoice2', 3,
                               ifelse(display == 'playerGuess', 2, NA))),
         classification = ifelse(classification == "To stop me from earning money", 'Competitive',
                          ifelse(classification == "To earn as much money for themselves as possible", 'Individualist',
                          ifelse(classification == "To share the money between us evenly", 'Prosocial', 'NA')))) %>%
  slice(-1) %>%
  dplyr::select(-display) %>%
    group_by(ID, Phase) %>%
  mutate(correctSum = ifelse(is.na(correct), 0, correct),
         correctSum = sum(correctSum)) %>%
  group_by(ID) %>%
  mutate(SI = ifelse(is.na(SI), 0, SI),
         SI = sum(SI),
         HI = ifelse(is.na(HI), 0, HI),
         HI = sum(HI),
         agency = ifelse(is.na(agency), 0, agency),
         agency = sum(agency),
         attention = ifelse(is.na(attention), 0, attention),
         attention = sum(attention),
         server_alpha_ppt= ifelse(is.na(server_alpha_ppt), 0, server_alpha_ppt),
         server_alpha_ppt = sum(server_alpha_ppt),
         server_beta_ppt= ifelse(is.na(server_beta_ppt), 0, server_beta_ppt),
         server_beta_ppt = sum(server_beta_ppt),
         server_alpha_par= ifelse(is.na(server_alpha_par), 0, server_alpha_par),
         server_alpha_par = sum(server_alpha_par),
         server_beta_par= ifelse(is.na(server_beta_par), 0, server_beta_par),
         server_beta_par = sum(server_beta_par))
 return(y)

}

# Function to filter rows past the first instance of 'end' in the 'display' column
filter_past_end <- function(df, column_name, entry) {
  # Find the index of the first occurrence of 'end'
  first_end_index <- which(df[[column_name]] == entry)[1]

  # Subset the data frame to include rows up to the first 'end'
  if (!is.na(first_end_index)) {
    df_filtered <- df[1:first_end_index, ]
  } else {
    df_filtered <- df
  }

  return(df_filtered)
}

simulate_simplex <- function(data=NULL, trials=100, res=0.25, v=1){

  int=res
  alpha = seq(0, 30, int)
  beta  = seq(-30, 30, int*2)

  if(is.null(data)){

    t=trials
    env = data.frame(
      S1 = sample(0:12, t, T),
      O1 = sample(0:12, t, T),
      S2 = sample(0:12, t, T),
      O2 = sample(0:12, t, T)
    )

    for(i in 1:nrow(env)){
      s1 = env[i,'S1'];
      s2 = env[i,'S2'];
      o1 = env[i,'O1'];
      o2 = env[i,'O2'];
      env[i,'S1'] <- ifelse(s1<o1, o1, s1); #print(s1<o1)
      env[i,'O1'] <- ifelse(s1<o1, s1, o1);
      env[i,'S2'] <- ifelse(s2<o2, o2, s2);
      env[i,'O2'] <- ifelse(s2<o2, s2, o2);
    }

  } else {

    first_id = unique(data$ID)[1]
    env   = data %>%
      filter(ID==first_id, Phase ==2) %>%
      na.omit()

  }

  grid       = meshgrid(alpha,beta);
  alpha_grid = grid$x
  beta_grid  = grid$y
  simpV  = matrix(NA, nrow = length(alpha), ncol = length(beta));
  rownames(simpV)=alpha;
  colnames(simpV)=beta
  simpS  = simpV
  simpO  = simpV

  vS = v # sd of prior

  for(i in 1:length(alpha)){
  a = alpha[i];
    for(j in 1:length(beta)){
      b = beta[j];
      prior = (dnorm(alpha_grid,a,vS)*dnorm(beta_grid,b,vS));
      prior = prior / sum(as.vector(prior))
      s=0
      o=0
      for(k in 1:nrow(env)){
        s1 = as.numeric(env[k,'S1'])
        s2 = as.numeric(env[k,'S2'])
        o1 = as.numeric(env[k,'O1'])
        o2 = as.numeric(env[k,'O2'])
        u1 = (a * s1) + (b * max(s1-o1, 0))
        u2 = (a * s2) + (b * max(s2-o2, 0))
        p1 = 1/(1+exp(-(u1-u2))) * sum(as.vector(prior))
        p2 = 1-p1; p2 = ifelse(p2 < 0, 0, p2)
        s  = s+sample(c(s1, s2), 1, prob = c(p1, p2))
        o  = o+sample(c(o1, o2), 1, prob = c(p1, p2))
      }
    simpV[i,j]   <- s-o
    simpS[i,j]   <- s
    simpO[i,j]   <- o
    }
  print(paste(i,' of ', length(alpha)))
  }

  simpLongV <- simpV %>%
    as.data.frame() %>%
    mutate(alpha = alpha) %>%
    pivot_longer(1:length(beta), names_to = 'beta', values_to = 'val') %>%
    mutate(beta = as.numeric(beta))
  simpLongS <- simpS %>%
    as.data.frame() %>%
    mutate(alpha = alpha) %>%
    pivot_longer(1:length(beta), names_to = 'beta', values_to = 'val') %>%
    mutate(beta = as.numeric(beta))
  simpLongO <- simpO %>%
    as.data.frame() %>%
    mutate(alpha = alpha) %>%
    pivot_longer(1:length(beta), names_to = 'beta', values_to = 'val') %>%
    mutate(beta = as.numeric(beta))
  return(list(diff=simpLongV, self=simpLongS, other=simpLongO))
}

extract_and_combine <- function(results1, column_name, n_length) {
  # Determine the number of vectors
  n_results1 <- length(results1)

  # Function to extract a specific column from each result and ensure no zeros
  extract_column <- function(results, i, column_name, n_length) {
    # Extract the vector
    column_vector <- results[i,][[1]][[1]][,,1][[column_name]]

    # Remove zeros
    column_vector_non_zero <- column_vector[column_vector != 0]

    # Ensure there are enough non-zero elements to fill the desired length
    if (length(column_vector_non_zero) < n_length) {
      warning(paste("Not enough non-zero elements to fill length", n_length, "for row", i))
      # Fill with non-zero elements followed by NAs
      result_vector <- c(column_vector_non_zero, rep(NA, n_length - length(column_vector_non_zero)))
    } else {
      # Take the first n_length non-zero elements
      result_vector <- column_vector_non_zero[1:n_length]
    }

    return(result_vector)
  }

  # Extract the specified column vectors for the result set
  results1_matrix <- t(sapply(1:n_results1, function(i) extract_column(results1, i, column_name, n_length)))

  return(results1_matrix)
}

calculate_bic <- function(LL, k, n) {
  BIC <- -LL + k * log(n)
  return(BIC)
}

extract_parameters <- function(x, ID) {
  # Check the best fitting model
  y <- which.max(x$cbm[,,1]$output[,,1]$exceedance.prob)

  # Determine the number of rows needed for x_parms
  num_rows <- nrow(x$cbm[,,1]$output[,,1]$responsibility)
  num_cols <- ncol(x$cbm[,,1]$output[,,1]$parameters[y,][[1]][[1]])

  # Extract the relevant parameters
  x_parms <- as.data.frame(matrix(NA, nrow = num_rows, ncol = num_cols+1))
  x_parms[,1:num_cols] <- x$cbm[,,1]$output[,,1]$parameters[y,][[1]][[1]]
  x_parms[,num_cols+1] <- deparse(substitute(x))

  # Create joint_parms data frame with the required transformations
  if (num_cols > 6) {
    joint_parms <- data.frame(
      ID        = ID,
      alpha     = (1 / (1 + exp(-x_parms[,1]))) * 30,
      beta      = x_parms[,2],
      alpha_v   = exp(x_parms[,3]),
      beta_v    = exp(x_parms[,4]),
      alpha_par = (1 / (1 + exp(-x_parms[,5]))) * 30,
      beta_par  = x_parms[,6],
      alpha_ref = exp(x_parms[,7]),
      beta_ref  = exp(x_parms[,8]),
      group     = x_parms[,num_cols + 1],
      Model     = y
    ) %>%
      distinct()
  } else {
    joint_parms <- data.frame(
      ID        = ID,
      alpha     = (1 / (1 + exp(-x_parms[,1]))) * 30,
      beta      = x_parms[,2],
      alpha_v   = exp(x_parms[,3]),
      beta_v    = exp(x_parms[,4]),
      alpha_par = NA,
      beta_par  = NA,
      alpha_ref = exp(x_parms[,5]),
      beta_ref  = exp(x_parms[,6]),
      group     = x_parms[,num_cols + 1],
      Model     = y
    ) %>%
      distinct()
  }

  return(joint_parms)
}

# Function to calculate KL divergence
calculate_KL_divergence <- function(prob_A, prob_B) {
  # Ensure that both distributions have the same support
  if (length(prob_A) != length(prob_B)) {
    stop("Probability distributions must have the same length.")
  }

  # Remove zeros from probabilities to avoid NaN values in the logarithm
  prob_A <- prob_A[prob_A > 0]
  prob_B <- prob_B[prob_B > 0]

  # Calculate KL divergence
  KL_divergence <- sum(prob_A * log(prob_A / prob_B))
  return(KL_divergence)
}

meshgrid <- function(x,y){
  mesh <- list()
  m=length(x);
  n=length(y);
  mesh$x <- matrix(rep(x,each=n),nrow=n);
  mesh$y <- matrix(rep(y,m),nrow=n)
  return(mesh)
}

create_output_file <- function(conv=T, res, data, inc, T1, T2, T3, parmsP){

  if(conv){
  simulatedDF <- list(

    marginals = data.frame(
    alpha_grid   = seq(0, res, inc),
    beta_grid    = seq(-res, res, inc*2),
    Phase1_alpha = rep(NA, length(seq(0, res, inc))),
    Phase1_beta  = rep(NA, length(seq(0, res, inc))),
    Phase2a_alpha= rep(NA, length(seq(0, res, inc))),
    Phase2a_beta = rep(NA, length(seq(0, res, inc))),
    Phase2_alpha = rep(NA, length(seq(0, res, inc))),
    Phase2_beta  = rep(NA, length(seq(0, res, inc))),
    Phase3_alpha = rep(NA, length(seq(0, res, inc))),
    Phase3_beta  = rep(NA, length(seq(0, res, inc)))),

    Participant  = data.frame(
        ID           = data[,1],
        Trial        = 1:(T1+T2+T3),
        s1           = data[,3],
        o1           = data[,4],
        s2           = data[,5],
        o2           = data[,6],
        PPTActions   = rep(NA, (T1+T2+T3)),
        Partner      = data[,8],
        Correct      = rep(NA, (T1+T2+T3)),
        FixActions   = rep(NA, (T1+T2+T3)),
        ProbAction   = rep(NA, (T1+T2+T3)),
        Phase        = c(rep(1, T1), rep(2, T2), rep(3, T3))
    ),

    Priors       = rep(NA, length(parmsP)),
    PartnerParms = rep(NA, 4),
    Posteriors   = rep(NA, 4),
    Shift        = rep(NA, 4),

    Plot         = NA,
    LL           = NA
  )
  } else {
  simulatedDF <- list(

    marginals = data.frame(
    alpha_grid   = seq(0, res, inc),
    beta_grid    = seq(-res, res, inc*2),
    Phase1_alpha = parmsP[1],
    Phase1_beta  = parmsP[2],
    Phase2a_alpha= rep(NA, length(seq(0, res, inc))),
    Phase2a_beta = rep(NA, length(seq(0, res, inc))),
    Phase2_alpha = rep(NA, length(seq(0, res, inc))),
    Phase2_beta  = rep(NA, length(seq(0, res, inc)))),

    Participant  = data.frame(
        ID           = data[,1],
        Trial        = 1:(T1+T2+T3),
        s1           = data[,3],
        o1           = data[,4],
        s2           = data[,5],
        o2           = data[,6],
        PPTActions   = rep(NA, (T1+T2+T3)),
        Partner      = data[,8],
        Correct      = rep(NA, (T1+T2+T3)),
        FixActions   = rep(NA, (T1+T2+T3)),
        ProbAction   = rep(NA, (T1+T2+T3)),
        Phase        = c(rep(1, T1), rep(2, T2), rep(3, T3))
    ),

    Priors       = rep(NA, length(parmsP)),
    PartnerParms = rep(NA, 4),
    Posteriors   = rep(NA, 4),

    Plot         = NA,
    LL           = NA
  )
  }

  return(simulatedDF)

}

std_dev_prob <- function(alpha_grid, beta_grid, mu_alpha, mu_beta, prob_alpha, prob_beta){

  sdalpha = rep(NA, length(alpha_grid[,1]))
  sdbeta  = sdalpha

  for (i in 1:length(sdalpha)){
    sdalpha[i] = ((alpha_grid[1, i] - mu_alpha)^2) * prob_alpha[i]
    sdbeta[i]  = ((beta_grid[i, 1]  - mu_beta)^2) * prob_beta[i]
  }

  sdalpha  <- sum(sdalpha)
  sdalpha  <- sqrt(sdalpha)
  sdbeta   <- sum(sdbeta)
  sdbeta   <- sqrt(sdbeta)

  return(c(sdalpha, sdbeta))

}

convolve_distributions <- function(self, self_v, self_v_ref, sigma_par, mu_par) {

  # Precompute frequently used terms
  self_v_dash <- self_v^2
  self_v_ref_dash <- self_v_ref^2
  sigma_par_dash <- sigma_par^2

  # Compute convolve_self_v
  self_var <- self_v_dash^-1 + ((2 * self_v_ref_dash + sigma_par_dash)^-1)
  convolve_self_v_v <- self_var^-1
  convolve_self_v <- sqrt(self_var^-1)

  # Compute convolve_self_m
  convolve_self_m <- convolve_self_v_v * ((self_v_dash^-1) * self + ((2 * self_v_ref_dash + sigma_par_dash)^-1) * mu_par)

  return(list(convolve_self_v = convolve_self_v, convolve_self_m = convolve_self_m))
}

simulate_phase_decisions <- function(parms, data, phase = 2){

  # Initialise

  res = 30; # resolution of belief grid
  T1  = length(data %>% filter(Phase == phase) %>% rownames());  # trials for phase 1

  #Phase 1
  alpha            = as.numeric(parms[1])
  beta             = as.numeric(parms[2])

  # initialised dummy values
  decisions        <- data.frame(
    ppt1 = rep(NA, T1),
    par1 = rep(NA, T1),
    ppt2 = rep(NA, T1),
    par2 = rep(NA, T1),
    Ac   = rep(NA, T1)
  )

  if(length(parms) == 3){decisions$type = parms[3]}

  # Phase 1 choices of the participant

    for (t in 1:T1){

    s1 = as.numeric(data[t, 3]/10);
    o1 = as.numeric(data[t, 4]/10);
    s2 = as.numeric(data[t, 5]/10);
    o2 = as.numeric(data[t, 6]/10);

    decisions[t, 1] = s1*10
    decisions[t, 2] = o1*10
    decisions[t, 3] = s2*10
    decisions[t, 4] = o2*10

    val1 = alpha*s1 + beta*max(s1-o1,0) ;
    val2 = alpha*s2 + beta*max(s2-o2,0) ;

    pchoose1=(1/(1+exp(-(val1 - val2)))); # probability of 1
    simA = sample(c(1,2),1, prob = c(pchoose1, 1-pchoose1));

    decisions[t, 5] = simA;

    }

      return(decisions)


}