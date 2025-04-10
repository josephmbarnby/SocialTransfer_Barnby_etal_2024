#Analysis script for:
#Transfer of social information reduces uncertainty about the self and others
#CC Joe Barnby 2024

# Libraries ---------------------------------------------------------------
rm(list=ls())
setwd('~/SocialTransfer_Barnby_etal_2024/')
library(R.matlab)  # Functions for handling mat files etc.
library(Hmisc)
library(ppcor)
library(tidyverse)
library(foreach)
library(doParallel)
library(tidyquant)
source('UtilityFunctions.R')
library(patchwork)
library(ggpubr)
library(BayesianFirstAid)

# Load Data ---------------------------------------------------------------

intent_dat <- read.csv('Data/intent_dat.csv') %>% dplyr::select(-X)

## Colours for groups in plots ------------------------------------------------------

colour_group <- c('#A81ADB', '#2BD8FF')

# Sample characteristics --------------------------------------------------

demos_dat <- intent_dat %>%
  mutate(group = ifelse(nchar(ID) < 10, 'BPD', 'CON'),
         across(Gender:CTQtot, ~ifelse(.x==9999, NA, .x))) %>%
  dplyr::select(ID, Gender:A_ETS_Credulity, group) %>%
  distinct()

demos_dat %>%
  split(.$group) %>%
  map(summary)

#Test demos
t.test(Gender ~ group, demos_dat)
t.test(Age ~ group, demos_dat)
t.test(YearsinEdu ~ group, demos_dat)
t.test(SocialDeprivationRank ~ group, demos_dat)

sd(na.omit(demos_dat$Age[demos_dat$group=='CON']))
sd(na.omit(demos_dat$YearsinEdu[demos_dat$group=='BPD']))
sd(na.omit(demos_dat$SocialDeprivationRank[demos_dat$group=='CON']))

sd(na.omit(demos_dat$A_ETS_Credulity[demos_dat$group=='CON']))
sd(na.omit(demos_dat$A_ETS_Credulity[demos_dat$group=='BPD']))

## Psychometric differences ------------------------------------------------

t.test(RGPTSB ~ group, demos_dat)

for(i in 14:(length(demos_dat)-1)){
  x <- scale(demos_dat[,i])
  demos_dat[,i] <- as.vector(x)
}

ggplot(demos_dat %>%
         pivot_longer(c(
                        MZQ_TotalScore, CAMSQ_Self, CAMSQ_Other, A_ETS_Trust, A_ETS_Mistrust, A_ETS_Credulity,
                        RGPTSA, RGPTSB, CTQtot),
                      names_to = 'Var', values_to = 'Val'),
       aes(Var, Val, colour = group))+
  stat_summary()+
  scale_color_manual(values=colour_group)+
  stat_compare_means(label='p.signif', label.y = 1, show.legend = F, paired = F, method = 't.test')+
  scale_x_discrete(labels = rev(c('Persecutory\nIdeation', 'Referential\nIdeation', 'Mentalising', 'Childhood\nTrauma', 'CMS\n(Self)', 'CMS\n(Other)',
                                  'ETS Trust', 'ETS Mistrust', 'ETS Cred.')))+
  coord_flip()+
  theme_bw(base_size=18)+
  theme(axis.title.y = element_blank())

# Behavioural analysis ----------------------------------------------------

beh_anal <- intent_dat %>%
  mutate(group=ifelse(nchar(ID)<10, 'BPD', 'CON'),
         distancebeta=abs(server_beta_par-server_beta_ppt),
         distancealpha=abs(server_alpha_par-server_alpha_ppt),
         across(Gender:CTQtot, ~ifelse(.x==9999, NA, .x))
         ) %>%
  filter(Phase==2)%>%
  dplyr::select(ID, SI:last_col(), -RT) %>%
  distinct()

## Psychometric curves and correctness --------

beh_anal %>%
  group_by(group) %>%
  summarise(cor = mean(correctSum)/54,
            sd = sd(correctSum)/54)

psych_curves <- intent_correctSumpsych_curves <- intent_dat %>%
  filter(Phase == 2) %>%
  group_by(ID) %>%
  mutate(group = ifelse(nchar(ID)<10,'BPD','CON'),
         corSum = sum(correct),
         corSum =corSum/54,
         across(Gender:A_ETS_Credulity, ~ifelse(.x==9999, NA, .x)))


summary(lm(correctSum ~ group,
            data = beh_anal))

summary(lme4::glmer(correct ~ trial + (1|ID),
            data = psych_curves,
            family = "binomial"))

summary(lme4::glmer(correct ~ trial + group + (1|ID),
            data = psych_curves,
            family = "binomial"))

BPD_curve <- glm(correct ~ trial,
    data = psych_curves %>% filter(group=='BPD'),
    family = "binomial")
CON_curve <- glm(correct ~ trial,
    data = psych_curves %>% filter(group=='CON'),
    family = "binomial")

pred_time_dat <- data.frame(trial=seq(0, 54, 0.01))

BPD_prob <- BPD_curve %>% predict(pred_time_dat, type = "response")
CON_prob <- CON_curve %>% predict(pred_time_dat, type = "response")

pred_time_dat$BPD_prob = BPD_prob
pred_time_dat$CON_prob = CON_prob
pred_time_dat <- pred_time_dat %>%
  pivot_longer(2:3, names_to = 'Group', values_to = 'p(Correct)') %>%
  mutate(Group=ifelse(Group=='BPD_prob', 'BPD', 'CON'))

log_corp1 <- ggplot(pred_time_dat,
       aes(trial, `p(Correct)`, colour = Group))+
  geom_smooth(se = T, size = 1) +  # Smoother lines without confidence intervals
  scale_color_manual(values=colour_group)+
  labs(x = 'Trial', y = 'P(Correct)') +  # Proper axis labels
  theme_bw() +
  theme(
    text = element_text(size = 16),
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(colour = 'black', fill = 'white'),
    legend.title = element_blank(),  # Improve legend title readability
    legend.text = element_text(size = 12)  # Improve legend text readability
  )+
ggplot(psych_curves %>%
         dplyr::select(group, ID, corSum) %>%
         distinct(),
       aes(group, corSum, fill = group))+
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, height = 0.2, size = 2) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0.5, size = 1.2)+
  scale_fill_manual(values=colour_group)+
  scale_y_continuous(expand=c(0,0))+
  #stat_compare_means(label.y = 0.25, size = 5)+
  coord_cartesian(y = c(0, 1))+
  labs(x = 'Group', y = '% Correct')+
  theme_bw()+
  theme(
      legend.position = 'none',
      text = element_text(size = 18),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5)
    )

log_corp1

## Behavioural choice ------------------------------------------------------

choices_ppts_1 <- intent_dat %>%
  filter(Phase==1) %>%
  mutate(group=ifelse(nchar(ID)<10, 'BPD', 'CON'),
         relation = case_when(
              (S1 == O1 & S2 > O2 & S2 > S1 & choice==1) | (S2 == O2 & S1 > O1 & S1 > S2 & choice == 2) ~ "Pi",
              (S1 == O1 & S2 > O2 & O1 > O2 & choice==1) | (S2 == O2 & S1 > O1 & O2 > O1 & choice == 2) ~ "Pc",
              (S1 >  O1 & S2 > O2 & S1 > S2 & choice==1) | (S1 >  O1 & S2 > O2 & S2 > S1 & choice == 2) ~ 'Ic',
              TRUE ~ NA
            )) %>%
  dplyr::select(relation, group, trial, ID) %>%
  distinct()

summary_choices_ppts_1 <- choices_ppts_1 %>%
  group_by(ID, relation, group) %>%
  summarise(total=n()) %>%
  pivot_wider(id_cols = c(ID, group), names_from = 'relation', values_from = 'total') %>%
  mutate(Pc=ifelse(is.na(Pc), 0, Pc),
         Ic=ifelse(is.na(Ic), 0, Ic),
         Pi=ifelse(is.na(Pi), 0, Pi)
         )

summary_choices_ppts_1 %>%
  group_by(group) %>%
  mutate(
         meanPc = mean(Pc),
         meanIc = mean(Ic),
         meanPi = mean(Pi),
         sdPc = sd(Pc),
         sdIc = sd(Ic),
         sdPi = sd(Pi),
         ) %>%
  dplyr::select(meanPc:sdPi, group) %>%
  distinct() %>%
  as.data.frame()

t.test(Pc ~ group, data=summary_choices_ppts_1)
t.test(Pi ~ group, data=summary_choices_ppts_1)
t.test(Ic ~ group, data=summary_choices_ppts_1)

choices_ppts_3 <- intent_dat %>%
  filter(Phase==3) %>%
  mutate(group=ifelse(nchar(ID)<10, 'BPD', 'CON'),
         relation = case_when(
              (S1 == O1 & S2 > O2 & S2 > S1 & choice==1) | (S2 == O2 & S1 > O1 & S1 > S2 & choice == 2) ~ "Pi",
              (S1 == O1 & S2 > O2 & O1 > O2 & choice==1) | (S2 == O2 & S1 > O1 & O2 > O1 & choice == 2) ~ "Pc",
              (S1 >  O1 & S2 > O2 & S1 > S2 & choice==1) | (S1 >  O1 & S2 > O2 & S2 > S1 & choice == 2) ~ 'Ic',
              TRUE ~ NA
            )) %>%
  dplyr::select(relation, group, trial, ID) %>%
  distinct()

summary_choices_ppts_3 <- choices_ppts_3 %>%
  group_by(ID, relation, group) %>%
  summarise(total=n()) %>%
  pivot_wider(id_cols = c(ID, group), names_from = 'relation', values_from = 'total') %>%
  mutate(
         Pc=ifelse(is.na(Pc), 0, Pc),
         Ic=ifelse(is.na(Ic), 0, Ic),
         Pi=ifelse(is.na(Pi), 0, Pi)
         )

summary_choices_ppts_3 %>%
  group_by(group) %>%
  mutate(
         meanPc = mean(Pc),
         meanIc = mean(Ic),
         meanPi = mean(Pi),
         sdPc = sd(Pc),
         sdIc = sd(Ic),
         sdPi = sd(Pi),
         ) %>%
  dplyr::select(meanPc:sdPi, group) %>%
  distinct() %>%
  as.data.frame()

t.test(Pc ~ group, data=summary_choices_ppts_3)
t.test(Pi ~ group, data=summary_choices_ppts_3)
t.test(Ic ~ group, data=summary_choices_ppts_3)

choices_ppts_diff <- intent_dat %>%
  filter(Phase%in%c(1,3)) %>%
  mutate(group=ifelse(nchar(ID)<10, 'BPD', 'CON'),
         relation = case_when(
              (S1 == O1 & S2 > O2 & S2 > S1 & choice==1) | (S2 == O2 & S1 > O1 & S1 > S2 & choice == 2) ~ "Pi",
              (S1 == O1 & S2 > O2 & O1 > O2 & choice==1) | (S2 == O2 & S1 > O1 & O2 > O1 & choice == 2) ~ "Pc",
              (S1 >  O1 & S2 > O2 & S1 > S2 & choice==1) | (S1 >  O1 & S2 > O2 & S2 > S1 & choice == 2) ~ 'Ic',
              TRUE ~ NA
            )) %>%
  dplyr::select(relation, group, trial, Phase, ID) %>%
  distinct()

summary_choices_ppts_diff <- choices_ppts_diff %>%
  group_by(ID, relation, group, Phase) %>%
  summarise(total=n()) %>%
  pivot_wider(id_cols = c(ID, group, Phase), names_from = 'relation', values_from = 'total') %>%
  mutate(
         Pc=ifelse(is.na(Pc), 0, Pc),
         Ic=ifelse(is.na(Ic), 0, Ic),
         Pi=ifelse(is.na(Pi), 0, Pi)
         )

summary_choices_ppts_diff %>%
  group_by(group, Phase) %>%
  mutate(
         meanPc = mean(Pc),
         meanIc = mean(Ic),
         meanPi = mean(Pi),
         sdPc = sd(Pc),
         sdIc = sd(Ic),
         sdPi = sd(Pi),
         ) %>%
  dplyr::select(meanPc:sdPi, group) %>%
  distinct() %>%
  as.data.frame()

t.test(Pc ~ Phase, data=summary_choices_ppts_diff %>% filter(group=='CON'))
t.test(Pi ~ Phase, data=summary_choices_ppts_diff %>% filter(group=='CON'))
t.test(Ic ~ Phase, data=summary_choices_ppts_diff %>% filter(group=='CON'))

t.test(Pc ~ Phase, data=summary_choices_ppts_diff %>% filter(group=='BPD'))
t.test(Pi ~ Phase, data=summary_choices_ppts_diff %>% filter(group=='BPD'))
t.test(Ic ~ Phase, data=summary_choices_ppts_diff %>% filter(group=='BPD'))

# Model based analysis ----------------------------------------------------

## Simplex ----------

simpLong <- simulate_simplex(intent_dat, res = 0.5, v = 0.1)

#place participants on grid
simpPPTs <- ggplot(simpLong$diff, aes(alpha, beta, fill = val)) +
  geom_tile() +
  scale_fill_gradient(low = 'black', high = 'white', name = expression(paste(Delta, 'Value'))) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = expression(paste(alpha)), y = expression(paste(beta))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_blank())
simpPPTs

## Model Simulations -------------------------------------------------------

s_uncertainty <- seq(1, 15, 1)
o_uncertainty <- seq(1, 15, 1)

loop_values <- expand.grid(s_uncertainty, o_uncertainty)
shift_beta  <- matrix(NA,
                      nrow = nrow(loop_values),
                      ncol = (ncol(loop_values) + ncol(loop_values) + 2)
                      )

for(i in 1:nrow(loop_values)){

  par_alpha <- intent_dat_126[intent_dat_126$ID=='PDC195',]$server_alpha_par[1]
  par_beta  <- intent_dat_126[intent_dat_126$ID=='PDC195',]$server_beta_par[2]

  print(loop_values[i,])
  parms <- c(par_alpha, -20, 1, loop_values[i,1], 1, loop_values[i,2])
  run1  <- ABA_shift_Gen(parms,
                        intent_dat_126[intent_dat_126$ID=='PDC195',],
                        sim=1,
                        plot = 0)
  ppt_act   <- simulate_phase_decisions(parms[1:2], intent_dat_126[intent_dat_126$ID=='PDC195',])
  par_act   <- simulate_phase_decisions(c(par_alpha, par_beta), intent_dat_126[intent_dat_126$ID=='PDC195',])

  match     <- ifelse(ppt_act==par_act, 1, 0); matched = sum(match[,5])

  shift_beta[i,1] <- run1$Shift[2]/5
  shift_beta[i,2] <- run1$Shift[4]
  shift_beta[i,3] <- matched
  shift_beta[i,4] <- loop_values[i,1]
  shift_beta[i,5] <- loop_values[i,2]

}

colnames(shift_beta)      <- c('m_shift', 'sd_shift', 'match', 's_u', 'o_u', 'disp')
shift_beta_l <- shift_beta %>%
  as.data.frame() %>%
  mutate(s_u = factor(s_u),
         o_u = factor(o_u)) %>%
  pivot_longer(1:2, names_to = 'parameter', values_to = 'shift')

ggplot(shift_beta_l %>%
         filter(parameter == 'sd_shift'),
       aes(s_u, o_u, fill = shift))+
  scale_fill_gradient(low = 'black',high = 'white',
                      name = expression(paste(Delta, theta[ppt]^sigma)))+
  geom_tile()+
  labs(y = expression(paste('Other Unc.')),
       x = expression(paste('Self Unc.'))
       ) +
  scale_y_discrete(breaks = c(1, 5, 10, 15)) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  theme_bw(base_size=18)+
ggplot(shift_beta_l %>%
         filter(parameter == 'm_shift'),
       aes(s_u, o_u, fill = shift))+
  scale_fill_gradient(low = 'white', high = 'black',
                      name = expression(paste(Delta, theta[ppt]^mu)))+
  geom_tile() +
  labs(y = expression(paste('Other Unc.')),
       x = expression(paste('Self Unc.'))
       ) +
  scale_y_discrete(breaks = c(1, 5, 10, 15)) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  theme_bw(base_size=18)


## Load fitted parameters --------------------------------------------------

#Fitted hierarchical Parameters
BPD_full      <- readMat('FittedModelFiles/hbi_g1_BPD.mat')
CON_full      <- readMat('FittedModelFiles/hbi_g2_BPD.mat')

#Simulated choices with best fitting model
BPD_sim       <- readMat('Data_Simulated/data_sim_g1.mat')
CON_sim       <- readMat('Data_Simulated/data_sim_g2.mat')

BPD_full_sim  <- readMat('Data_Simulated/full_sim_g1.mat')
CON_full_sim  <- readMat('Data_Simulated/full_sim_g2.mat')

#Recovered models
BPD_hbi_recov <- readMat('FittedModelFiles/hbi_g1_BPD_recovery.mat')
CON_hbi_recov <- readMat('FittedModelFiles/hbi_g2_BPD_recovery.mat')

## Load & transform individual parameters ----------------------------------------------

ID1 <- beh_anal[beh_anal$group=='BPD',]$ID
ID2 <- beh_anal[beh_anal$group=='CON',]$ID

BPD_parms <- extract_parameters(BPD_full, ID1)
CON_parms <- extract_parameters(CON_full, ID2)

joint_parms <- rbind(BPD_parms, CON_parms) %>%
               plyr::join(., intent_dat %>%
                            filter(Phase ==2) %>%
                            dplyr::select(ID, SI:A_ETS_Credulity, -RT), by = 'ID') %>%
  distinct()

## Responsibility ----------------------------------------------------------

mod_col <- c('#CCDBDC', '#2BD8FF', '#80CED7', '#020100', '#A81ADB')

respon <- BPD_full$cbm[,,1]$output[,,1]$responsibility %>%
  as.data.frame()%>%
  dplyr::select(1:5)%>%
  mutate(ID = as.vector(ID1), group = 'BPD', ID_n = 1:length(ID)) %>%
  rbind(.,CON_full$cbm[,,1]$output[,,1]$responsibility %>%
          as.data.frame()%>%
          dplyr::select(1:5)%>%
          mutate(ID = as.vector(ID2), group = 'CON', ID_n = 1:length(ID))) %>%
  mutate(ID_n = c(1:50, 1:53)) %>%
  dplyr::rename('M2'=2, 'M3'=4, 'M1'=3, 'M4'=5, 'Beta'=1) %>%
  pivot_longer(1:5, names_to = 'Model', values_to = 'Respon')

freq <- BPD_full$cbm[,,1]$output[,,1]$model.frequency %>%
  as.data.frame()%>%
  dplyr::select(1:5)%>%
  mutate(group = 'BPD') %>%
  rbind(.,CON_full$cbm[,,1]$output[,,1]$model.frequency %>%
          as.data.frame()%>%
          dplyr::select(1:5)%>%
          mutate(group = 'CON')) %>%
  dplyr::rename('M2'=2, 'M3'=4, 'M1'=3, 'M4'=5, 'Beta'=1) %>%
  pivot_longer(1:5, names_to = 'Model', values_to = 'Freq')

exprob <-rbind(BPD_full$cbm[,,1]$output[,,1]$exceedance.prob,
               CON_full$cbm[,,1]$output[,,1]$exceedance.prob) %>%
          as.data.frame()%>%
  dplyr::select(1:5)%>%
  mutate(group = c('BPD', 'CON')) %>%
  dplyr::rename('M2'=2, 'M3'=4, 'M1'=3, 'M4'=5, 'Beta'=1) %>%
  pivot_longer(1:5, names_to = 'Model', values_to = 'ExProb')

model_comparison_outcomes <- respon %>%
  plyr::join(., freq, by = c('group', 'Model')) %>%
  plyr::join(., exprob, by = c('group', 'Model'))

library(patchwork)

rplot <- ggplot(respon,
       aes(ID_n, Respon, fill = Model))+
  geom_col(colour = 'black')+
  facet_wrap(~group, scales = 'free_x')+
  scale_fill_manual(values = mod_col)+
  labs(x = 'ID', y = 'Respon.')+
  coord_cartesian(ylim = c(0,1))+
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50), expand = c(0,0))+
  scale_y_continuous(breaks = c(0.01, 0.5, 1), labels = c('0','0.5','1'), expand = c(0,0))+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'none',
        axis.text.x = element_blank())

fplot <- ggplot(freq,
       aes(Model, ifelse(Freq<0.01, 0+0.01, Freq), fill = Model))+
  geom_col(position = 'dodge', width = 0.5, colour = 'black')+
  labs(x = 'Model', y = 'Frequency')+
  scale_fill_manual(values=mod_col)+
  coord_cartesian(ylim = c(0,1))+
  facet_wrap(~group)+
  scale_y_continuous(breaks = c(0.01, 0.5, 1), labels = c('0','0.5','1'), expand = c(0,0))+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'none')

explot <- ggplot(exprob,
       aes(Model, ifelse(ExProb<0.01, 0+0.01, ExProb), fill = Model))+
  geom_col(position = 'dodge', width = 0.5, colour = 'black')+
  scale_fill_manual(values=mod_col,
                    name='Model')+
  facet_wrap(~group)+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = c(0.01, 0.5, 1), labels = c('0','0.5','1'), expand = c(0,0))+
  labs(y = 'Ex. Prob')+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'top')

comp_plot <- explot/fplot/rplot &
  theme(plot.margin = margin(rep(1, 4), unit = 'mm'),
        text = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
  )
comp_plot

## Recovery ------------------------------------------------

### Model Recovery ------------------------------------------------

respon_r <- BPD_hbi_recov$cbm[,,1]$output[,,1]$responsibility %>%
  as.data.frame()%>%
  dplyr::select(1:5)%>%
  mutate(ID = as.vector(ID1), Group = 'BPD', ID_n = 1:length(ID)) %>%
  rbind(.,CON_hbi_recov$cbm[,,1]$output[,,1]$responsibility %>%
          as.data.frame()%>%
          dplyr::select(1:5)%>%
          mutate(ID = as.vector(ID2), Group = 'CON', ID_n = 1:length(ID))) %>%
  mutate(ID_n = c(1:50, 1:53)) %>%
  rename('M2'=2, 'M3'=4, 'M1'=3, 'M4'=5, 'Beta'=1) %>%
  pivot_longer(1:5, names_to = 'Model', values_to = 'Respon')

freq_r <- BPD_hbi_recov$cbm[,,1]$output[,,1]$model.frequency %>%
  as.data.frame()%>%
  dplyr::select(1:5)%>%
  mutate(Group = 'BPD') %>%
  rbind(.,CON_hbi_recov$cbm[,,1]$output[,,1]$model.frequency %>%
          as.data.frame()%>%
          dplyr::select(1:5)%>%
          mutate(Group = 'CON')) %>%
  rename('M2'=2, 'M3'=4, 'M1'=3, 'M4'=5, 'Beta'=1) %>%
  pivot_longer(1:5, names_to = 'Model', values_to = 'Freq')

exprob_r <-rbind(BPD_hbi_recov$cbm[,,1]$output[,,1]$exceedance.prob,
               CON_hbi_recov$cbm[,,1]$output[,,1]$exceedance.prob) %>%
          as.data.frame()%>%
  dplyr::select(1:5)%>%
  mutate(Group = c('BPD', 'CON')) %>%
  rename('M2'=2, 'M3'=4, 'M1'=3, 'M4'=5, 'Beta'=1) %>%
  pivot_longer(1:5, names_to = 'Model', values_to = 'ExProb')

library(patchwork)

rplot_r <- ggplot(respon_r,
       aes(ID_n, Respon, fill = Model))+
  geom_col(colour = 'black')+
  facet_wrap(~Group, scales = 'free_x')+
  scale_fill_manual(values = mod_col)+
  labs(x = 'ID', y = 'Respon')+
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50), expand = c(0,0))+
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0))+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fplot_r <- ggplot(freq_r,
       aes(Model, Freq, fill = Model))+
  geom_col(position = 'dodge', width = 0.5, colour = 'black')+
  labs(x = 'Model', y = 'Frequency')+
  scale_fill_manual(values = mod_col)+
  facet_wrap(~Group)+
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0), limits=c(0,1))+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'none')

explot_r <- ggplot(exprob_r,
       aes(Model, ifelse(ExProb<0.01, 0+0.01, ExProb), fill = Model))+
  geom_col(position = 'dodge', width = 0.5, colour = 'black')+
  scale_fill_manual(values = mod_col)+
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0), limits=c(0,1))+
  facet_wrap(~Group)+
  labs(y = 'Ex. Prob')+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'top')

model_rec <- explot_r/fplot_r/rplot_r &
  theme(plot.margin = margin(rep(1, 4), unit = 'mm'),
        text = element_text(size = 20),
        axis.title.x = element_blank()
  )

model_rec

### Model Confusion --------------------------

confus <- cbind(
  rbind(BPD_full$cbm[,,1]$output[,,1]$responsibility,
        CON_full$cbm[,,1]$output[,,1]$responsibility),
  rbind(BPD_hbi_recov$cbm[,,1]$output[,,1]$responsibility,
        CON_hbi_recov$cbm[,,1]$output[,,1]$responsibility)
      )

confus_r   <- rcorr(confus, type = 'spearman')
confus_r$r <- confus_r$r[1:5, 6:10]
confus_r$P <- confus_r$P[1:5, 6:10]

colnames(confus_r$r) <- c('Beta', 'M1', 'M2', 'M3', 'M4')
rownames(confus_r$r) <- c('Beta', 'M1', 'M2', 'M3', 'M4')

ggcorrplot::ggcorrplot(confus_r$r,
                       hc.order = F,
   outline.col = "white",
   colors = c('#291720','white', '#007A5A'),
   p.mat = confus_r$P,
   lab = T,
   legend.title = 'Pearson\nR',
   insig='blank',
   sig.level = 0.01) +
  labs(x = 'Recovered', y = 'Real')+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme_bw()+
  theme(axis.title = element_text(color = 'black', size = 16),
        axis.text  = element_text(size=15),
        panel.grid = element_blank(),
        legend.position = 'none')

### Parameter Recovery --------------------------

BPD_hbi_recov$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]]
CON_hbi_recov$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]]

rec_BPD_parms <- as.data.frame(matrix(NA, nrow = nrow(BPD_hbi_recov$cbm[,,1]$output[,,1]$responsibility), 9))
rec_CON_parms <- as.data.frame(matrix(NA, nrow = nrow(CON_hbi_recov$cbm[,,1]$output[,,1]$responsibility), 7))

rec_BPD_parms[,1:8] <- BPD_hbi_recov$cbm[,,1]$output[,,1]$parameters[5,][[1]][[1]]
rec_BPD_parms[,9]   <- 'BPD_r'
rec_CON_parms[,1:6] <- CON_hbi_recov$cbm[,,1]$output[,,1]$parameters[3,][[1]][[1]]
rec_CON_parms[,7]   <- 'CON_r'

rec_joint_parms <- data.frame(ID      = c(ID1, ID2),
                           rec_alpha   = c((1/(1+exp(-rec_BPD_parms[,1])))*30,
                                           (1/(1+exp(-rec_CON_parms[,1])))*30),
                           rec_beta    = c(rec_BPD_parms[,2],
                                           rec_CON_parms[,2]),
                           rec_alpha_v= c(exp(rec_BPD_parms[,3]),
                                          exp(rec_CON_parms[,3])),
                           rec_beta_v = c(exp(rec_BPD_parms[,4]),
                                          exp(rec_CON_parms[,4])),
                           rec_alpha_par= c((1/(1+exp(-rec_BPD_parms[,5])))*30,
                                            rep(NA,53)),
                           rec_beta_par = c(rec_BPD_parms[,6],
                                            rep(NA,53)),
                           rec_alpha_ref= c(exp(rec_BPD_parms[,7]),
                                            exp(rec_CON_parms[,5])),
                           rec_beta_ref = c(exp(rec_BPD_parms[,8]),
                                        exp(rec_CON_parms[,6])),
                           group   = c(rec_BPD_parms[,9],
                                       rec_CON_parms[,7])) %>%
  distinct()


recovered_ps <- plyr::join(rec_joint_parms %>% dplyr::select(-group),
                           joint_parms %>% dplyr::select(ID, alpha:group),
                           by = 'ID')

rec_cor_parms  <- rcorr(as.matrix(recovered_ps %>%
                                       dplyr::select(c(rec_alpha:beta_ref,
                                                       -group,
                                                       -rec_alpha_par,-rec_beta_par,
                                                       -alpha_par, -beta_par))
                                   ))

ggcorrplot::ggcorrplot(rec_cor_parms$r[1:6, 7:12],
                       p.mat = rec_cor_parms$P[1:6, 7:12],
                       lab = T,
                       colors = c( '#291720','white','#007A5A'),
                       legend.title = 'Pearson\nR', sig.level = 0.01,insig='blank',
                       pch.cex = 15,
                       pch.col = 1) +
  scale_x_discrete(labels = c(
    expression(paste(alpha['ppt']^m)),
    expression(paste(beta['ppt']^m)),
    expression(paste(alpha['ppt']^sigma)),
    expression(paste(beta['ppt']^sigma)),
    expression(paste(alpha['par']^ref)),
    expression(paste(beta['par']^ref))
  ),
  expand=c(0,0)) +
    scale_y_discrete(labels = c(
    expression(paste(alpha['ppt']^m)),
    expression(paste(beta['ppt']^m)),
    expression(paste(alpha['ppt']^sigma)),
    expression(paste(beta['ppt']^sigma)),
    expression(paste(alpha['par']^ref)),
    expression(paste(beta['par']^ref))
  ),
  expand=c(0,0)) +
  labs(x = 'Recovered', y = 'Real')+
  theme_bw()+
  theme(axis.title = element_text(color = 'black', size = 16),
        axis.text  = element_text(size=15),
        legend.position = 'none')

## Generative Ability --------------------------------------------------------

#### Congruency of matching between participant and partner (SimAFix/SimA) ------------------------------------------------

n_vectorsBPD        <- length(BPD_full_sim$results)
result_matrix_BPDFIX<- matrix(NA, nrow = 54, ncol = n_vectorsBPD)
result_matrix_BPDACT<- matrix(NA, nrow = 54, ncol = n_vectorsBPD)
n_vectorsCON        <- length(CON_full_sim$results)
result_matrix_CONFIX<- matrix(NA, nrow = 54, ncol = n_vectorsCON)
result_matrix_CONACT<- matrix(NA, nrow = 54, ncol = n_vectorsCON)

for(i in 1:n_vectorsBPD){
  result_matrix_BPDFIX[,i] <- BPD_full_sim$results[i,][[1]][[1]][,,1]$simAFix[BPD_full_sim$results[i,][[1]][[1]][,,1]$simAFix!=0]
  result_matrix_BPDACT[,i] <- BPD_full_sim$results[i,][[1]][[1]][,,1]$simA[37:90]
}
for(i in 1:n_vectorsCON){
  result_matrix_CONFIX[,i] <- CON_full_sim$results[i,][[1]][[1]][,,1]$simAFix[CON_full_sim$results[i,][[1]][[1]][,,1]$simAFix!=0]
  result_matrix_CONACT[,i] <- CON_full_sim$results[i,][[1]][[1]][,,1]$simA[37:90]
}

FixCong <- cbind(result_matrix_BPDFIX, result_matrix_CONFIX)
ActCong <- cbind(result_matrix_BPDACT, result_matrix_CONACT)

PPT_PARcong <- matrix(NA, nrow = 54, ncol = ncol(FixCong))

for (i in 1:ncol(FixCong)){
  match_ppt_par  <- ifelse(FixCong[,i]==ActCong[,i], 1, 0)
  PPT_PARcong[,i]<- match_ppt_par
}

PPT_PARcong <- as.data.frame(PPT_PARcong)
names(PPT_PARcong)[1:50]  <- ID1$ID
names(PPT_PARcong)[51:103] <- ID2$ID

PPT_PARcong <- PPT_PARcong %>%
  pivot_longer(1:103, names_to = 'ID', values_to = 'match_ppt_par') %>%
  group_by(ID) %>%
  mutate(matched_tot = sum(match_ppt_par)/54,
         Group = ifelse(nchar(ID)<10, 'BPD', 'CON')) %>%
  dplyr::select(ID, matched_tot, Group) %>%
  distinct()

ppt_par_congp1 <- ggplot(PPT_PARcong, aes(matched_tot, Group, fill = Group))+
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, height = 0.2)+
  geom_boxplot(width = 0.2, alpha = 0.7)+
  coord_cartesian(x = c(0, 1))+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, '.25','.5', '.75', '1'),
                     expand = c(0,0))+
  scale_fill_manual(values=colour_group)+
  labs(x = 'PPT-PAR Similarity')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

#### LL ------------------------------------------------

LL <- matrix(NA, nrow = 103, ncol = 6)
LL[,1] <- rbind(BPD_full_sim$F,CON_full_sim$F)
LL[,2] <- c(extract_and_combine(BPD_full_sim$results, 'lik1', 1), extract_and_combine(CON_full_sim$results, 'lik1', 1))
LL[,3] <- c(extract_and_combine(BPD_full_sim$results, 'lik2', 1), extract_and_combine(CON_full_sim$results, 'lik2', 1))
LL[,4] <- c(extract_and_combine(BPD_full_sim$results, 'lik3', 1), extract_and_combine(CON_full_sim$results, 'lik3', 1))
LL[,5] <- as.matrix(rbind(intent_dat %>% filter(nchar(ID)<10) %>% dplyr::select(ID) %>% unique(),
                          intent_dat %>% filter(!nchar(ID)<10) %>% dplyr::select(ID) %>% unique()))
LL[,6] <- ifelse(nchar(LL[,5])<10, 'BPD', 'CON')

LL  <- as.data.frame(LL) %>% mutate(across(1:4, as.numeric))
names(LL) <- c('LL', 'lik1', 'lik2', 'lik3', 'ID', 'group')

#### Probabilities ------------------------------------------------

prob1 <- rbind(extract_and_combine(BPD_full_sim$results, 'prob1', 36), extract_and_combine(CON_full_sim$results, 'prob1', 36))
prob2 <- rbind(extract_and_combine(BPD_full_sim$results, 'prob2', 54), extract_and_combine(CON_full_sim$results, 'prob2', 54))
prob3 <- rbind(extract_and_combine(BPD_full_sim$results, 'prob3', 36), extract_and_combine(CON_full_sim$results, 'prob3', 36))

probabilities <- matrix(NA, nrow = 126, ncol = nrow(prob1))

for (i in 1:nrow(prob1)){
  probabilities[1:36, i]   <- prob1[i,]
  probabilities[37:90, i]  <- prob2[i,]
  probabilities[91:126, i] <- prob3[i,]
}

probabilities <- as.data.frame(probabilities)
names(probabilities)[1:50]  <- ID1
names(probabilities)[51:103] <- ID2

probabilities <- probabilities %>%
  pivot_longer(1:103, names_to = 'ID', values_to = 'Probs') %>%
  group_by(ID) %>%
  mutate(Trial = 1:126,
         Group = ifelse(nchar(ID)<10, 'BPD', 'CON'),
         Phase = ifelse(Trial %in% 1:36, 1, ifelse(Trial %in% 37:90, 2, 3))) %>%
  arrange(ID, Trial)

#### Check action congruency ------------------------------------------------

sim_actions      <- BPD_sim$data.sim[,1][[1]][[1]]
sim_actions      <- as_tibble(sim_actions)

for (i in 2:length(unique(intent_dat$ID))){
  if(i %in% 2:50){
  sim_actions         <- rbind(sim_actions, as.data.frame(BPD_sim$data.sim[,i][[1]][[1]]))
  }
  if(i %in% 51:103){
  sim_actions         <- rbind(sim_actions, as.data.frame(CON_sim$data.sim[,i-50][[1]][[1]]))
  }
}

sim_actions[,1]              <- rbind(intent_dat %>% filter(nchar(ID)<10) %>% dplyr::select(ID),
                                      intent_dat %>% filter(!nchar(ID)<10) %>% dplyr::select(ID))
names(sim_actions)[c(1:2,7)] <- c('ID', 'trial', 'sim_choice')
sim_actions <- sim_actions %>% dplyr::select(1,2,7)
sim_actions <- sim_actions %>% group_by(ID) %>% mutate(trial = 1:126)
intent_dat_126_check <- intent_dat %>% group_by(ID) %>% mutate(trial = 1:126)

action_cong <- plyr::join(sim_actions, intent_dat_126_check, by = c('ID', 'trial')) %>%
  mutate(cong = ifelse(sim_choice == choice, 1, 0)) %>%
  group_by(ID) %>%
  mutate(sumCong = sum(cong)/126,
         group = ifelse(nchar(ID)<10, 'BPD', 'CON')) %>%
  group_by(ID, Phase) %>%
  mutate(sumPhaseCong = ifelse(Phase %in% c(1,3), sum(cong)/36, sum(cong)/54))

#### Plot all ------------------------------------------------

probp1 <- ggplot(probabilities, aes(Probs, fill = Group))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=colour_group)+
  labs(x = expression(paste('p(Choice | ', theta, ')')))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(legend.position = 'none',
        text = element_text(size = 16),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

congphase <- ggplot(action_cong %>%
         dplyr::select(ID, sumPhaseCong, group, Phase) %>% distinct(),
       aes(sumPhaseCong, group, fill = group))+
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, height = 0.2)+
  geom_boxplot(width = 0.2, alpha = 0.7)+
  coord_cartesian(x = c(0.5, 1))+
  scale_x_continuous(breaks = c(0.5, 0.75, 1), labels = c('.5', '.75', '1'), expand = c(0,0))+
  #scale_y_continuous(expand=c(0,0))+
  facet_wrap(~Phase, scales = 'free_x')+
  labs(x = 'Model Accuracy by Phase')+
  theme_bw()+
  theme(legend.position = 'none')

congall <- ggplot(action_cong %>%
                  ungroup() %>%
                  dplyr::select(ID, sumCong, group) %>%
                  distinct(),
       aes(sumCong, group, fill = group))+
  #geom_density(alpha = 0.5)+
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, height = 0.2)+
  geom_boxplot(width = 0.2, alpha = 0.7)+
  geom_vline(aes(xintercept = mean(sumCong)), size = 1.2)+
  coord_cartesian(x = c(0.5, 1))+
  geom_vline(xintercept = c(0.5, 1), colour = c('grey', 'darkgreen'), size = 3)+
  coord_cartesian(xlim = c(0.5,1))+
  scale_x_continuous(breaks = c(seq(0.5, 1, 0.1)), labels = c('.5', '.6', '.7', '.8', '.9', '1'), expand = c(0,0))+
  labs(x = 'Overall Model Accuracy')+
  theme_bw()+
  theme(legend.position = c(0.15, 0.65),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black'))

LLall <- ggplot(LL, aes(LL, fill = group))+
  geom_density(alpha = 0.5)+
  geom_vline(aes(xintercept = mean(LL)), size = 1.2)+
  geom_vline(xintercept = c(log(0.5)*126, 0), colour = c('grey', 'darkgreen'), size = 3)+
  scale_fill_manual(values=colour_group)+
  scale_x_continuous(expand = c(0,0), limits = c(-90, 1))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = 'Total LL Overall')+
  theme_bw()+
  theme(legend.position = 'none',
        text = element_text(size = 16),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

Accp1 <- ((ppt_par_congp1 | congall)/congphase) &
  scale_fill_manual(values=colour_group) &
  theme(text = element_text(size = 16),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

Accp1/(LLall | probp1)

#### Check behavioural comparison --------------------------------------------
sumCbpd <- rep(NA, length(BPD_sim$data.sim))
sumCcon <- rep(NA, length(CON_sim$data.sim))

for(i in 1:length(BPD_sim$data.sim)){
 y           <- as.data.frame(BPD_sim$data.sim[,i][[1]][[1]])
 y           <- y[y$V9==2,]
 y$cor       <- ifelse(y$V7==y$V8, 1, 0)
 sumCbpd[i]  <- sum(y[,10])
 if(i>2){predCorbpd <- rbind(predCorbpd, y[,c(1,2,10)])}else{predCorbpd <- y[,c(1,2,10)]}
}

for(i in 1:length(CON_sim$data.sim)){
 y           <- as.data.frame(CON_sim$data.sim[,i][[1]][[1]])
 y           <- y[y$V9==2,]
 y$cor       <- ifelse(y$V7==y$V8, 1, 0)
 sumCcon[i]  <- sum(y[,10])
 if(i>2){predCorcon <- rbind(predCorcon, y[,c(1,2,10)])}else{predCorcon <- y[,c(1,2,10)]}
}

BPD_curveSim <- glm(cor~V2,
    data = predCorbpd,
    family = "binomial")
CON_curveSim <- glm(cor~V2,
    data = predCorcon,
    family = "binomial")

pred_time_datSim <- data.frame(V2=seq(0, 54, 0.01))

BPD_probSim <- BPD_curveSim %>% predict(pred_time_datSim, type = "response")
CON_probSim <- CON_curveSim %>% predict(pred_time_datSim, type = "response")

pred_time_datSim$BPD_probSim = BPD_probSim
pred_time_datSim$CON_probSim = CON_probSim
pred_time_datSim <- pred_time_datSim %>%
  pivot_longer(2:3, names_to = 'Group', values_to = 'p(Correct)') %>%
  mutate(Group=ifelse(Group=='BPD_probSim', 'BPD', 'CON'))

ggplot(pred_time_datSim,
       aes(V2, `p(Correct)`, colour = Group))+
  geom_smooth(se = T, size = 1) +  # Smoother lines without confidence intervals
  scale_color_manual(values=colour_group)+
  labs(x = 'Trial', y = 'P(Correct)') +  # Proper axis labels
  theme_bw() +
  theme(
    text = element_text(size = 16),
    legend.position = 'none',
    legend.background = element_rect(colour = 'black', fill = 'white'),
    legend.title = element_blank(),  # Improve legend title readability
    legend.text = element_text(size = 12)  # Improve legend text readability
  )+
ggplot(data.frame(sumCor=c(sumCbpd,sumCcon)/54, group=c(rep('BPD', 50), rep('CON', 53))),
       aes(group, sumCor, fill = group))+
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, height = 0.2, size = 2) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0.5, size = 1.2)+
  scale_fill_manual(values=colour_group)+
  scale_y_continuous(expand=c(0,0))+
  stat_compare_means(label.y = 0.25, size = 5)+
  coord_cartesian(y = c(0, 1))+
  labs(x = 'Group', y = '% Correct')+
  theme_bw()+
  theme(
      legend.position = 'none',
      text = element_text(size = 18),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5)
    )

#Plot together:
library(ggpattern)
ModCorSumComp <- rbind(psych_curves %>%
         dplyr::select(group, corSum) %>%
         distinct() %>%
           mutate(type='Obs'),
      data.frame(corSum=c(sumCbpd,sumCcon)/54, group=c(rep('BPD', 50), rep('CON', 53))) %>%
        mutate(type='Model'))

ggplot(rbind(pred_time_datSim %>% mutate(type='Model') %>% rename(trial=V2),
             pred_time_dat    %>% mutate(type='Obs.')),
       aes(trial, `p(Correct)`, colour = Group))+
  geom_smooth(se = T, size = 1, aes(linetype=type)) +  # Smoother lines without confidence intervals
  scale_color_manual(values=colour_group)+
  labs(x = 'Trial', y = 'P(Correct)') +  # Proper axis labels
  theme_bw() +
  theme(
    text = element_text(size = 16),
    legend.background = element_rect(colour = 'black', fill = 'white'),
    legend.title = element_blank(),  # Improve legend title readability
    legend.text = element_text(size = 12)  # Improve legend text readability
  )+
ggplot(ModCorSumComp,
       aes(group, corSum, fill = group, pattern = type, shape=type, alpha = type))+
  geom_jitter(alpha = 0.1, width = 0.1, height = 0.2, size = 2) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  #geom_boxplot_pattern() +
  geom_hline(yintercept = 0.5, size = 1.2)+
  scale_fill_manual(values=colour_group)+
  scale_alpha_manual(name = "type", values = c(1, 0.1))+
  scale_shape_manual(values=c(21, 22))+
  scale_y_continuous(expand=c(0,0))+
  stat_compare_means(label.y = 0.25, size = 5, label='p.signif')+
  coord_cartesian(y = c(0, 1))+
  labs(x = 'Group', y = '% Correct')+
  theme_bw()+
  theme(
      legend.position = 'none',
      text = element_text(size = 18),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5)
    )+

ggplot(plyr::join(rbind(psych_curves %>%
       dplyr::select(group, corSum) %>%
       distinct() %>%
       rename(Obs=corSum)),
       data.frame(Model=c(sumCbpd,sumCcon)/54, ID = c(ID1, ID2)),
       by = 'ID'),
       aes(Obs, Model, fill = group))+
  geom_point(shape=21)+
  geom_smooth(method='lm', colour = 'black')+
  stat_cor(aes(label = ..p.label.., colour = group), label.x = c(0.5, 0.5), label.y = c(0.25, 0.35))+
  scale_fill_manual(values=colour_group)+
  scale_colour_manual(values=colour_group)+
  scale_shape_manual(values=c(21, 22))+
  scale_y_continuous(expand=c(0,0), limits=c(0.25,1), breaks=c(0, 0.5, 1))+
  scale_x_continuous(expand=c(0,0), limits=c(0.25,1), breaks=c(0, 0.5, 1))+
  coord_cartesian(y = c(0, 1))+
  labs(x = 'Obs.', y = 'Model')+
  theme_bw()+
  theme(
      legend.position = 'none',
      text = element_text(size = 18),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5)
    )

### Compare with choice behaviour -----

choice_beh_corr <- plyr::join(summary_choices_ppts_1 %>% dplyr::select(ID, Ic, Pc, Pi), joint_parms_ex, by = 'ID')

choice_par <- choice_beh_corr %>%
  ungroup() %>%
  dplyr::select(Pc, Ic, Pi, beta, alpha) %>%
  as.data.frame() %>%
  as.matrix() %>%
  pcor()

ggcorrplot::ggcorrplot(choice_par$estimate[1:3, 4:5],
                       p.mat = choice_par$p.value[1:3, 4:5],
                       lab = T,
                       sig.level = 0.01,
                       insig = 'blank'
                        )

ggplot(choice_beh_corr, aes(Pc, beta))+
ggplot(choice_beh_corr, aes(Ic, alpha))+
ggplot(choice_beh_corr, aes(Pi, alpha))+
ggplot(choice_beh_corr, aes(Ic, beta))&
  geom_point()&
  geom_smooth(method = 'lm') &
  stat_cor()

## Calculate Belief Distributions -----------------------------------------------------------

lim = 30
res = 0.25

for (i in 1:length(ID1)){
  xdist <-   data.frame(
    bP1  = as.vector(BPD_full_sim$results[i,][[1]][[1]][,,1]$beta.marg1),
    bP2a = as.vector(BPD_full_sim$results[i,][[1]][[1]][,,1]$beta.marg2a),
    bP2b = as.vector(BPD_full_sim$results[i,][[1]][[1]][,,1]$beta.marg2b),
    bP3  = as.vector(BPD_full_sim$results[i,][[1]][[1]][,,1]$beta.marg3),
    aP1  = as.vector(BPD_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg1),
    aP2a = as.vector(BPD_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg2a),
    aP2b = as.vector(BPD_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg2b),
    aP3  = as.vector(BPD_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg3),
    beta = seq(-lim, lim, res),
    alpha = seq(0, lim, res/2),
    ID = ID1[i],
    group = 'BPD',
    Model = 'M4'
  )
  if(i>1){BPDdist <- rbind(BPDdist, xdist)} else {BPDdist <- xdist}
}
for (i in 1:length(ID2)){
  zdist <-   data.frame(
  bP1   = as.vector(CON_full_sim$results[i,][[1]][[1]][,,1]$beta.marg1),
  bP2a  = as.vector(CON_full_sim$results[i,][[1]][[1]][,,1]$beta.marg2a),
  bP2b  = as.vector(CON_full_sim$results[i,][[1]][[1]][,,1]$beta.marg2b),
  bP3   = as.vector(CON_full_sim$results[i,][[1]][[1]][,,1]$beta.marg3),
  aP1   = as.vector(CON_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg1),
  aP2a  = as.vector(CON_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg2a),
  aP2b  = as.vector(CON_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg2b),
  aP3   = as.vector(CON_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg3),
  beta  = seq(-lim, lim, res),
  alpha = seq(0, lim, res/2),
  ID = ID2[i],
  group = 'CON',
  Model = 'M1'
  )
  if(i>1){CONdist <- rbind(CONdist, zdist)} else {CONdist <- zdist}
}

int_dist <- rbind(BPDdist, CONdist)

## Individual Distribution Examples For Each Model -----------

BPDplotindiv <- ggplot(int_dist %>% filter(ID == ID1[2], group == 'BPD'))+
  geom_line(aes(beta, bP1, colour = '1 - P1'), size = 1.2)+
  geom_line(aes(beta, bP2a, colour = '2A - Prior'), size = 1.2)+
  geom_line(aes(beta, bP2b, colour = '2B - Post.'), size = 1.2)+
  #geom_vline(xintercept = intent_dat[intent_dat$ID=='PDA130',]$server_beta_par[1],
  #           linetype = 2, colour = 'grey')+
  labs(title = '', x = expression(paste(beta)), y = expression(paste('p(',beta,')')))+
  scale_colour_manual(values = c('black','#FBD19D', '#F7941D', 'grey'))+
  theme(legend.position = 'none')
CONplotindiv <- ggplot(int_dist %>% filter(ID == ID2[17], group == 'CON'))+
  geom_line(aes(beta, bP1, colour = '1 - P1'), size = 1.2)+
  geom_line(aes(beta, bP2a, colour = '2A - Prior'), size = 1.2)+
  geom_line(aes(beta, bP2b, colour = '2B - Post.'), size = 1.2)+
  geom_line(aes(beta, bP3, colour = '3 - P3'), size = 1.2)+
  #geom_vline(xintercept = intent_dat[intent_dat$ID=='PD061123SJ',]$server_beta_par[1],
  #           linetype = 2, colour = 'grey')+
  labs(title = '', x = expression(paste(beta)), y = expression(paste('p(',beta,')')))+
  scale_colour_manual(values = c('black','#FBD19D', '#F7941D', 'grey'))+
  theme(legend.position = 'none')

(CONplotindiv/
  BPDplotindiv) &
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin=margin(1,1,1,1),
        axis.title.y =element_text(vjust=-10))

## Parameter analysis ----------------------------------------------------------

#Calculate final distribution at end of Phase 2

par_parms <- int_dist %>%
  group_by(ID) %>%
  mutate(beta_par_m = sum(beta*bP2b),
         alpha_par_m = sum(alpha*aP2b),
         beta_par_pri = sum(beta*bP2a),
         alpha_par_pri = sum(alpha*aP2a),
         shift_beta_p2 = abs(beta_par_pri-beta_par_m),
         shift_alpha_p2 = abs(alpha_par_pri-alpha_par_m),
         beta_hat = sum(beta*bP3),
         alpha_hat= sum(alpha*aP3),
         beta_hat_sd = sqrt(sum((beta-beta_hat)^2 * bP3)),
         alpha_hat_sd = sqrt(sum((alpha-alpha_hat)^2 * aP3))) %>%
  dplyr::select(ID, group,
                beta_par_m, alpha_par_m,
                beta_par_pri, alpha_par_pri,
                beta_hat, alpha_hat,
                beta_hat_sd, alpha_hat_sd,
                shift_beta_p2, shift_alpha_p2) %>%
  distinct() %>%
  plyr::join(., intent_dat %>%
         dplyr::select(ID, server_beta_par, server_alpha_par)) %>%
  mutate(disp_beta_par = abs(server_beta_par - beta_par_m),
         disp_alpha_par = abs(server_alpha_par - alpha_par_m)) %>%
  dplyr::select(ID,
                beta_par_m, alpha_par_m,
                beta_par_pri, alpha_par_pri,
                disp_beta_par, disp_alpha_par,
                beta_hat_sd, alpha_hat_sd,
                beta_hat, alpha_hat,
                shift_beta_p2, shift_alpha_p2) %>%
  distinct()

joint_parms_ex <- plyr::join(par_parms, joint_parms, by = 'ID') %>%
  mutate(disp_beta_ppt = abs(beta_hat - beta),
         disp_alpha_ppt = abs(alpha_hat - alpha),
         disp_beta_hat_sd = abs(beta_hat_sd - beta_v),
         disp_alpha_hat_sd = abs(alpha_hat_sd - alpha_v),
         distancebeta=abs(server_beta_par-server_beta_ppt),
         distancealpha=abs(server_alpha_par-server_alpha_ppt))


### Partner median priors in BPD ---------------------------------------------------

beta_bpd <- joint_parms_ex %>% filter(group=='BPD_full')

ggplot(beta_bpd,
       aes(beta_par))+
  geom_density(aes(beta),
               fill=colour_group[1], alpha = 0.2)+
  geom_density(fill=colour_group[1], size=1)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(x=expression(paste(beta['par']^m)),
       y=expression(paste('p(',beta['par']^m, ')')))+
  theme_bw(base_size=24)+
  theme(panel.grid = element_blank())+

ggplot(beta_bpd,
       aes(alpha_par))+
  geom_density(aes(alpha),
               fill=colour_group[1], alpha = 0.2)+
  geom_density(data=beta_bpd, fill=colour_group[1], size=1)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(x=expression(paste(alpha['par']^m)),
       y=expression(paste('p(',alpha['par']^m, ')')))+
  theme_bw(base_size=24)+
  theme(panel.grid = element_blank())

bayes.t.test(beta_bpd$beta,
             beta_bpd$beta_par,
             paired = T)

bayes.t.test(beta_bpd$beta_par, mu = 0)

bayes.t.test(beta_bpd$alpha,
             beta_bpd$alpha_par,
             paired = T)

### Group Differences -----------------

bayes_list <- c('alpha',
                'beta',
                'alpha_v',
                'beta_v',
                'alpha_ref',
                'beta_ref',
                'alpha_hat',
                'beta_hat',
                'disp_alpha_par',
                'disp_beta_par',
                'shift_beta_p2',
                'shift_alpha_p2'
                )

library(BayesianFirstAid)
for(i in bayes_list){

  x_parms <- joint_parms_ex

  print(paste('Now running ', i, sep = ''))

  c1 <- bayes.t.test(x_parms[,i][x_parms$group=='BPD_full',],
                     x_parms[,i][x_parms$group=='CON_full',],
                     paired = F)

  x1 <- as.data.frame(t(c1$stats[5,])) %>%
        mutate(parameter = i,
               sig = ifelse((HDIlo > 0 & HDIup > 0) | (HDIlo < 0 & HDIup < 0), 'Yes', 'No'),
               eff = c1$stats[5,1],
               effl= c1$stats[5,5],
               effh= c1$stats[5,6],
               eff_size=c1$stats[8,1],
               eff_sizel= c1$stats[8,5],
               eff_sizeh= c1$stats[8,6],
               sd_diff = c1$stats[6,1],
               sdHDIlo = c1$stats[6,5],
               sdHDIup = c1$stats[6,6],
               sigSD = ifelse((sdHDIlo > 0 & sdHDIup > 0) | (sdHDIlo < 0 & sdHDIup < 0), 'Yes', 'No'))
  x3 <- rbind(c1$mcmc_samples[[1]] %>% as.data.frame(),
              c1$mcmc_samples[[2]] %>% as.data.frame(),
              c1$mcmc_samples[[3]] %>% as.data.frame()) %>%
        mutate(parameter = i)

  if(i == bayes_list[1]){
    x2 <- x1
    xmcmc <- x3
  } else {
    x2 <- rbind(x2, x1)
    xmcmc <- rbind(xmcmc, x3)
  }
}

xmcmc1 <- plyr::join(xmcmc, x2 %>% dplyr::select(sig, sigSD, parameter), by = 'parameter')
ggmcmc <- xmcmc1 %>% mutate(parameter = factor(parameter, levels = bayes_list))
ggx2   <- x2 %>% mutate(var = bayes_list) %>% dplyr::select(var, everything())

median(ggmcmc[ggmcmc$parameter=='beta',]$mu_x)
median(ggmcmc[ggmcmc$parameter=='beta',]$mu_y)
median(ggmcmc[ggmcmc$parameter=='alpha',]$mu_x)
median(ggmcmc[ggmcmc$parameter=='alpha',]$mu_y)

ggplot(ggmcmc %>%
  mutate(category = ifelse(str_detect(parameter, "beta"), "beta", "alpha")) %>%
  filter(parameter %in% c('alpha', 'beta')),
  aes(parameter, mu_diff, fill = sig)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  tidybayes::stat_dist_pointinterval(aes(colour = sig)) +
  scale_colour_manual(values = c('grey', colour_group[1]), name = 'HDI Crosses Zero') +
  scale_x_discrete(labels = function(x) {
    ifelse(str_detect(x, "beta"),
           c(expression(paste(beta['ppt']^m))
             #expression(paste(beta['ppt']^sigma))
             ),
           c(expression(paste(alpha['ppt']^m))
             #expression(paste(alpha['ppt']^sigma))
             ))
  }) +
  labs(y = expression(paste(Delta, mu, ' [BPD - CON]'))) +
  facet_wrap(~category, scales = "free", nrow = 2)  +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),
        panel.grid = element_blank())+

ggplot(ggmcmc %>%
  mutate(category = ifelse(str_detect(parameter, "beta"), "beta", "alpha")) %>%
  filter(parameter %in% c('alpha_v', 'beta_v')),
  aes(parameter, mu_diff, fill = sig)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  tidybayes::stat_dist_pointinterval(aes(colour = sig)) +
  scale_colour_manual(values = c(colour_group[1]), name = 'HDI Crosses Zero') +
  scale_x_discrete(labels = function(x) {
    ifelse(str_detect(x, "beta"),
           c(#expression(paste(beta['ppt']^m)),
             expression(paste(beta['ppt']^sigma))),
           c(#expression(paste(alpha['ppt']^m)),
             expression(paste(alpha['ppt']^sigma))))
  }) +
  labs(y = expression(paste(Delta, mu, ' [BPD - CON]'))) +
  facet_wrap(~category, scales = "free", nrow = 2)  +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        strip.text = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),
        panel.grid = element_blank()) &
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=22))

ggplot(ggmcmc %>%
  mutate(category = ifelse(str_detect(parameter, "beta"), "beta", "alpha")) %>%
  filter(parameter %in% c('alpha_ref','beta_ref', 'shift_alpha_p2', 'shift_beta_p2')),
  aes(parameter, mu_diff, fill = sig)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  tidybayes::stat_dist_pointinterval(aes(colour = sig)) +
  scale_colour_manual(values = c('grey',colour_group[1]), name = 'HDI Crosses Zero') +
  scale_x_discrete(labels = function(x) {
    ifelse(str_detect(x, "beta"),
           c(expression(paste(beta['par']^ref)),
             expression(paste(Delta, beta['par']^m))),
           c(expression(paste(alpha['par']^ref)),
             expression(paste(Delta, alpha['par']^m))))
  }) +
  labs(y = expression(paste(Delta, mu, ' [BPD - CON]'))) +
  facet_wrap(~category, scales = "free", nrow = 2) +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        strip.text = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=22))

### One sided test of difference in beta values for BPD -----------

one_sided_joint <- joint_parms_ex %>%
  filter(group=='BPD_full') %>%
  mutate(beta_shift_prior = abs(beta-beta_par),
         alpha_shift_prior = abs(alpha-alpha_par))
bayes.t.test(x=one_sided_joint$beta_shift_prior, y=NULL)
bayes.t.test(x=one_sided_joint$alpha_shift_prior, y=NULL)
bayes.t.test(x=one_sided_joint$beta_par, y=NULL)

## Mean value shift ------

mean_val_shift <- int_dist %>%
  group_by(ID) %>%
  mutate(B1  = sum(beta*bP1),
         B2a = sum(beta*bP2a),
         B2b = sum(beta*bP2b),
         B3  = sum(beta*bP3),
         A1  = sum(alpha*aP1),
         A2a = sum(alpha*aP2a),
         A2b = sum(alpha*aP2b),
         A3  = sum(alpha*aP3)) %>%
  dplyr::select(ID: A3) %>%
  distinct()

p1 <- ggplot(mean_val_shift %>%
               filter(group=='BPD') %>%
               mutate(BinComp = ifelse(B1 > 0, 'Competitive', 'Prosocial')) %>%
               pivot_longer(4:5, names_to = 'Phase', values_to = 'MeanBeta'),
             aes(Phase, MeanBeta, fill = BinComp)) +
  geom_hline(yintercept = 0, colour = 'grey', linetype = 2) +
  geom_jitter(aes(colour = BinComp), alpha = 0.1) +
  stat_summary(geom = 'line', colour = 'black', aes(group = BinComp), size = 1.1) +
  stat_summary(shape = 21, colour = 'black', size = 1) +
  coord_cartesian(ylim = c(-25, 20)) +
  scale_x_discrete(labels = c('Phase 1\nPreferences', 'Phase 2\nPrior'))+
  scale_fill_brewer(palette = 'Dark2') +
  scale_colour_brewer(palette = 'Dark2')+
  labs(x = 'Phase', y = expression(paste(beta^m)))

p2 <- ggplot(mean_val_shift %>%
               filter(group == 'CON') %>%
               mutate(BinComp = ifelse(B1 > 0, 'C', 'P'),
                      BDelta = B3 - B1),
             aes(BinComp, BDelta, fill = BinComp)) +
  geom_hline(yintercept = 0, colour = 'black') +
  geom_jitter(aes(colour = BinComp), alpha = 0.5, width = 0.2) +
  stat_summary(geom = 'bar', colour = 'black', size = 0.5, alpha = 0.7) +
  coord_cartesian(ylim = c(-25, 20)) +
  labs(x = '', y = expression(paste(Delta, beta[ppt]^m, ' [B1-B3]'))) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_colour_brewer(palette = 'Dark2')

# Combine the plots side-by-side
p1 + p2 &
  theme_bw(base_size = 24) &
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        plot.margin = margin(1,1,1,1)) &
  plot_layout(widths = c(5, 1))

ggplot(mean_val_shift %>%
               filter(group=='BPD') %>%
               mutate(BinComp = ifelse(B1 > 0, 'Competitive', 'Prosocial')) %>%
               pivot_longer(4:5, names_to = 'Phase', values_to = 'MeanBeta'),
             aes(Phase, MeanBeta, group = ID, colour = MeanBeta)) +
  geom_hline(yintercept = 0, colour = 'grey', linetype = 2) +
  geom_line(colour = 'grey', alpha = 0.1)+
  geom_jitter(width = 0.1)+
  scale_colour_gradient(low = "#D95F02", high = "#1B9E77")+
  labs(x = 'Phase', y = expression(paste(beta[ppt]^m)))+
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.title = element_blank(),
        legend.position = 'none')

mean_val_shift %>%
  mutate(BinComp = ifelse(B1 > 0, 'Competitive', 'Prosocial')) %>%
  pivot_longer(B1:B2a, names_to = 'Phase', values_to = 'Mean') %>%
  filter(BinComp=='Prosocial') %>%
  dplyr::select(Phase, Mean, BinComp, ID) %>%
  distinct() %>%
  glm(Mean ~ Phase, data = .) %>%
  summary()

mean_val_shift %>%
  mutate(BinComp = ifelse(B1 > 0, 'Competitive', 'Prosocial')) %>%
  pivot_longer(B1:B2a, names_to = 'Phase', values_to = 'Mean') %>%
  filter(BinComp=='Competitive') %>%
  dplyr::select(Phase, Mean, BinComp, ID) %>%
  distinct() %>%
  glm(Mean ~ Phase, data = .) %>%
  confint()

mean_val_shift %>%
  pivot_longer(A1:A2a, names_to = 'Phase', values_to = 'Mean') %>%
  dplyr::select(Phase, Mean, ID) %>%
  distinct() %>%
  glm(Mean ~ Phase, data = .) %>%
  confint()

plot_shift_model_b <- mean_val_shift %>%
  mutate(B2a = ifelse(Model=='M1', B1, B2a)) %>%
  pivot_longer(4:7, names_to = 'Phase', values_to = 'Val') %>%
  filter(Phase%in%c('B1', 'B2a', 'B2b'))
plot_shift_model_a <- mean_val_shift %>%
  mutate(A2a = ifelse(Model=='M1', A1, A2a)) %>%
  pivot_longer(8:11, names_to = 'Phase', values_to = 'Val') %>%
  filter(Phase%in%c('A1', 'A2a', 'A2b'))

b_shift_p <- ggplot(plot_shift_model_b, aes(Phase, Val, group = group))+
  #geom_line(aes(group=ID, colour = group), alpha = 0.1)+
  stat_summary(geom = 'line', colour = 'grey', size = 1)+
  stat_summary(aes(fill = group, shape = group), colour = 'black', size = 1)+
  scale_fill_manual(values=colour_group)+
  scale_colour_manual(values=colour_group)+
  scale_shape_manual(values = c(21, 22))+
  scale_x_discrete(expand = c(0,0.05))+
  labs(y = expression(beta))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

a_shift_p <- ggplot(plot_shift_model_a, aes(Phase, Val, group = group))+
  #geom_line(aes(group=ID, colour = group), alpha = 0.1)+
  stat_summary(geom = 'line', colour = 'grey', size = 1)+
  stat_summary(aes(fill = group, shape = group), colour = 'black', size = 1)+
  scale_fill_manual(values=colour_group)+
  scale_colour_manual(values=colour_group)+
  scale_shape_manual(values = c(21, 22))+
  scale_x_discrete(expand = c(0,0.05))+
  labs(y = expression(alpha))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

a_shift_p/b_shift_p

### Belief Updates In Phase 2 ------------

b_updates <- intent_dat %>%
  filter(Phase==2) %>%
  dplyr::select(ID, trial, correctSum, RT) %>%
  mutate(group = ifelse(nchar(ID) < 10, 'BPD', 'CON'),
         kl_div_a = 0,
         kl_div_b = 0)

#sanity check for kl divs
beta_check <- data.frame(
  beta = seq(-30, 30, 0.25),
  t1 = as.numeric(BPD_full_sim$results[1,][[1]][[1]][,,1]$beta.marg1),
  t2 = as.numeric(BPD_full_sim$results[1,][[1]][[1]][,,1]$beta.cont[,38])
)

#loops for all
for(k in 1:2){
  if(k == 1){x = BPD_full_sim$results; group = 'BPD'; ID = ID1}
  if(k == 2){x = CON_full_sim$results; group = 'CON'; ID = ID2}
  for(j in 1:length(ID)){
    kl_divsa = rep(NA, 54)
    kl_divsb = kl_divsa
    for(i in 1:54){
      x[j,][[1]][[1]][,,1]$alpha.cont[,36] = x[j,][[1]][[1]][,,1]$alpha.marg1
      x[j,][[1]][[1]][,,1]$beta.cont[,36]  = x[j,][[1]][[1]][,,1]$beta.marg1

      bsa            = x[j,][[1]][[1]][,,1]$alpha.cont[,36:90]
      b_t2a          = bsa[,i+1]
      b_t1a          = bsa[,i]
      bsb            = x[j,][[1]][[1]][,,1]$beta.cont[,36:90]
      b_t2b          = bsb[,i+1]
      b_t1b          = bsb[,i]
      kl_divsa[i]    = calculate_KL_divergence(b_t2a, b_t1a)
      kl_divsb[i]    = calculate_KL_divergence(b_t2b, b_t1b)
    }
    b_updates[b_updates$group==group & b_updates$ID==ID[j],'kl_div_a'][1:54] <- kl_divsa
    b_updates[b_updates$group==group & b_updates$ID==ID[j],'kl_div_b'][1:54] <- kl_divsb
  }
}

b_updates <- b_updates %>%
  mutate(kl_div_aroll = rollapply(kl_div_a, width = 5, FUN = mean, align = "right", fill = NA),
         kl_div_broll = rollapply(kl_div_b, width = 5, FUN = mean, align = "right", fill = NA))

facet_labels <- c('alpha', 'beta')
names(facet_labels) <- c('kl_div_a', 'kl_div_b')
ggplot(b_updates %>%
         filter(trial != 1) %>%
         pivot_longer(6:7, names_to = 'Parameter', values_to = 'KL_Div'),
       aes(trial, KL_Div, fill = group, colour = group))+
  geom_smooth(alpha = 0.2)+
  scale_fill_manual(values=colour_group)+
  scale_colour_manual(values=colour_group)+
  scale_x_continuous(breaks = c(2, 25, 50))+
  facet_wrap(~Parameter, nrow=2) +
  labs(y = expression(paste('D'[KL])), x = 'Trial')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=18),
        legend.position = 'none',
        legend.title = element_blank(),
        strip.text.x = element_blank())

summary(lm(scale(kl_div_a) ~ trial*group, data = b_updates %>% filter(trial%in%c(2:54))))
summary(lm(scale(kl_div_b) ~ trial*group, data = b_updates %>% filter(trial%in%c(2:54))))

summary(lm(scale(kl) ~ trial*type,
           data = b_updates %>%
             filter(trial%in%c(2:54)) %>%
             pivot_longer(kl_div_a:kl_div_b, names_to = 'type', values_to = 'kl')
           ))

## Check Model 3 for Self Change Confirmation --------------------------------

BPD_full_sim_M3_d <- readMat('Data_Simulated/full_sim_g1_M3_only.mat')
CON_full_sim_M3_d <- readMat('Data_Simulated/full_sim_g2_M3_only.mat')

lim = 30
res = 0.25

for (i in 1:length(ID1)){
  xdist <-   data.frame(
    bP1  = as.vector(BPD_full_sim_M3_d$results[i,][[1]][[1]][,,1]$beta.marg1),
    bP2a = as.vector(BPD_full_sim_M3_d$results[i,][[1]][[1]][,,1]$beta.marg2a),
    bP2b = as.vector(BPD_full_sim_M3_d$results[i,][[1]][[1]][,,1]$beta.marg2b),
    bP3  = as.vector(BPD_full_sim_M3_d$results[i,][[1]][[1]][,,1]$beta.marg3),
    aP1  = as.vector(BPD_full_sim_M3_d$results[i,][[1]][[1]][,,1]$alpha.marg1),
    aP2a = as.vector(BPD_full_sim_M3_d$results[i,][[1]][[1]][,,1]$alpha.marg2a),
    aP2b = as.vector(BPD_full_sim_M3_d$results[i,][[1]][[1]][,,1]$alpha.marg2b),
    aP3  = as.vector(BPD_full_sim_M3_d$results[i,][[1]][[1]][,,1]$alpha.marg3),
    beta = seq(-lim, lim, res),
    alpha = seq(0, lim, res/2),
    ID = ID1[i],
    group = 'BPD',
    Model = 'M3'
  )
  if(i>1){BPDdist_M3 <- rbind(BPDdist_M3, xdist)} else {BPDdist_M3 <- xdist}
}
for (i in 1:length(ID2)){
  zdist <-   data.frame(
    bP1  = as.vector(CON_full_sim_M3_d$results[i,][[1]][[1]][,,1]$beta.marg1),
    bP2a = as.vector(CON_full_sim_M3_d$results[i,][[1]][[1]][,,1]$beta.marg2a),
    bP2b = as.vector(CON_full_sim_M3_d$results[i,][[1]][[1]][,,1]$beta.marg2b),
    bP3  = as.vector(CON_full_sim_M3_d$results[i,][[1]][[1]][,,1]$beta.marg3),
    aP1  = as.vector(CON_full_sim_M3_d$results[i,][[1]][[1]][,,1]$alpha.marg1),
    aP2a = as.vector(CON_full_sim_M3_d$results[i,][[1]][[1]][,,1]$alpha.marg2a),
    aP2b = as.vector(CON_full_sim_M3_d$results[i,][[1]][[1]][,,1]$alpha.marg2b),
    aP3  = as.vector(CON_full_sim_M3_d$results[i,][[1]][[1]][,,1]$alpha.marg3),
  beta  = seq(-lim, lim, res),
  alpha = seq(0, lim, res/2),
  ID = ID2[i],
  group = 'CON',
  Model = 'M3'
  )
  if(i>1){CONdist_M3 <- rbind(CONdist_M3, zdist)} else {CONdist_M3 <- zdist}
}

int_dist_M3 <- rbind(BPDdist_M3, CONdist_M3)

mean_val_shift_M3 <- int_dist_M3 %>%
  group_by(ID) %>%
  mutate(B1  = sum(beta*bP1),
         B2a = sum(beta*bP2a),
         B2b = sum(beta*bP2b),
         B3  = sum(beta*bP3),
         A1  = sum(alpha*aP1),
         A2a = sum(alpha*aP2a),
         A2b = sum(alpha*aP2b),
         A3  = sum(alpha*aP3),
         B1_sd = sqrt(sum((beta-B1)^2 * bP1)),
         B2a_sd= sqrt(sum((beta-B2a)^2 * bP2a)),
         B3_sd = sqrt(sum((beta-B3)^2 * bP3)),
         A1_sd = sqrt(sum((alpha-A1)^2* aP1)),
         A2a_sd= sqrt(sum((alpha-A2a)^2 *aP2a)),
         A3_sd = sqrt(sum((alpha-A3)^2* aP3)),
         Delta_B    = B3-B1,
         Delta_A    = A3-A1,
         Delta_B_sd = B3_sd - B1_sd,
         Delta_A_sd = A3_sd - A1_sd) %>%
  dplyr::select(ID: Delta_A_sd) %>%
  pivot_longer(Delta_B:Delta_A_sd, names_to = 'Delta', values_to = 'Shift')  %>%
  distinct()

ggplot(mean_val_shift_M3 %>% filter(Delta%in%c('Delta_B', 'Delta_A')),
       aes(Delta, Shift, fill = group))+
  geom_point(shape=21, alpha=0.3, size=3,position = position_dodge(width=0.5))+
  geom_boxplot(width=0.5, outliers = F, alpha = 0.7)+
  labs(x='', y=expression(paste(Delta, theta['ppt']^m)))+
  scale_x_discrete(labels = c(expression(paste(alpha['ppt']^m)), expression(paste(beta['ppt']^m))))+
  scale_fill_manual(values=colour_group)+
  stat_compare_means(label = 'p.signif', size = 8, label.y = 15)+
  theme_bw(base_size=22)+
ggplot(mean_val_shift_M3 %>% filter(Delta%in%c('Delta_B_sd', 'Delta_A_sd')),
       aes(Delta, Shift, fill = group))+
  geom_point(shape=21, alpha=0.3, size=3,position = position_dodge(width=0.5))+
  geom_boxplot(width=0.5, outliers = F, alpha = 0.7)+
  labs(x='', y=expression(paste(Delta, theta['ppt']^sigma)))+
  scale_x_discrete(labels = c(expression(paste(alpha['ppt']^sigma)), expression(paste(beta['ppt']^sigma))))+
  scale_y_continuous(labels = c(0, -2, -4, -6), breaks = c(0, -2, -4, -6), limits = c(-6, 0))+
  scale_fill_manual(values=colour_group)+
  stat_compare_means(label = 'p.signif', size = 8, label.y = -6)+
  theme_bw(base_size=22)&
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=22),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),
        axis.title.y = element_text(vjust=-1))

summary(lm(abs(Shift) ~ group + Delta, mean_val_shift_M3%>% filter(Delta%in%c('Delta_B', 'Delta_A'))))
confint(lm(abs(Shift) ~ group + Delta, mean_val_shift_M3%>% filter(Delta%in%c('Delta_B', 'Delta_A'))))

summary(lm(abs(Shift) ~ group + Delta, mean_val_shift_M3%>% filter(Delta%in%c('Delta_B_sd', 'Delta_A_sd'))))
confint(lm(abs(Shift) ~ group + Delta, mean_val_shift_M3%>% filter(Delta%in%c('Delta_B_sd', 'Delta_A_sd'))))

summary(lm(abs(Shift) ~ group * Delta, mean_val_shift_M3%>% filter(Delta%in%c('Delta_B', 'Delta_A'))))
confint(lm(abs(Shift) ~ group * Delta, mean_val_shift_M3%>% filter(Delta%in%c('Delta_B', 'Delta_A'))))
summary(lm(abs(Shift) ~ group * Delta, mean_val_shift_M3%>% filter(Delta%in%c('Delta_B_sd', 'Delta_A_sd'))))
confint(lm(abs(Shift) ~ group * Delta, mean_val_shift_M3%>% filter(Delta%in%c('Delta_B_sd', 'Delta_A_sd'))))

### Against ideal outcomes --------------------------------------------------

m_shift <- ggplot(shift_beta_l %>%
                  filter(parameter == 'm_shift'),
              aes(s_u, o_u, fill = shift))+
  scale_fill_gradient(high = 'black', low = 'white',
                      name = expression(paste(Delta, beta[ppt]^mu)),
                      breaks = c(5, 20, 40))+
  geom_tile() +
  geom_point(data = joint_parms_ex %>%
                 filter(group=='CON_full',
                        beta_v < 15, beta_ref < 15),
             aes(x=beta_v, y = beta_ref),
             shape=21,fill=colour_group[2],size=4, alpha = 0.4)+
  geom_smooth(data = joint_parms_ex %>%
                 filter(group=='CON_full',
                        beta_v < 15, beta_ref < 15),
             aes(x=beta_v, y = beta_ref),
             fill = 'grey', formula = y ~ x + I(x^2),
             method = "lm", se = FALSE, colour = "black") &
  labs(y = expression(paste('Other Unc.')),
       x = expression(paste('Self Unc.'))
       )
sd_shift <- ggplot(shift_beta_l %>%
         filter(parameter == 'sd_shift'),
       aes(s_u, o_u, fill = shift))+
  scale_fill_gradient(low = 'black', high = 'white',
                      name = expression(paste(Delta, beta[ppt]^sigma)),
                      breaks = c(0, -5, -10))+
  geom_tile() +
  geom_point(data = joint_parms_ex %>%
                 filter(group=='CON_full',
                        beta_v < 15, beta_ref < 15),
             aes(x=beta_v, y = beta_ref),
             shape=21,fill=colour_group[2],size=4, alpha = 0.4)+
  geom_smooth(data = joint_parms_ex %>%
                 filter(group=='CON_full',
                        beta_v < 15, beta_ref < 15),
             aes(x=beta_v, y = beta_ref),
             fill = 'grey', formula = y ~ x + I(x^2),
             method = "lm", se = FALSE, colour = "black")

ideal_plot <- (sd_shift|m_shift)&
  labs(y = expression(paste('Other Unc.')),
       x = expression(paste('Self Unc.'))
       ) &
  scale_y_discrete(breaks = c(0, 5, 10, 15), expand = c(0,0)) &
  scale_x_discrete(breaks = c(0, 5, 10, 15), expand = c(0,0)) &
  theme_bw(base_size=18)&
  theme(legend.position = 'top')

ideal_plot

# Check against model predictions

delta_shift_real <- joint_parms_ex %>%
  filter(group == 'CON_full',
         beta_v < 15,
         beta_ref < 15) %>%
  mutate(delta_beta_ppt_m = abs(beta - beta_hat),
         delta_beta_ppt_sd= -abs(beta_v - beta_hat_sd),
         beta_v = round(beta_v),
         beta_ref = round(beta_ref)) %>%
  dplyr::select(delta_beta_ppt_m, delta_beta_ppt_sd, beta_v, beta_ref)

delta_shift_model <- shift_beta_l %>%
  pivot_wider(id_cols = c(s_u, o_u), names_from = parameter, values_from = shift) %>%
  rename(beta_v = s_u,
         beta_ref = o_u,
         delta_beta_ppt_m = m_shift,
         delta_beta_ppt_sd = sd_shift) %>%
  mutate(beta_v = as.numeric(beta_v),
         beta_ref = as.numeric(beta_ref))

delta_shift <- left_join(delta_shift_real, delta_shift_model,
                         by = c("beta_v", "beta_ref"),
                         suffix = c("_real", "_simulated"))

comparison_plot <- ggplot(delta_shift,
       aes(delta_beta_ppt_sd_real, delta_beta_ppt_sd_simulated)) +
  geom_point(shape=21, size = 3, fill = 'grey')+
  geom_smooth(method='lm', colour = 'black')+
  scale_x_continuous(breaks = c(0, -5, -10))+
  scale_y_continuous(breaks = c(0, -5, -10))+
  scale_alpha_continuous(range = c(0.2,1))+
  theme_bw(base_size=18)+
ggplot(delta_shift,
       aes(delta_beta_ppt_m_real, delta_beta_ppt_m_simulated)) +
  geom_point(shape=21, size= 3, fill = 'grey')+
  geom_smooth(method='lm', colour = 'black')+
  scale_x_continuous(breaks = c(0, 5, 10))+
  scale_y_continuous(breaks = c(0, 25, 50))+
  theme_bw(base_size=18)&
  labs(x=expression(paste('Real ',Delta, beta[ppt])),
       y=expression(paste('Sim. ',Delta, beta[ppt])))&
  theme(legend.position = 'none')

ideal_plot/comparison_plot

# Correlations ------------------------------------------------------------

## Atttributions ---------

# Transform the data
transformed_data <- joint_parms_ex %>%
  ungroup() %>%
  pivot_longer(cols = c(distancealpha, distancebeta, alpha_ref, beta_ref), names_to = 'Parm', values_to = 'Val1') %>%
  pivot_longer(cols = c(HI, SI), names_to = 'Att', values_to = 'Val2') %>%
  dplyr::select(Parm, Att, Val1, Val2, group) %>%
  distinct()

# Initialize a data frame to store the results
results <- data.frame(
  Parm = character(),
  Att = character(),
  Rho = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  sig=numeric(),
  p=numeric(),
  stringsAsFactors = FALSE
)

# Loop through each combination of Parm and Att to perform Spearman correlations
for (parm in unique(transformed_data$Parm)) {
  for (att in unique(transformed_data$Att)) {
    # Filter data for the current Parm and Att
    data_subset <- transformed_data %>%
      filter(Parm == parm, Att == att)

    x_val1 <- lm(Val1~group, data_subset)
    x_val2 <- lm(Val2~group, data_subset)

    # Perform the Spearman correlation
    corr_result <- cor.test(x_val1$residuals, x_val2$residuals, type = "spearman")

    # Append the results to the results data frame
    results <- rbind(results, data.frame(Parm = parm,
                                         Att = att,
                                         Rho = as.numeric(corr_result$estimate),
                                         CI_Lower = corr_result$conf.int[1],
                                         CI_Upper = corr_result$conf.int[2],
                                         sig = ifelse((corr_result$conf.int[1] < 0 & corr_result$conf.int[2] < 0) |
                                                        (corr_result$conf.int[1] > 0 & corr_result$conf.int[2] > 0),
                                                      'Yes', 'No'),
                                         p=corr_result$p.value))
  }
}

p.adjust(results$p[1:4], method = 'fdr')
p.adjust(results$p[5:8], method = 'fdr')

# Create the ggplot2 object
ggplot(results, aes(x = Att, y = Rho,  colour = sig)) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper, colour = sig), width = 0.1) +
  geom_point(size=3) +
  geom_hline(yintercept = 0)+
  labs(y = expression(paste(rho, ' ± 95% CI'))) +
  facet_wrap(~Parm, scales = 'free_x', nrow = 1)+
  scale_x_discrete(labels = c(
    'HI', 'SI'
    ))+
  scale_colour_manual(values = c('grey','#BA2D0B'),
                      name = '95%CI\nCrosses\nZero')+
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust=-1),
        strip.text.x = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank())

residual_HI <-  lm(HI ~ shift_beta_p2, data = joint_parms_ex)
cor.test(joint_parms_ex$beta_ref, residual_HI$residuals, method = 'spearman')
cor.test(joint_parms_ex$distancebeta, residual_HI$residuals, method = 'spearman')

## Psychometrics ---------

supp_8_test <- plyr::join(
  mean_val_shift_M3,
  intent_dat %>% dplyr::select(-trial:-correct, -RT, -Phase, -correctSum) %>% distinct(),
  by = c('ID')
) %>%
  filter(Delta %in% c('Delta_B', 'Delta_A')) %>%
  dplyr::select(ID, group, Model, Delta, Shift, MZQ_TotalScore, CTQtot, RGPTSA, RGPTSB, B1:A3_sd) %>%
  pivot_wider(id_cols = c(ID, group, MZQ_TotalScore:RGPTSB, B1:A3_sd), names_from = Delta, values_from = Shift) %>%
  ungroup() %>%
  mutate(
    Delta_A = abs(Delta_A),
    Delta_B = abs(Delta_B)
  ) %>%
  distinct()

network_pcor <- estimateNetwork(
                         as.matrix(supp_8_test %>%
                            ungroup() %>%
                            dplyr::select(MZQ_TotalScore, CTQtot, RGPTSB, Delta_A, Delta_B, -ID)
                          ),
                default = 'pcor', checkNumeric = T)

network_pcor$results
plot(network_pcor)
boot_pcor <- bootnet(network_pcor, nBoots = 5000)
summary(boot_pcor) %>% filter((CIlower < 0 & CIupper < 0) | (CIlower>0&CIupper>0))

### Permutation test ----
perms <- 5000
perm_graph <- list()

for(i in 1:perms){
network_pcor_perm <- estimateNetwork(
                          as.matrix(supp_8_test %>%
                            ungroup() %>%
                            dplyr::select(MZQ_TotalScore, CTQtot, RGPTSB, Delta_A, Delta_B, -ID)
                            ) %>%
                            apply(., 2, sample),
                  default = 'pcor', checkNumeric = T)

perm_graph[[i]] <- network_pcor_perm$graph %>%
  as.data.frame() %>%
  mutate(cross_c = as.vector(rownames(network_pcor_perm$graph)),
         perm = i)

}

# Initialize a matrix to store p-values
p_value_matrix <- matrix(NA, nrow = 5, ncol = 5)
plot_data <- list()

# Calculate p-values for each cell
for (i in 1:5) {
  for (j in 1:5) {
    # Extract the distribution of the cell (i, j) from permuted matrices
    permuted_values <- sapply(perm_graph, function(mat) mat[i, j])

    # Two-tailed p-value: proportion of permuted values more extreme than the true value
    true_value <- network_pcor$graph[i, j]
    if(true_value!=0){
    p_value_matrix[i, j] <- mean(abs(permuted_values - mean(permuted_values)) >= abs(true_value - mean(permuted_values)))
    plot_data[[paste(rownames(network_pcor$graph)[i],colnames(network_pcor$graph)[j], sep = "_")]] <- data.frame(
      Cell = paste(rownames(network_pcor$graph)[i],colnames(network_pcor$graph)[j], sep = ""),
      Permuted = permuted_values,
      TrueValue = true_value
      )
    }
  }
}

colnames(p_value_matrix) <- colnames(network_pcor$graph)
rownames(p_value_matrix) <- rownames(network_pcor$graph)
p_value_matrix

# Combine all cell data into a single data frame
plot_data_df <- bind_rows(plot_data, .id = "Cell")

# Plot the distributions with ggplot2
ggplot(plot_data_df, aes(x = Permuted)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = TrueValue), color = "red", linetype = "dashed") +
  facet_wrap(~ Cell, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Permuted Distributions with True Values",
    x = "Value",
    y = "Frequency"
  )

### Plot surviving edges -----

#plot actual outcomes:
network_pcor_a <- network_pcor

for (i in 1:5) {
  for (j in 1:5) {
    if(is.na(p_value_matrix[i,j]) | p_value_matrix[i,j] > 0.05){
    network_pcor_a$graph[i,j] <- 0
    }
  }
}

library(qgraph)
qgraph(
  network_pcor_a$graph,
  layout = 'spring',
  labels = c('MZQ', 'CTQ', 'RGPTS',
             expression(paste(Delta, alpha)),
             expression(paste(Delta,beta))
             ),
  label.cex = 1,
  vsize = 10,
  color = c(rep('black', 3), rep('white', 2)),
  posCol = 'blue',
  edge.labels = T,
  edge.label.bg = T,
  edge.label.margin = 0.01,
  edge.label.cex = 2,
  repulsion = 1.1
               )

### Control for group status ------------------------------------------------

supp_8_group <- plyr::join(
  mean_val_shift_M3,
  intent_dat %>% dplyr::select(-trial:-correct, -RT, -Phase, -correctSum) %>% distinct(),
  by = c('ID')
) %>%
  filter(Delta %in% c('Delta_B', 'Delta_A')) %>%
  dplyr::select(ID, group, Model, Delta, Shift, MZQ_TotalScore, CTQtot, RGPTSA, RGPTSB, B1:A3_sd) %>%
  pivot_wider(id_cols = c(ID, group, MZQ_TotalScore:RGPTSB, B1:A3_sd), names_from = Delta, values_from = Shift) %>%
  ungroup() %>%
  mutate(
    Delta_A = abs(Delta_A),
    Delta_B = abs(Delta_B)
  ) %>%
  distinct() %>%

  ## Regress each specified variable against group and take residuals
  mutate(
    MZQ_TotalScore   = resid(lm(MZQ_TotalScore ~ group, na.action = na.exclude)),
    CTQtot = resid(lm(CTQtot ~ group, na.action = na.exclude)),
    RGPTSB = resid(lm(RGPTSB ~ group, na.action = na.exclude))
  )

network_pcor_group <- estimateNetwork(
                         as.matrix(supp_8_group %>%
                            ungroup() %>%
                            dplyr::select(MZQ_TotalScore, CTQtot, RGPTSB, Delta_A, Delta_B, -ID)
                          ),
                default = 'pcor', checkNumeric = T)

network_pcor_group$results
plot(network_pcor_group)
boot_pcor_group <- bootnet(network_pcor_group, nBoots = 5000)
summary(boot_pcor_group) %>% filter((CIlower < 0 & CIupper < 0) | (CIlower>0&CIupper>0))

# Supplementary Figure 1 --------------------------------------------------

bpd_group <- data.frame(
alpha   = c((1/(1+exp(-BPD_full$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]][,1])))*30,
            (1/(1+exp(-BPD_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,5][[1]][[1]][,1])))),
beta    = c(BPD_full$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]][,2],
            BPD_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,5][[1]][[1]][,2]),
alpha_v = c(exp(BPD_full$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]][,3]),
            exp(BPD_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,5][[1]][[1]][,3])),
beta_v  = c(exp(BPD_full$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]][,4]),
            exp(BPD_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,5][[1]][[1]][,4])),
alpha_p = c((1/(1+exp(-BPD_full$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]][,5])))*30,
            (1/(1+exp(-BPD_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,5][[1]][[1]][,5])))),
beta_p  = c(BPD_full$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]][,6],
            BPD_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,5][[1]][[1]][,6]),
alpha_r = c(exp(BPD_full$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]][,7]),
            exp(BPD_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,5][[1]][[1]][,7])),
beta_r  = c(exp(BPD_full$cbm[,,1]$output[,,1]$group.mean[,5][[1]][[1]][,8]),
            exp(BPD_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,5][[1]][[1]][,8])),
group   = 'BPD',
type    = c('mean', 'error')
)

con_group <- data.frame(
alpha   = c((1/(1+exp(-CON_full$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,1])))*30,
            (1/(1+exp(-CON_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,1])))),
beta    = c(CON_full$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,2],
            CON_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,2]),
alpha_v = c(exp(CON_full$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,3]),
            exp(CON_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,3])),
beta_v  = c(exp(CON_full$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,4]),
            exp(CON_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,4])),
alpha_p = NA,
beta_p  = NA,
alpha_r = c(exp(CON_full$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,5]),
            exp(CON_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,5])),
beta_r  = c(exp(CON_full$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,6]),
            exp(CON_full$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,6])),
group   = 'CON',
type    = c('mean', 'error')
)

parm_group <- rbind(bpd_group, con_group) %>%
  pivot_longer(1:8, names_to = 'Parm', values_to = 'Val') %>%
  mutate(category = ifelse(str_detect(Parm, "beta"), "beta", "alpha"))

parm_group$Parm <- factor(parm_group$Parm,
                          levels = c("alpha", "beta", "alpha_v", "beta_v",
                                     "alpha_p", "beta_p", "alpha_r", "beta_r"))

# Separate the data into mean and error datasets
means  <- parm_group %>% filter(type == "mean")
errors <- parm_group %>% filter(type == "error")

# Merge the datasets on group, Parm, and category
parm_group_full <- means %>%
  left_join(errors, by = c("group", "Parm", "category"), suffix = c("_mean", "_error"))

ggplot(parm_group_full,
       aes(Parm, Val_mean, fill = group)) +
  geom_col(position = 'dodge')+
  geom_errorbar(aes(ymin = Val_mean - Val_error, ymax = Val_mean + Val_error),
                width = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=colour_group)+
  facet_wrap(~category, scale = 'free')+
  scale_y_continuous(expand=c(0,0.2))+
  scale_x_discrete(labels = function(x) {
    ifelse(str_detect(x, "beta"),
           c(expression(paste(beta['ppt']^m)),
             expression(paste(beta['ppt']^sigma)),
             expression(paste(beta['par']^m)),
             expression(paste(beta['par']^sigma))
             ),
           c(expression(paste(alpha['ppt']^m)),
             expression(paste(alpha['ppt']^sigma)),
             expression(paste(alpha['par']^m)),
             expression(paste(alpha['par']^sigma))
             ))
  },
  expand = c(0,0.5)) +
  labs(y = 'Value\n(Transformed Space)')+
  theme_bw(base_size=18)+
  theme(legend.position = c(0.9, 0.23),
        legend.background = element_rect(colour = 'black'),
        strip.text.x = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank())

# Supplementary Figure 2 --------------------------------------------------

plot_parms <- rbind(BPD_parms,
                    CON_parms) %>%
         pivot_longer(2:9, names_to = 'parms', values_to = 'est')

ggplot(plot_parms,
       aes(est, fill = group))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=colour_group)+
  facet_wrap(~parms, scale = 'free', nrow = 4,
             labeller = labeller(
               parms = c(
                      alpha=expression(paste(alpha['ppt']^m)),
                      beta=expression(paste(beta['ppt']^m)),
                      alpha_v=expression(paste(alpha['ppt']^sigma)),
                      beta_v=expression(paste(beta['ppt']^sigma)),
                      alpha_ref=expression(paste(alpha['ref']^sigma)),
                      beta_ref=expression(paste(beta['ref']^sigma))
                      )
               )
             )+
  labs(y = 'Density', x = 'Value\n(Transformed Space)')+
  theme_bw()+
  theme(legend.position = 'none',
        text = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(t=1,r=1,b=1,l=1),
        panel.grid = element_blank())

# Supplementary Figure 3 --------------------------------------------------

range_w <- seq(0, 1, 0.1)
int_test <- list()

for(i in 1:length(range_w)){

priors_int <- c(5, -15, 4, 4, 5, 15, 4, 4, range_w[i])
int_test_1 <- ABA_shift_Gen_integratedprior(priors_int,
                                           intent_dat_126[intent_dat_126$ID=='PDC195',],
                                           sim=1)
int_test[[i]] <- int_test_1$marginals %>%
  as.data.frame() %>%
  mutate(w = range_w[i])

}

int_plot <- bind_rows(int_test)

ggplot(int_plot, aes(beta_grid, Phase2a_beta, group = w))+
  geom_line(alpha=0.5, aes(colour = w))+
  geom_line(data=int_plot %>% filter(w==0),   aes(beta_grid, Phase2a_beta), size = 1, colour = '#F6AA1C')+
  #geom_line(data=int_plot %>% filter(w==0.5), aes(beta_grid, Phase2a_beta), size = 1, colour = '#941B0C')+
  geom_line(data=int_plot %>% filter(w==1),   aes(beta_grid, Phase2a_beta), size = 1, colour = '#220901')+
  scale_colour_gradient2(mid='#941B0C', high='#220901', low='#F6AA1C', midpoint = 0.5, guide = 'colourbar',
                         name = expression(omega))+
  labs(x=expression(beta), y = expression(paste('p(', beta,')')))+
  theme_bw(base_size=20)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y =  element_blank(),
        legend.position = 'top',
        legend.key.width = unit(1.5,'cm'),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_blank(),
        plot.margin = margin(0,0,0,0),
        panel.border = element_blank())

# Supplementary Figure 4 --------------------------------------------------

ggplot(int_dist %>% filter(ID %in% joint_parms_ex[joint_parms_ex$beta>1,]$ID[6], group == 'BPD'))+
  geom_line(aes(beta, bP1, colour = '1 - P1', group =ID), size = 1.2)+
  geom_line(aes(beta, bP2a, colour = '2A - Prior', group =ID), size = 1.2)+
  geom_line(aes(beta, bP2b, colour = '2B - Post.', group =ID), size = 1.2)+
  geom_vline(xintercept = joint_parms_ex[joint_parms_ex$beta>1,]$server_beta_par[6])+
  labs(title = '', x = expression(paste(beta)), y = expression(paste('p(',beta,')')))+
  scale_colour_manual(values = c('black','#FBD19D', '#F7941D', 'grey'), name='Type')+
  theme(axis.title.y = element_blank(),
        legend.position = c(0.25, 0.75))&
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16))

# Supplementary Figure 5 --------------------------------------------------

t.test(sumCbpd, sumCcon)
predCor <- rbind(predCorbpd %>% mutate(group='BPD'),
                 predCorcon %>% mutate(group='CON'))

psych_curves %>%
  group_by(trial, group) %>%
  mutate(corSumT = sum(correct)/54)%>%
  ggplot(aes(trial, corSumT, colour = group, fill = group))+
  geom_hline(yintercept = 0.5)+
  stat_summary(geom='ribbon')+
  stat_summary(geom='line')+
  scale_y_continuous(limits=c(0.3, 1))+
  scale_color_manual(values = colour_group)+
  scale_fill_manual(values = colour_group)+
  labs(x = 'Trial', y = 'p(correct)', title = 'Real Data')+
  theme_bw(base_size=18)+
  theme(legend.position = c(0.7, 0.15),
        legend.title = element_blank(),
        panel.grid = element_blank())+
predCor%>%
  group_by(V2, group) %>%
  mutate(corSumT = sum(cor)/54)%>%
  ggplot(aes(V2, corSumT, colour = group, fill = group))+
  geom_hline(yintercept = 0.5)+
  stat_summary(geom='ribbon')+
  stat_summary(geom='line')+
  scale_y_continuous(limits=c(0.3, 1))+
  scale_color_manual(values = colour_group)+
  scale_fill_manual(values = colour_group)+
  labs(x = 'Trial', y = 'p(correct)', title = 'Simulated')+
  theme_bw(base_size=18)+
  theme(legend.position = 'none',
        panel.grid = element_blank())


# Supplementary Figure 6 --------------------------------------------------

parm_vals <- intent_dat_126 %>%
  dplyr::select(ID, server_alpha_ppt, server_alpha_par) %>%
  distinct() %>%
  pivot_longer(c(server_alpha_ppt, server_alpha_par), names_to = 'player', values_to = 'alpha_val') %>%
  mutate(player = ifelse(player == 'server_alpha_ppt', 'PPT', 'PAR')) %>%
  plyr::join(.,
        intent_dat_126 %>%
  dplyr::select(ID, server_beta_ppt, server_beta_par) %>%
  distinct() %>%
  pivot_longer(c(server_beta_ppt, server_beta_par), names_to = 'player', values_to = 'beta_val') %>%
  mutate(player = ifelse(player == 'server_beta_ppt', 'PPT', 'PAR'),
         group=ifelse(nchar(ID)<7, 'BPD', 'CON')) %>%
  dplyr::select(ID, beta_val, player, group),
  by = c('ID', 'player')
  )

density_plot <- ggplot(parm_vals %>% filter(player == 'PPT'), aes(alpha_val, beta_val))+
  geom_density_2d(contour_var = 'ndensity',
                  colour = "black",
                  bins = 10)+
  geom_jitter(aes(fill = group), shape=21,size=2)+
  labs(x = expression(alpha), y = expression(beta))+
  scale_x_continuous(expand=c(0,1), limits = c(0, 30))+
  scale_y_continuous(expand=c(0,1), limits = c(-20, 20))+
  scale_fill_manual(values=colour_group)+
  theme_bw() +
ggplot(parm_vals %>% filter(player == 'PAR'), aes(alpha_val, beta_val))+
  geom_density_2d(contour_var = 'ndensity',
                  colour = "black",
                  bins = 10)+
  geom_jitter(aes(fill = group), shape=21,size=2)+
  labs(x = expression(alpha), y = expression(beta))+
  scale_x_continuous(expand=c(0,1), limits = c(0, 30))+
  scale_y_continuous(expand=c(0,1), limits = c(-20, 20))+
  scale_fill_manual(values=colour_group)+
  theme_bw() &
theme(legend.position = c(0.8, 0.85),
      legend.title = element_blank(),
        text = element_text(size=18),
        panel.grid = element_blank())
density_plot

# Supplementary Figure 7 --------------------------------------------------

real_cor_parms_metric  <- rcorr(as.matrix(joint_parms_ex[,
                                                   c('MZQ_TotalScore' ,
                                                     'CAMSQ_Self', 'CAMSQ_Other',
                                                     'RGPTSA', 'RGPTSB', 'CTQtot',
                                                     #'A_ETS_Mistrust', 'A_ETS_Trust', 'A_ETS_Credulity',
                                                     'alpha', 'beta', 'alpha_v', 'beta_v', 'alpha_ref', 'beta_ref', 'shift_alpha_p2', 'shift_beta_p2'
                                                     )]
                                   ), type = 'spearman')

ggcorrplot::ggcorrplot(real_cor_parms_metric$r[1:13, 1:13],
                       p.mat = na.fill(real_cor_parms_metric$P[1:13, 1:13], 0),
                       lab = T,
                       colors = c( '#291720','white','#007A5A'),
                       legend.title = 'Pearson\nR',
                       sig.level = 0.05,
                       type = 'upper',
                       show.diag = F,
                       insig='blank',
                       pch.cex = 15,
                       pch.col = 1) +
   scale_x_discrete(labels = c(
     'MZQ',
     'CAMSQ (Self)',
     'CAMSQ (Other)',
     'RGPTSA',
     'RGPTSB',
     'CTQ',
     'Mistrust',
     'Trust',
     'Credulity',
     expression(paste(alpha['ppt']^m)),
     expression(paste(beta['ppt']^m)),
     expression(paste(alpha['ppt']^sigma))
   ),
   expand=c(0,0)) +
     scale_y_discrete(labels = c(
     'CAMSQ (Self)',
     'CAMSQ (Other)',
     'RGPTSA',
     'RGPTSB',
     'CTQ',
     'Mistrust',
     'Trust',
     'Credulity',
     expression(paste(alpha['ppt']^m)),
     expression(paste(beta['ppt']^m)),
     expression(paste(alpha['ppt']^sigma)),
     expression(paste(beta['ppt']^sigma))
   ),
    expand=c(0,0)) +
  theme_bw(base_size=18)+
  theme(axis.title = element_blank(),
        axis.text  = element_text(size=18),
        legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45))

ggcorrplot::ggcorrplot(real_cor_parms_metric$r[c(1:9, 14:17), c(1:9, 14:17)],
                       p.mat = na.fill(real_cor_parms_metric$P[c(1:9, 14:17), c(1:9, 14:17)], 0),
                       lab = T,
                       colors = c( '#291720','white','#007A5A'),
                       legend.title = 'Pearson\nR',
                       sig.level = 0.05,
                       type = 'upper',
                       show.diag = F,
                       insig='blank',
                       pch.cex = 15,
                       pch.col = 1) +
   scale_x_discrete(labels = c(
     'MZQ',
     'CAMSQ (Self)',
     'CAMSQ (Other)',
     'RGPTSA',
     'RGPTSB',
     'CTQ',
     'Mistrust',
     'Trust',
     'Credulity',
     expression(paste(alpha['par']^ref)),
     expression(paste(beta['par']^ref)),
     expression(paste(Delta, alpha['par']^m))
   ),
   expand=c(0,0)) +
     scale_y_discrete(labels = c(
     'CAMSQ (Self)',
     'CAMSQ (Other)',
     'RGPTSA',
     'RGPTSB',
     'CTQ',
     'Mistrust',
     'Trust',
     'Credulity',
     expression(paste(alpha['par']^ref)),
     expression(paste(beta['par']^ref)),
     expression(paste(Delta, alpha['par']^m)),
     expression(paste(Delta, beta['par']^m))
   ),
    expand=c(0,0)) +
  theme_bw(base_size=18)+
  theme(axis.title = element_blank(),
        axis.text  = element_text(size=18),
        legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45))


# Supplementary Figure 8 --------------------------------------------------

# Transform the data
transformed_data2 <- joint_parms_ex %>%
  ungroup() %>%
  mutate(
    CTQtot         = resid(lm(CTQtot ~ group, na.action = na.exclude)),
    RGPTSB         = resid(lm(RGPTSB ~ group, na.action = na.exclude)),
    MZQ_TotalScore = resid(lm(MZQ_TotalScore ~ group, na.action = na.exclude)),
    CAMSQ_Self     = resid(lm(CAMSQ_Self ~ group, na.action = na.exclude)),
    CAMSQ_Other    = resid(lm(CAMSQ_Other ~ group, na.action = na.exclude))
  ) %>%
  pivot_longer(cols = c(shift_alpha_p2, shift_beta_p2, alpha_ref, beta_ref), names_to = 'Parm', values_to = 'Val1') %>%
  pivot_longer(cols = c(CTQtot, RGPTSB, MZQ_TotalScore, CAMSQ_Self, CAMSQ_Other), names_to = 'Att', values_to = 'Val2') %>%
  dplyr::select(Parm, Att, Val1, Val2, group) %>%
  distinct()

# Initialize a data frame to store the results
results2 <- data.frame(
  Parm = character(),
  Att = character(),
  Rho = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  sig=numeric(),
  p=numeric(),
  stringsAsFactors = FALSE
)

# Loop through each combination of Parm and Att to perform Spearman correlations
for (parm in unique(transformed_data2$Parm)) {
  for (att in unique(transformed_data2$Att)) {
    # Filter data for the current Parm and Att
    data_subset2 <- transformed_data2 %>%
      filter(Parm == parm, Att == att)

    y                 <- lm(data_subset2$Val1~joint_parms_ex$distancealpha+joint_parms_ex$distancebeta)
    data_subset2$Val1 <- y$residuals

    # Perform the Spearman correlation
    corr_result2 <- cor.test(data_subset2$Val1, data_subset2$Val2, type = "spearman")

    # Append the results to the results data frame
    results2 <- rbind(results2, data.frame(Parm = parm,
                                            Att = att,
                                            Rho = as.numeric(corr_result2$estimate),
                                            CI_Lower = corr_result2$conf.int[1],
                                            CI_Upper = corr_result2$conf.int[2],
                                            sig = ifelse((corr_result2$conf.int[1] < 0 & corr_result2$conf.int[2] < 0) |
                                                         (corr_result2$conf.int[1] > 0 & corr_result2$conf.int[2] > 0),
                                                       'Yes', 'No'),
                                            p=corr_result2$p.value))
  }
}

# Create the ggplot2 object
ggplot(results2 %>%
         filter(Att %in% c('CTQtot', 'RGPTSB')),
       aes(x = Att, y = Rho,  colour = sig)) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper, colour = sig), width = 0.1) +
  geom_point(size=3) +
  geom_hline(yintercept = 0)+
  labs(y = expression(paste(rho, ' ± 95% CI'))) +
  facet_wrap(~Parm, scales = 'free_x', nrow=1)+
  scale_x_discrete(labels = c(
   'CTQ',
   'R-GPTSB'
    ))+
  scale_colour_manual(values = c('grey','#BA2D0B'),
                      name = '95%CI\nCrosses\nZero')+
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust=-1),
        strip.text.x = element_blank(),
        legend.position = 'none')

# Supplementary Figure 9 --------------------------------------------------

RT_anal <- intent_dat %>%
  mutate(group=ifelse(nchar(ID)<10, 'BPD', 'CON'),
         distancebeta=abs(server_beta_par-server_beta_ppt),
         distancealpha=abs(server_alpha_par-server_alpha_ppt),
         across(Gender:CTQtot, ~ifelse(.x==9999, NA, .x)),
         total_distance = distancebeta+distancealpha,
         total_distance_q = ifelse(total_distance < 24.75, 'low', 'high'),
         total_distance_q = factor(total_distance_q, levels = c('low', 'high'), ordered = T)
         ) %>%
  distinct()

ggplot(RT_anal %>%
         filter(RT<10000, Phase ==2),
       aes(trial, RT, colour= total_distance_q))+
  #stat_summary(geom='line', alpha = 0.01, aes(group=ID))+
  stat_summary(geom='line', alpha = 0.2, size=2)+
  geom_smooth(method='lm', formula = y ~ x + I(x^2))+
  scale_color_manual(name='Distance', values = c('#9A031E', 'black'))+
  coord_cartesian(ylim = c(1800, 4500))+
  stat_cor(label.y = c(4100, 3800), label.x = c(20), size=4, method = 'spearman')+
  labs(x='Trial',y='RT (ms)')+
  theme_bw(base_size=18)+
  theme(legend.position = 'none',
        panel.grid = element_blank())

#Phase 1
RT_model4 <- lme4::lmer(RT ~ group + (1|ID), data = RT_anal %>% filter(Phase%in%c(1)))
summary(RT_model4)
confint(RT_model4)
anova(RT_model4)

#Phase 2
RT_model1a <- lme4::lmer(RT ~ trial +  (1|ID), data = RT_anal %>% filter(Phase==2))
summary(RT_model1a)
confint(RT_model1a)

RT_model1b <- lme4::lmer(RT ~ total_distance +(1|ID), data = RT_anal %>% filter(Phase==2))
summary(RT_model1b)
confint(RT_model1b)

RT_model1c <- lme4::lmer(RT ~ trial:total_distance + (1|ID), data = RT_anal %>% filter(Phase==2))
summary(RT_model1c)
confint(RT_model1c)

RT_model2 <- lme4::lmer(RT ~ trial*group + (1|ID), data = RT_anal %>% filter(Phase==2))
summary(RT_model2)
confint(RT_model2)

#Phase 1-3
RT_model3a <- lme4::lmer(RT ~ Phase + group + (1|ID), data = RT_anal %>% filter(Phase%in%c(1,3)))
summary(RT_model3a)
confint(RT_model3a)

RT_model3b <- lme4::lmer(RT ~ group + (1|ID), data = RT_anal %>% filter(Phase%in%c(3)))
summary(RT_model3b)
confint(RT_model3b)
anova(RT_model3b)

# Supplementary Figure 10 --------------------------------------------------

library(rstanarm)
library(bayesplot)

## Load M1 only ----

M1_full <- readMat('FittedModelFiles/hbi_full_M1.mat')
M1_BPD  <- readMat('FittedModelFiles/hbi_g1_BPD_M1.mat')
M1_CON  <- readMat('FittedModelFiles/hbi_g2_BPD_M1.mat')

M1_BPD_sim  <- readMat('Data_Simulated/full_sim_g1_M1_only.mat')
M1_CON_sim  <- readMat('Data_Simulated/full_sim_g2_M1_only.mat')

M1_intent <- read_csv("Data/formatlab_BPD_allsamp.csv") %>%
  dplyr::select(-`...1`)

par_full <- as.data.frame(M1_full$cbm[,,1]$output[,,1]$parameters[1,])
colnames(par_full) <- c('alpha', 'beta', 'alpha_v', 'beta_v', 'alpha_ref', 'beta_ref')
par_full <- par_full %>%
  mutate(ID = unique(intent_dat_126$ID),
         type = 'full')

par_sep  <- rbind(as.data.frame(M1_BPD$cbm[,,1]$output[,,1]$parameters[1,]),
                  as.data.frame(M1_CON$cbm[,,1]$output[,,1]$parameters[1,]))
colnames(par_sep) <- c('alpha', 'beta', 'alpha_v', 'beta_v', 'alpha_ref', 'beta_ref')
par_sep <- par_sep %>%
  mutate(ID = c(ID1$ID, ID2$ID),
         type = 'sep')

M1_w_ID <- rbind(par_full,par_sep) %>%
           mutate(group = ifelse(nchar(ID)<10, 'BPD', 'CON'))

ggplot(M1_w_ID, aes(beta_ref, fill = group))+
  geom_density()+
  facet_wrap(~type)

# Process and summarize data
test_parms_m1 <- M1_w_ID %>%
  filter(type == 'sep') %>%
  mutate(across(alpha_v:beta_ref, ~ exp(.)),
         alpha = (1 / (1 + exp(-alpha))) * 30)

test_parms_m1 %>%
  group_by(group) %>%
  summarise(
    across(alpha:beta_ref,
           list(mean = ~ mean(.), sd = ~ sd(.)),
           .names = "{.col}_{.fn}")
  ) %>%
  as.data.frame()

summary(stan_glm(alpha_v ~ group, data = test_parms_m1 ))

## Belief Distributions ------

lim = 30
res = 0.25

for (i in 1:length(ID1$ID)){
  xdist <-   data.frame(
    bP1  = as.vector(M1_BPD_sim$results[i,][[1]][[1]][,,1]$beta.marg1),
    bP2a = as.vector(M1_BPD_sim$results[i,][[1]][[1]][,,1]$beta.marg2a),
    bP2b = as.vector(M1_BPD_sim$results[i,][[1]][[1]][,,1]$beta.marg2b),
    bP3  = as.vector(M1_BPD_sim$results[i,][[1]][[1]][,,1]$beta.marg3),
    aP1  = as.vector(M1_BPD_sim$results[i,][[1]][[1]][,,1]$alpha.marg1),
    aP2a = as.vector(M1_BPD_sim$results[i,][[1]][[1]][,,1]$alpha.marg2a),
    aP2b = as.vector(M1_BPD_sim$results[i,][[1]][[1]][,,1]$alpha.marg2b),
    aP3  = as.vector(M1_BPD_sim$results[i,][[1]][[1]][,,1]$alpha.marg3),
    beta = seq(-lim, lim, res),
    alpha = seq(0, lim, res/2),
    ID = ID1$ID[i],
    group = 'BPD',
    Model = 'M1'
  )
  if(i>1){BPDdist <- rbind(BPDdist, xdist)} else {BPDdist <- xdist}
}
for (i in 1:length(ID2$ID)){
  zdist <-   data.frame(
  bP1   = as.vector(M1_CON_sim$results[i,][[1]][[1]][,,1]$beta.marg1),
  bP2a  = as.vector(M1_CON_sim$results[i,][[1]][[1]][,,1]$beta.marg2a),
  bP2b  = as.vector(M1_CON_sim$results[i,][[1]][[1]][,,1]$beta.marg2b),
  bP3   = as.vector(M1_CON_sim$results[i,][[1]][[1]][,,1]$beta.marg3),
  aP1   = as.vector(M1_CON_sim$results[i,][[1]][[1]][,,1]$alpha.marg1),
  aP2a  = as.vector(M1_CON_sim$results[i,][[1]][[1]][,,1]$alpha.marg2a),
  aP2b  = as.vector(M1_CON_sim$results[i,][[1]][[1]][,,1]$alpha.marg2b),
  aP3   = as.vector(M1_CON_sim$results[i,][[1]][[1]][,,1]$alpha.marg3),
  beta  = seq(-lim, lim, res),
  alpha = seq(0, lim, res/2),
  ID = ID2$ID[i],
  group = 'CON',
  Model = 'M1'
  )
  if(i>1){CONdist <- rbind(CONdist, zdist)} else {CONdist <- zdist}
}

int_dist_m1 <- rbind(BPDdist, CONdist)

par_parms_m1 <- int_dist_m1 %>%
  group_by(ID) %>%
  mutate(beta_par_m = sum(beta*bP2b),
         alpha_par_m = sum(alpha*aP2b),
         beta_par_pri = sum(beta*bP2a),
         alpha_par_pri = sum(alpha*aP2a),
         B1_sd = sqrt(sum((beta-bP1)^2 * bP1)),
         B2a_sd= sqrt(sum((beta-bP2a)^2 * bP2a)),
         B2b_sd= sqrt(sum((beta-bP2b)^2 * bP2b)),
         B3_sd = sqrt(sum((beta-bP3)^2 * bP3)),
         A1_sd = sqrt(sum((alpha-aP1)^2* aP1)),
         A2a_sd= sqrt(sum((alpha-aP2a)^2 *aP2a)),
         A2b_sd= sqrt(sum((alpha-aP2b)^2 * aP2b)),
         A3_sd = sqrt(sum((alpha-aP3)^2* aP3)),
         shift_beta_p2 = abs(beta_par_pri-beta_par_m),
         shift_alpha_p2 = abs(alpha_par_pri-alpha_par_m),
         beta_hat = sum(beta*bP3),
         alpha_hat= sum(alpha*aP3),
         beta_hat_sd = sqrt(sum((beta-beta_hat)^2 * bP3)),
         alpha_hat_sd = sqrt(sum((alpha-alpha_hat)^2 * aP3))) %>%
  dplyr::select(ID, group,
                beta_par_m, alpha_par_m,
                beta_par_pri, alpha_par_pri,
                beta_hat, alpha_hat,
                beta_hat_sd, alpha_hat_sd,
                shift_beta_p2, shift_alpha_p2,
                B2a_sd, B2b_sd, A2a_sd, A2b_sd) %>%
  distinct() %>%
  plyr::join(., intent_dat_126 %>%
         dplyr::select(ID, server_beta_par, server_alpha_par)) %>%
  mutate(disp_beta_par = abs(server_beta_par - beta_par_m),
         disp_alpha_par = abs(server_alpha_par - alpha_par_m)) %>%
  dplyr::select(ID,
                beta_par_m, alpha_par_m,
                beta_par_pri, alpha_par_pri,
                disp_beta_par, disp_alpha_par,
                beta_hat_sd, alpha_hat_sd,
                beta_hat, alpha_hat,
                shift_beta_p2, shift_alpha_p2,
                B2a_sd, B2b_sd, A2a_sd, A2b_sd) %>%
  distinct()

M1_only_parameters <- plyr::join(
  M1_w_ID %>% filter(type =='sep'),
  par_parms_m1,
  by = 'ID'
  )

## Mean val shift -----

mean_val_shift_m1 <- int_dist_m1 %>%
  group_by(ID) %>%
  mutate(B1  = sum(beta*bP1),
         B2a = sum(beta*bP2a),
         B2b = sum(beta*bP2b),
         B3  = sum(beta*bP3),
         A1  = sum(alpha*aP1),
         A2a = sum(alpha*aP2a),
         A2b = sum(alpha*aP2b),
         A3  = sum(alpha*aP3),
         B1_sd = sqrt(sum((beta-bP1)^2 * bP1)),
         B2a_sd= sqrt(sum((beta-bP2a)^2 * bP2a)),
         B2b_sd= sqrt(sum((beta-bP2b)^2 * bP2b)),
         B3_sd = sqrt(sum((beta-bP3)^2 * bP3)),
         A1_sd = sqrt(sum((alpha-aP1)^2* aP1)),
         A2a_sd= sqrt(sum((alpha-aP2a)^2 *aP2a)),
         A2b_sd= sqrt(sum((alpha-aP2b)^2 * aP2b)),
         A3_sd = sqrt(sum((alpha-aP3)^2* aP3)),
         delta_b = abs(B3-B1),
         delta_a = abs(A3-A1),
         delta_b_sd = abs(B3_sd-B1_sd),
         delta_a_sd = abs(A3_sd-A1_sd)) %>%
  dplyr::select(ID: delta_a %>%
  distinct()

sd(mean_val_shift_m1[mean_val_shift_m1$group=='CON',]$delta_b_sd)

summary(stan_glm(delta_b~group, data = mean_val_shift_m1))
summary(stan_glm(delta_a~group, data = mean_val_shift_m1))
summary(stan_glm(delta_b_sd~group, data = mean_val_shift_m1))
summary(stan_glm(delta_a_sd~group, data = mean_val_shift_m1))

## Value shift M4$2a - M4$2b, M1$ref ---------------------------------------------

beta_cor_test <- mean_val_shift %>%
  filter(group=='BPD') %>%
  ungroup() %>%
  plyr::join(.,
             M1_only_parameters %>%
               dplyr::select(ID, beta_ref, alpha_ref),
             by = 'ID')

cor.test(exp(beta_cor_test$beta_ref), abs(beta_cor_test$B1-beta_cor_test$B2a))
cor.test(exp(beta_cor_test$alpha_ref), abs(beta_cor_test$A1-beta_cor_test$A2a))

## Replicate Phase 2 analysis -----------

# Process and summarize data
M1_only_parameters %>%
  filter(type == 'sep') %>%
  group_by(group) %>%
  summarise(
    across(c(shift_beta_p2, shift_alpha_p2,A2b_sd, B2b_sd),
           list(mean = ~ mean(.), sd = ~ sd(.)),
           .names = "{.col}_{.fn}")
  ) %>%
  as.data.frame()

summary(stan_glm(shift_beta_p2 ~ group, M1_only_parameters, family = gaussian()))
summary(stan_glm(shift_alpha_p2 ~ group, M1_only_parameters, family = gaussian()))

summary(stan_glm(B2b_sd ~ group, M1_only_parameters, family = gaussian()))
summary(stan_glm(A2b_sd ~ group, M1_only_parameters, family = gaussian()))

## Belief Updates In Phase 2 ------------

b_updates_m1 <- intent_dat_126 %>%
  filter(Phase==2) %>%
  dplyr::select(ID, trial, correctSum) %>%
  mutate(group = ifelse(nchar(ID) < 10, 'BPD', 'CON'),
         kl_div_a = 0,
         kl_div_b = 0)

#sanity check for kl divs
beta_check_m1 <- data.frame(
  beta = seq(-30, 30, 0.25),
  t1 = as.numeric(M1_BPD_sim$results[5,][[1]][[1]][,,1]$beta.marg1),
  t2 = as.numeric(M1_BPD_sim$results[5,][[1]][[1]][,,1]$beta.cont[,38])
)

ggplot(beta_check_m1, aes(beta, t1))+
  geom_line(colour = 'red')+
  geom_line(data = beta_check_m1, aes(beta, t2))

#loops for all
for(k in 1:2){
  if(k == 1){x = M1_BPD_sim$results; group = 'BPD'; ID = ID1$ID}
  if(k == 2){x = M1_CON_sim$results; group = 'CON'; ID = ID2$ID}
  for(j in 1:length(ID)){
    kl_divsa = rep(NA, 54)
    kl_divsb = kl_divsa
    for(i in 1:54){
      x[j,][[1]][[1]][,,1]$alpha.cont[,36] = x[j,][[1]][[1]][,,1]$alpha.marg1
      x[j,][[1]][[1]][,,1]$beta.cont[,36]  = x[j,][[1]][[1]][,,1]$beta.marg1

      bsa            = x[j,][[1]][[1]][,,1]$alpha.cont[,36:90]
      b_t2a          = bsa[,i+1]
      b_t1a          = bsa[,i]
      bsb            = x[j,][[1]][[1]][,,1]$beta.cont[,36:90]
      b_t2b          = bsb[,i+1]
      b_t1b          = bsb[,i]
      kl_divsa[i]    = calculate_KL_divergence(b_t2a, b_t1a)
      kl_divsb[i]    = calculate_KL_divergence(b_t2b, b_t1b)
    }
    b_updates_m1[b_updates_m1$group==group & b_updates_m1$ID==ID[j],'kl_div_a'][1:54,] <- kl_divsa
    b_updates_m1[b_updates_m1$group==group & b_updates_m1$ID==ID[j],'kl_div_b'][1:54,] <- kl_divsb
  }
}

b_updates_m1 <- b_updates_m1 %>%
  mutate(kl_div_aroll = rollapply(kl_div_a, width = 5, FUN = mean, align = "right", fill = NA),
         kl_div_broll = rollapply(kl_div_b, width = 5, FUN = mean, align = "right", fill = NA))

facet_labels <- c('alpha', 'beta')
names(facet_labels) <- c('kl_div_a', 'kl_div_b')
ggplot(b_updates_m1 %>%
         filter(trial != 1) %>%
         pivot_longer(5:6, names_to = 'Parameter', values_to = 'KL_Div'),
       aes(trial, KL_Div, fill = group, colour = group))+
  geom_smooth(alpha = 0.2)+
  scale_fill_manual(values=colour_group)+
  scale_colour_manual(values=colour_group)+
  scale_x_continuous(breaks = c(2, 25, 50))+
  facet_wrap(~Parameter, nrow=2) +
  labs(y = expression(paste('D'[KL])), x = 'Trial')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=18),
        legend.position = 'none',
        legend.title = element_blank(),
        strip.text.x = element_blank())

summary(lmerTest::lmer(scale(kl_div_a) ~ trial*group+(1|ID), data = b_updates_m1 %>% filter(trial%in%c(2:54))))
summary(lmerTest::lmer(scale(kl_div_b) ~ trial*group+(1|ID), data = b_updates_m1 %>% filter(trial%in%c(2:54))))

summary(stan_glm(scale(kl) ~ trial*type,
           data = b_updates_m1 %>%
             filter(trial%in%c(2:54)) %>%
             pivot_longer(kl_div_a:kl_div_b, names_to = 'type', values_to = 'kl')
           ))

## Groups fitted together ----

sanity_beta_ref_tog <- stan_glm(exp(beta_ref) ~ group,
                                M1_w_ID %>% filter(type == 'full'),
                                family = gaussian())

summary(sanity_beta_ref_tog)
sanity_beta_post_tog <- as.matrix(sanity_beta_ref_tog) %>%
  as.data.frame() %>%
  mutate(type = 'FULL')

sanity_alpha_ref_tog <- stan_glm(exp(alpha_ref) ~ group,
                                M1_w_ID %>% filter(type == 'full'),
                                family = gaussian())

summary(sanity_alpha_ref_tog)
sanity_alpha_post_tog <- as.matrix(sanity_alpha_ref_tog) %>%
  as.data.frame() %>%
  mutate(type = 'FULL')

## Groups fitted seperately ----

sanity_beta_ref_sep <- stan_glm(exp(beta_ref) ~ group,
                                M1_w_ID %>% filter(type =='sep'),
                                family = gaussian())

summary(sanity_beta_ref_sep)
sanity_beta_post_sep  <- as.matrix(sanity_beta_ref_sep)%>%
  as.data.frame() %>%
  mutate(type = 'SEP')

sanity_alpha_ref_sep <- stan_glm(exp(alpha_ref) ~ group,
                                M1_w_ID %>% filter(type == 'sep'),
                                family = gaussian())

summary(sanity_alpha_ref_sep)
sanity_alpha_post_sep <- as.matrix(sanity_alpha_ref_sep) %>%
  as.data.frame() %>%
  mutate(type = 'SEP')

## Plot ----

library(ggridges)
library(tidybayes)

rbind(sanity_beta_post_tog,
      sanity_beta_post_sep) %>%
  ggplot(aes(x = groupCON, y = type)) +
  stat_interval(.width = c(.90, .8, .5)) +
  ggdist::stat_halfeye(.width = c(.90, .8, .5))+
  geom_vline(xintercept = 0, size = 1)+
  scale_fill_brewer() +
  labs(x = expression(paste(Delta, mu, '[CON-BPD]'))) +
  coord_cartesian(xlim = c(-6, 4))+
  theme_bw(base_size = 24) +
  theme(
    panel.border  = element_blank(),
    legend.position = 'none',
    panel.grid    = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.title.y  = element_blank()
  ) +

ggplot(M1_w_ID_diff, aes(x = exp(full), y = exp(sep))) +
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(x = expression(paste(beta[par]^ref, ' [FULL]')),
       y = expression(paste(beta[par]^ref, ' [SEP]')), ) +
  theme_bw(base_size = 24) +
  theme(
    panel.border  = element_blank(),
    legend.position = 'none',
    panel.grid    = element_blank(),
    axis.ticks.y  = element_blank()
  ) +

rbind(sanity_alpha_post_tog,
      sanity_alpha_post_sep) %>%
  ggplot(aes(x = groupCON, y = type)) +
  stat_interval(.width = c(.90, .8, .5)) +
  ggdist::stat_halfeye(.width = c(.90, .8, .5))+
  geom_vline(xintercept = 0, size = 1)+
  scale_fill_brewer() +
  labs(x = expression(paste(Delta, mu, '[CON-BPD]'))) +
  coord_cartesian(xlim = c(-6, 4))+
  theme_bw(base_size = 24) +
  theme(
    panel.border  = element_blank(),
    legend.position = 'none',
    panel.grid    = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.title.y  = element_blank()
  ) +

ggplot(M1_w_ID_diff_a, aes(x = exp(full), y = exp(sep))) +
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(x = expression(paste(alpha[par]^ref, ' [FULL]')),
       y = expression(paste(alpha[par]^ref, ' [SEP]')), ) +
  theme_bw(base_size = 24) +
  theme(
    panel.border  = element_blank(),
    legend.position = 'none',
    panel.grid    = element_blank(),
    axis.ticks.y  = element_blank()
  )

## Decide whether beta or alpha more important ----

#first check within M4 to see which median difference is largest
stan_glm(val ~ parm,
         joint_parms_ex %>%
           pivot_longer(c(beta, beta_par), names_to = 'parm', values_to = 'val') %>%
           filter(group == 'BPD_full') %>%
           dplyr::select(parm, val, ID) %>%
           distinct(),
         family = gaussian()) %>%
  summary()

stan_glm(val ~ parm,
         joint_parms_ex %>%
           pivot_longer(c(alpha, alpha_par), names_to = 'parm', values_to = 'val') %>%
           filter(group == 'BPD_full') %>%
           dplyr::select(parm, val, ID) %>%
           distinct(),
         family = gaussian()) %>%
  summary()

#then check in M1 forced fits if alpha_ref value is larger
stan_glm(alpha_ref ~ group,
         M1_w_ID %>% filter(type == 'full'),
         family = gaussian()) %>%
  summary()

stan_glm(alpha_ref ~ group,
         M1_w_ID %>% filter(type == 'sep'),
         family = gaussian()) %>%
  summary()


top_median_parm <- joint_parms_ex %>%
  filter(group=='BPD_full') %>%
  ggplot()+
  geom_density(aes(beta), fill = 'darkred', alpha = 0.5)+
  stat_central_tendency(aes(beta), type = 'median', colour = 'darkred')+
  geom_density(aes(beta_par), fill = 'grey', alpha = 0.5)+
  stat_central_tendency(aes(beta_par), type = 'median', colour = 'black')+
  labs(x = expression(beta))+
  coord_cartesian()
bottom_median_parm <- joint_parms_ex %>%
  filter(group=='BPD_full') %>%
  ggplot()+
  geom_density(aes(alpha), fill = 'darkred', alpha = 0.5)+
  stat_central_tendency(aes(alpha), type = 'median', colour = 'darkred')+
  geom_density(aes(alpha_par), fill = 'grey', alpha = 0.5)+
  stat_central_tendency(aes(alpha_par), type = 'median', colour = 'black')+
  labs(x = expression(alpha))+
  coord_cartesian()

top_median_parm/bottom_median_parm&
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank())

# Supplementary Figure 11 --------------------------------------------------

# Initialize an empty list to store processed data frames
model_data_list <- list()

# Loop through models M1-M4
for (model in c("M1", "M2", "M3", "M4")) {
  # Load .mat files for BPD and CON groups
  BPD_mat <- readMat(paste0('MatlabFitting/HBI_fullBPD/hbi_g1_BPD_', model, '.mat'))
  CON_mat <- readMat(paste0('MatlabFitting/HBI_fullBPD/hbi_g2_BPD_', model, '.mat'))

  # Extract and bind parameters from both groups
  model_bind <- rbind(
    as.data.frame(BPD_mat$cbm[,,1]$output[,,1]$parameters[1,]),
    as.data.frame(CON_mat$cbm[,,1]$output[,,1]$parameters[1,])
  )

  # Rename columns
  if (model %in% c("M1", "M2")) {
    # Add placeholder columns for M1 and M2
    model_bind <- model_bind %>%
      mutate(alpha_par = NA, beta_par = NA) %>%
      dplyr::select(1:4, alpha_par, beta_par, everything())
  }

  colnames(model_bind) <- c('alpha', 'beta', 'alpha_v', 'beta_v',
                            'alpha_par', 'beta_par', 'alpha_ref', 'beta_ref')

  # Add ID and group columns
  model_bind <- model_bind %>%
    mutate(ID = c(ID1$ID, ID2$ID),
           group = ifelse(nchar(ID) < 10, 'BPD', 'CON'),
           Model = model) # Add the model name (M1-M4)

  # Append to the list
  model_data_list[[model]] <- model_bind
}

# Combine all models into a single data frame
final_data <- bind_rows(model_data_list)

## Correlations between parameters -----

# Create a helper function to calculate correlations
calculate_correlation <- function(x, y, label_x, label_y) {
  cor_result <- cor.test(x, y)
  data.frame(
    x_label = label_x,
    y_label = label_y,
    R = cor_result$estimate,
    p_value = cor_result$p.value
  )
}

# Prepare data for correlations
correlations <- rbind(
  calculate_correlation(
    abs(final_data[final_data$Model == "M4", ]$beta -
          final_data[final_data$Model == "M4", ]$beta_par),
    exp(final_data[final_data$Model == "M1", ]$beta_ref),
    "abs(beta_M4 - beta_par_M4)", "exp(beta_ref_M1)"
  ),
  calculate_correlation(
    abs(final_data[final_data$Model == "M4", ]$alpha -
          final_data[final_data$Model == "M4", ]$alpha_par),
    exp(final_data[final_data$Model == "M1", ]$alpha_ref),
    "abs(alpha_M4 - alpha_par_M4)", "exp(alpha_ref_M1)"
  ),
  calculate_correlation(
    abs(final_data[final_data$Model == "M3", ]$beta -
          final_data[final_data$Model == "M3", ]$beta_par),
    exp(final_data[final_data$Model == "M1", ]$beta_ref),
    "abs(beta_M3 - beta_par_M3)", "exp(beta_ref_M1)"
  ),
  calculate_correlation(
    abs(final_data[final_data$Model == "M3", ]$alpha -
          final_data[final_data$Model == "M3", ]$alpha_par),
    exp(final_data[final_data$Model == "M1", ]$alpha_ref),
    "abs(alpha_M3 - alpha_par_M3)", "exp(alpha_ref_M1)"
  ),
  calculate_correlation(
    abs(final_data[final_data$Model == "M4", ]$beta -
          final_data[final_data$Model == "M4", ]$beta_par),
    exp(final_data[final_data$Model == "M2", ]$beta_ref),
    "abs(beta_M4 - beta_par_M4)", "exp(beta_ref_M2)"
  ),
  calculate_correlation(
    abs(final_data[final_data$Model == "M4", ]$alpha -
          final_data[final_data$Model == "M4", ]$alpha_par),
    exp(final_data[final_data$Model == "M2", ]$alpha_ref),
    "abs(alpha_M4 - alpha_par_M4)", "exp(alpha_ref_M2)"
  ),
  calculate_correlation(
    abs(final_data[final_data$Model == "M3", ]$beta -
          final_data[final_data$Model == "M3", ]$beta_par),
    exp(final_data[final_data$Model == "M2", ]$beta_ref),
    "abs(beta_M3 - beta_par_M3)", "exp(beta_ref_M2)"
  ),
  calculate_correlation(
    abs(final_data[final_data$Model == "M3", ]$alpha -
          final_data[final_data$Model == "M3", ]$alpha_par),
    exp(final_data[final_data$Model == "M2", ]$alpha_ref),
    "abs(alpha_M3 - alpha_par_M3)", "exp(alpha_ref_M2)"
  ),
  calculate_correlation(
    exp(final_data[final_data$Model == "M2", ]$beta_ref),
    exp(final_data[final_data$Model == "M1", ]$beta_ref),
    "beta_ref M2", "beta_ref M1"
  ),
  calculate_correlation(
    final_data[final_data$Model == "M3", ]$beta_par,
    final_data[final_data$Model == "M4", ]$beta_par,
    "beta_par M3", "beta_par M4"
  ),
  calculate_correlation(
    exp(final_data[final_data$Model == "M2", ]$alpha_ref),
    exp(final_data[final_data$Model == "M1", ]$alpha_ref),
    "alpha_ref M2", "alpha_ref M1"
  ),
  calculate_correlation(
    final_data[final_data$Model == "M3", ]$alpha_par,
    final_data[final_data$Model == "M4", ]$alpha_par,
    "alpha_par M3", "alpha_par M4"
  ),
  calculate_correlation(
    abs(mean_val_shift_M3[mean_val_shift_M3$Delta=='Delta_A',]$Shift),
    mean_val_shift_m1$delta_a,
    "delta a M3", "delta a M1"
  ),
  calculate_correlation(
    abs(mean_val_shift_M3[mean_val_shift_M3$Delta=='Delta_B',]$Shift),
    mean_val_shift_m1$delta_b,
    "delta b M3", "delta b M1"
  )
)

# Plot the correlations
ggplot(correlations[1:8,], aes(x = x_label, y = y_label, fill = R)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(0, 1), name = "Correlation (R)"
  ) +
  geom_text(aes(label = sprintf("R=%.2f\np=%.2f", R, p_value)), color = "black", size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    text = element_text(size = 18),
    legend.position = 'top',
    legend.key.width = unit(1.5, 'cm'),
    legend.key.spacing = unit(0.5, 'cm'),
    legend.justification.top = c(0,0)
  ) +
  scale_x_discrete(labels = c(
    expression(paste(alpha[ppt]^m,'-',alpha[par]^m, '(M3)')),
    expression(paste(alpha[ppt]^m,'-',alpha[par]^m, '(M4)')),
    expression(paste(beta[ppt]^m,'-',beta[par]^m, '(M3)')),
    expression(paste(beta[ppt]^m,'-',beta[par]^m, '(M4)'))
  ))+
  scale_y_discrete(labels = c(
    expression(paste(alpha[ppt]^ref, '(M1)')),
    expression(paste(alpha[ppt]^ref, '(M2)')),
    expression(paste(beta[ppt]^ref,  '(M1)')),
    expression(paste(beta[ppt]^ref,  '(M2)'))
  ))+

ggplot(correlations[9:14,], aes(x = x_label, y = y_label, fill = R)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(0, 1), name = "Correlation (R)"
  ) +
  geom_text(aes(label = sprintf("R=%.2f\np=%.2f", R, p_value)), color = "black", size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    text = element_text(size = 18),
    legend.key.width = unit(1, 'cm'),
    legend.position = 'none'
  ) +
  scale_x_discrete(labels = c(
    expression(paste(alpha[par]^m, '(M3)')),
    expression(paste(alpha[ppt]^ref, '(M2)')),
    expression(paste(beta[par]^m, '(M3)')),
    expression(paste(beta[ppt]^ref, '(M2)')),
    expression(paste(Delta, alpha[ppt]^m, '(M3)')),
    expression(paste(Delta, beta[ppt]^m, '(M3)'))
  ))+
  scale_y_discrete(labels = c(
    expression(paste(alpha[par]^m, '(M4)')),
    expression(paste(alpha[ppt]^ref, '(M1)')),
    expression(paste(beta[par]^m, '(M4)')),
    expression(paste(beta[ppt]^ref, '(M1)')),
    expression(paste(Delta, alpha[ppt]^m, '(M1)')),
    expression(paste(Delta, beta[ppt]^m, '(M1)'))
  ))


# Supplementary Table S2 --------------------------------------------------

## Directly Approximated -------------------------------------------------

# Initialize an empty list to store processed data frames
model_data_list <- list()

# Loop through models M1-M4
for (model in c("M1", "M2", "M3", "M4")) {
  # Load .mat files for BPD and CON groups
  BPD_mat <- readMat(paste0('FittedModelFiles/hbi_g1_BPD_', model, '.mat'))
  CON_mat <- readMat(paste0('FittedModelFiles/hbi_g2_BPD_', model, '.mat'))

  # Extract and bind parameters from both groups
  model_bind <- rbind(
    as.data.frame(BPD_mat$cbm[,,1]$output[,,1]$parameters[1,]),
    as.data.frame(CON_mat$cbm[,,1]$output[,,1]$parameters[1,])
  )

  # Rename columns
  if (model %in% c("M1", "M2")) {
    # Add placeholder columns for M1 and M2
    model_bind <- model_bind %>%
      mutate(alpha_par = NA, beta_par = NA) %>%
      dplyr::select(1:4, alpha_par, beta_par, everything())
  }

  colnames(model_bind) <- c('alpha', 'beta', 'alpha_v', 'beta_v',
                            'alpha_par', 'beta_par', 'alpha_ref', 'beta_ref')

  # Add ID and group columns
  model_bind <- model_bind %>%
    mutate(ID = c(ID1$ID, ID2$ID),
           group = ifelse(nchar(ID) < 10, 'BPD', 'CON'),
           Model = model) # Add the model name (M1-M4)

  # Append to the list
  model_data_list[[model]] <- model_bind
}

# Combine all models into a single data frame
final_data <- bind_rows(model_data_list)

final_data %>%
  group_by(Model, group) %>%
  filter(Model != 'M1') %>%
  summarise(
    alpha_mean = mean((1 / (1 + exp(-alpha))) * 30, na.rm = TRUE),
    alpha_sd = sd((1 / (1 + exp(-alpha))) * 30, na.rm = TRUE),
    beta_mean = mean(beta, na.rm = TRUE),
    beta_sd = sd(beta, na.rm = TRUE),
    alpha_v_mean = mean(exp(alpha_v), na.rm = TRUE),
    alpha_v_sd = sd(exp(alpha_v), na.rm = TRUE),
    beta_v_mean = mean(exp(beta_v), na.rm = TRUE),
    beta_v_sd = sd(exp(beta_v), na.rm = TRUE),
    alpha_par_mean = mean((1 / (1 + exp(-alpha_par))) * 30, na.rm = TRUE),
    alpha_par_sd = sd((1 / (1 + exp(-alpha_par))) * 30, na.rm = TRUE),
    beta_par_mean = mean(beta_par, na.rm = TRUE),
    beta_par_sd = sd(beta_par, na.rm = TRUE),
    alpha_ref_mean = mean(exp(alpha_ref), na.rm = TRUE),
    alpha_ref_sd = sd(exp(alpha_ref), na.rm = TRUE),
    beta_ref_mean = mean(exp(beta_ref), na.rm = TRUE),
    beta_ref_sd = sd(exp(beta_ref), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  filter(Model == 'M4') %>%
  t()

summary(stan_glm(exp(beta_ref) ~ group, data = final_data %>% filter(Model=='M4')))

## Derived -------------------------------------------------

# Define constants
lim <- 30
res <- 0.25

# Initialize an empty list to store model data
all_model_data_mean_val <- list()
LL_sum <- matrix(NA, nrow = 4, ncol = 2)

# Loop through models M1 to M4
models <- c("M1", "M2", "M3", "M4")
rownames(LL_sum) <- models

for (model in models) {
  # Load corresponding .mat files for BPD and CON
  BPD_sim_sep <- readMat(paste0('Data_Simulated/full_sim_g1_', model, '_only.mat'))
  CON_sim_sep <- readMat(paste0('Data_Simulated/full_sim_g2_', model, '_only.mat'))

  # Initialize data frames for this model
  BPDdist <- NULL
  CONdist <- NULL

  # Process BPD Group
  for (i in 1:length(ID1)) {
    xdist <- data.frame(
      bP1  = as.vector(BPD_sim_sep$results[i,][[1]][[1]][,,1]$beta.marg1),
      bP2a = as.vector(BPD_sim_sep$results[i,][[1]][[1]][,,1]$beta.marg2a),
      bP2b = as.vector(BPD_sim_sep$results[i,][[1]][[1]][,,1]$beta.marg2b),
      bP3  = as.vector(BPD_sim_sep$results[i,][[1]][[1]][,,1]$beta.marg3),
      aP1  = as.vector(BPD_sim_sep$results[i,][[1]][[1]][,,1]$alpha.marg1),
      aP2a = as.vector(BPD_sim_sep$results[i,][[1]][[1]][,,1]$alpha.marg2a),
      aP2b = as.vector(BPD_sim_sep$results[i,][[1]][[1]][,,1]$alpha.marg2b),
      aP3  = as.vector(BPD_sim_sep$results[i,][[1]][[1]][,,1]$alpha.marg3),
      beta = seq(-lim, lim, res),
      alpha = seq(0, lim, res/2),
      ID = ID1[i],
      group = 'BPD',
      Model = model
    )
    BPDdist <- if (is.null(BPDdist)) xdist else rbind(BPDdist, xdist)
  }

  # Process CON Group
  for (i in 1:length(ID2)) {
    zdist <- data.frame(
      bP1  = as.vector(CON_sim_sep$results[i,][[1]][[1]][,,1]$beta.marg1),
      bP2a = as.vector(CON_sim_sep$results[i,][[1]][[1]][,,1]$beta.marg2a),
      bP2b = as.vector(CON_sim_sep$results[i,][[1]][[1]][,,1]$beta.marg2b),
      bP3  = as.vector(CON_sim_sep$results[i,][[1]][[1]][,,1]$beta.marg3),
      aP1  = as.vector(CON_sim_sep$results[i,][[1]][[1]][,,1]$alpha.marg1),
      aP2a = as.vector(CON_sim_sep$results[i,][[1]][[1]][,,1]$alpha.marg2a),
      aP2b = as.vector(CON_sim_sep$results[i,][[1]][[1]][,,1]$alpha.marg2b),
      aP3  = as.vector(CON_sim_sep$results[i,][[1]][[1]][,,1]$alpha.marg3),
      beta = seq(-lim, lim, res),
      alpha = seq(0, lim, res/2),
      ID = ID2[i],
      group = 'CON',
      Model = model
    )
    CONdist <- if (is.null(CONdist)) zdist else rbind(CONdist, zdist)
  }

  # Combine BPD and CON data for this model
  int_dist_all_mod <- rbind(BPDdist, CONdist)

  # Compute mean and standard deviations for key variables
  mean_val_shift_all_mod <- int_dist_all_mod %>%
    group_by(ID) %>%
    mutate(
      B1  = sum(beta * bP1),
      B2a = sum(beta * bP2a),
      B2b = sum(beta * bP2b),
      B3  = sum(beta * bP3),
      A1  = sum(alpha * aP1),
      A2a = sum(alpha * aP2a),
      A2b = sum(alpha * aP2b),
      A3  = sum(alpha * aP3),
      B1_sd = sqrt(sum((beta - B1)^2 * bP1)),
      B2a_sd = sqrt(sum((beta - B2a)^2 * bP2a)),
      B2b_sd = sqrt(sum((beta - B2b)^2 * bP2b)),
      B3_sd = sqrt(sum((beta - B3)^2 * bP3)),
      A1_sd = sqrt(sum((alpha - A1)^2 * aP1)),
      A2a_sd = sqrt(sum((alpha - A2a)^2 * aP2a)),
      A2b_sd = sqrt(sum((alpha - A2b)^2 * aP2b)),
      A3_sd = sqrt(sum((alpha - A3)^2 * aP3)),
      delta_a_par = abs(A2b - A2a),
      delta_b_par = abs(B2b - B2a),
      delta_b = abs(B3 - B1),
      delta_a = abs(A3 - A1),
      delta_b_sd = abs(B3_sd - B1_sd),
      delta_a_sd = abs(A3_sd - A1_sd)
    ) %>%
    dplyr::select(ID, group, Model, B1:delta_a_sd) %>%
    distinct()

  # Store results for this model
  all_model_data_mean_val[[model]] <- mean_val_shift_all_mod

}

# Combine data across all models
combined_all_mod_mean_val <- bind_rows(all_model_data_mean_val)

combined_all_mod_mean_val %>%
  filter(Model == 'M3') %>%
  group_by(group) %>%
  summarise(mean = mean(delta_a_par),
            sd = sd(delta_a_par))

rstanarm::stan_glm(B2b_sd ~ group,
                   data = combined_all_mod_mean_val %>%
                     filter(Model == 'M3')) %>%
  summary()

summary(lm(Shift ~ group * Delta,
           data = combined_all_mod_mean_val %>%
                     filter(Model == 'M3') %>%
             pivot_longer(delta_b_sd:delta_a_sd, names_to = 'Delta', values_to = 'Shift') %>%
             distinct()))

summary_data <- combined_all_mod_mean_val %>%
  group_by(Model) %>%
  summarise(
    mean_bP2b_BPD = mean(B2b_sd[group == "BPD"], na.rm = TRUE),
    sd_bP2b_BPD = sd(B2b_sd[group == "BPD"], na.rm = TRUE),
    mean_bP2b_CON = mean(B2b_sd[group == "CON"], na.rm = TRUE),
    sd_bP2b_CON = sd(B2b_sd[group == "CON"], na.rm = TRUE),
    ttest_bP2b_p = t.test(B2b_sd[group == "BPD"], B2b_sd[group == "CON"])$p.value,

    mean_aP2b_BPD = mean(A2b_sd[group == "BPD"], na.rm = TRUE),
    sd_aP2b_BPD = sd(A2b_sd[group == "BPD"], na.rm = TRUE),
    mean_aP2b_CON = mean(A2b_sd[group == "CON"], na.rm = TRUE),
    sd_aP2b_CON = sd(A2b_sd[group == "CON"], na.rm = TRUE),
    ttest_aP2b_p = t.test(A2b_sd[group == "BPD"], A2b_sd[group == "CON"])$p.value,

    .groups = "drop"
  )

# Print summary
print(summary_data)

# Supplementary Table S3 -------------------------------------

regression_results <- list()

# Define the models to iterate over
models <- c("M1", "M2", "M3", "M4")

# Loop through models M1 to M4
for (model in models) {
  # Load corresponding .mat files for BPD and CON
  BPD_sim_sep <- readMat(paste0('Data_Simulated/full_sim_g1_', model, '_only.mat'))
  CON_sim_sep <- readMat(paste0('Data_Simulated/full_sim_g2_', model, '_only.mat'))

  b_updates_sep <- intent_dat %>%
    filter(Phase==2) %>%
    dplyr::select(ID, trial, correctSum) %>%
    mutate(group = ifelse(nchar(ID) < 10, 'BPD', 'CON'),
           kl_div_a = 0,
           kl_div_b = 0)

  #loops for all
  for(k in 1:2){
    if(k == 1){x = BPD_sim_sep$results; group = 'BPD'; ID = ID1}
    if(k == 2){x = CON_sim_sep$results; group = 'CON'; ID = ID2}
    for(j in 1:length(ID)){
      kl_divsa = rep(NA, 54)
      kl_divsb = kl_divsa
      for(i in 1:54){
        x[j,][[1]][[1]][,,1]$alpha.cont[,36] = x[j,][[1]][[1]][,,1]$alpha.marg1
        x[j,][[1]][[1]][,,1]$beta.cont[,36]  = x[j,][[1]][[1]][,,1]$beta.marg1

        bsa            = x[j,][[1]][[1]][,,1]$alpha.cont[,36:90]
        b_t2a          = bsa[,i+1]
        b_t1a          = bsa[,i]
        bsb            = x[j,][[1]][[1]][,,1]$beta.cont[,36:90]
        b_t2b          = bsb[,i+1]
        b_t1b          = bsb[,i]
        kl_divsa[i]    = calculate_KL_divergence(b_t2a, b_t1a)
        kl_divsb[i]    = calculate_KL_divergence(b_t2b, b_t1b)
      }
      b_updates_sep[b_updates_sep$group==group & b_updates_sep$ID==ID[j],'kl_div_a'][1:54] <- kl_divsa
      b_updates_sep[b_updates_sep$group==group & b_updates_sep$ID==ID[j],'kl_div_b'][1:54] <- kl_divsb
    }
  }

    model_kl_grp_a   <- summary(lmerTest::lmer(scale(kl_div_a) ~ group + (1 | ID),
                                     data = b_updates_sep %>% filter(trial %in% c(2:54))))
    model_kl_grp_b   <- summary(lmerTest::lmer(scale(kl_div_b) ~ group + (1 | ID),
                                       data = b_updates_sep %>% filter(trial %in% c(2:54))))
    b_updates_sep <- b_updates_sep %>%
      pivot_longer(kl_div_a:kl_div_b, names_to = 'kl', values_to = 'div')

    # Regression analysis for KL divergence (alpha)
    model_kl_trial   <- summary(lmerTest::lmer(scale(div) ~ trial + (1 | ID),
                                         data = b_updates_sep %>% filter(trial %in% c(2:54))))
    model_kl_type    <- summary(lmerTest::lmer(scale(div) ~ kl + (1 | ID),
                                         data = b_updates_sep %>% filter(trial %in% c(2:54))))
    model_kl_int     <- summary(lmerTest::lmer(scale(div) ~ kl * trial + (1 | ID),
                                         data = b_updates_sep %>% filter(trial %in% c(2:54))))
    # Store results in a list
    regression_results[[model]] <- list(
      kl =  model_kl_type,
      tr =  model_kl_trial,
      int=  model_kl_int,
      grpa= model_kl_grp_a,
      grpb= model_kl_grp_b
    )
}

print(regression_results[[3]])
