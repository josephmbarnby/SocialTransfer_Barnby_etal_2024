# SocialTransfer (Barnby et al., 2024)
Repo for the paper: "Self-Other Generalisation Shapes Social Interaction and Is Disrupted in Borderline Personality Disorder"

* `Data` <br />
#real behavioural data <br />
     ├── intent_dat.csv # cleaned behavioural data from the Intentions Game <br />
* `Data_Simulated` <br />
#model-simulated behavioural data for each group {g1=BPD, g2=CON} <br />
     ├── full_sim_g1.csv #full model results (M4) <br />
     ├── full_sim_g2.csv #full model results (M1) <br />
     ├── data_sim_g1.csv #model-simulated actions (M4) <br />
     ├── data_sim_g2.csv #model-simulated actions (M1) <br />
     ├── full_sim_g1_M3.csv #full model results for M3 <br />
     ├── full_sim_g2_M3.csv #full model results for M3 <br />
* `Models` <br />
#Model files to fit and simulate data.  <br />
     ├── M1.m <br />
     ├── M2.m <br />
     ├── M3.m <br />
     ├── M4.m <br />
     ├── Beta.m <br />
     ├── Utilities <br />
* `Analysis` <br />
# Run this to reproduce all results from the paper <br />
     ├── BPD_CON_Analysis.R # Core analysis file to produce outcomes in the paper <br/>
