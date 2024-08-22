%%
% Barnby, J.M. (2024) Inequality Aversion and Paranoia
%
% Joe Barnby joe.barnby@rhul.ac.uk 2022
%
%% MASTER FILE - BPD DATA sim and refit %%
%% Load data (BPDDat)

clear;
addpath('FittedFiles/CBMCode', 'Data', 'FittedFiles', 'Models', 'Models/Utilities')

%Global values
Group = ["g1", "g2"];
trials=126;

%Create data container
DFS = cell(length(Group), 1);

for j = 1:length(Group)
        %Load Dat
        Load  = readtable(strcat('Data/intent_dat_',Group(j),'.csv'));
        Loadi = table2array(Load(:,contains(Load.Properties.VariableNames,...
                {'ID', 'trial', 'S1', 'O1', 'S2', 'O2', 'choice', 'real_answer', 'Phase'})));
     
        n = size(Loadi)/trials;
        DFS{j} = cell(n(1), 1);
        for i = 1:n(1)
            DFS{j}{i} = Loadi(1+(i-1)*trials:trials+(i-1)*trials,:);
        end
end

%% Simulate the model with real parameters

% ensure the model is set to output probabilities of the model
% as well as likelihood for each phase, simulated action, and reward
% contingencies

mkdir('Data_Simulated')

Models = [
    "M1",...
    "Beta",...
    "M2",...    
    "M3",...
    "M4"
    ];

ModelsFunc = { 
    @ABA_M1,...
    @ABA_Beta,... 
    @ABA_M2,...
    @ABA_M3,... 
    @ABA_M4  
};

%%

for j=1:length(DFS)

cbmReal = load(strcat('FittedFiles/hbi_',Group(j),'_BPD.mat'));

[M,I]            = max(cbmReal.cbm.output.exceedance_prob, [], 2);
optim_parms      = cbmReal.cbm.output.parameters;

F        = zeros(length(DFS{j}),1);
results  = cell(length(DFS{j}),1);

data_sim = cell(1, length(DFS{j}));
for k=1:length(DFS{j}) 
    data_sim{k} = zeros(126, 9);
end

            for i=1:length(DFS{j})
                   
                disp(strcat('Now fitting model=', Models(I(i))))
                disp(i)
                [F(i), results{i}] = ...
                ModelsFunc{I(i)}(optim_parms{I(i)}(i, :), DFS{j}{i}, 1);
        
                data_sim{i}(:,1)=i;
                data_sim{i}(:,2)=DFS{j}{i}(:,2);
                data_sim{i}(:,3)=DFS{j}{i}(:,3);
                data_sim{i}(:,4)=DFS{j}{i}(:,4);
                data_sim{i}(:,5)=DFS{j}{i}(:,5);
                data_sim{i}(:,6)=DFS{j}{i}(:,6);
                data_sim{i}(:,7)=results{i}.simA;
                data_sim{i}(:,8)=DFS{j}{i}(:,8);
                data_sim{i}(:,9)=DFS{j}{i}(:,9);

            end

        filenom_sim_dat = strcat('Data_Simulated/data_sim_',Group(j),'.mat');
        save(filenom_sim_dat, 'data_sim')
        disp('Group (action sim) saved')

        filenom_full_sim_dat = strcat('Data_Simulated/full_sim_',Group(j),'.mat');
        save(filenom_full_sim_dat, 'F', "results")
        disp('Group (full sim) saved')

end
