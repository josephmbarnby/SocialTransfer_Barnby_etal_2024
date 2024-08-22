%%
% Barnby, J.M. (2024) Inequality Aversion and Paranoia
%
% Joe Barnby joe.barnby@rhul.ac.uk 2022
%
%% MASTER FILE %%
%% Load data (BPDDat)

addpath('CBMCode', 'Data', 'Models', 'Models/Utilities')

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

%% Set WDs

mkdir('FittedFiles')

Models = [
    "M1",...
    "Beta",...
    "M2",...
    "M3",...
    "M4"...
    ];

Folder = 'FittedFiles/Lap_Subj_BPD';

%BB

for j = 1:length(Group)
        for q = 1:length(Models)

        S = Group(j); M = Models(q);
        mkdir(fullfile(Folder, strcat('lap_Subj_',M,'_',S)));
        end
end

%% Set Priors

v = 6.5;

ModelsFunc = {    
@ABA_M1,...  
@ABA_Beta,... 
@ABA_M2,...
@ABA_M3,... 
@ABA_M4
};

prior = {...
struct('mean', zeros(6,1), 'variance', v),...
struct('mean', zeros(3,1), 'variance', v),...
struct('mean', zeros(6,1), 'variance', v),...
struct('mean', zeros(8,1), 'variance', v),...
struct('mean', zeros(8,1), 'variance', v)
};

%% Run Lap
%% BB Models (Extended)

mkdir('FittedFiles/LaplaceFit_BPD')
MetaFol    = 'FittedFiles/Lap_Subj_BPD';

for j = 1:length(Group) % run a loop for each subsample
        for q = 1:length(Models)

        DatUsei = DFS{j}; 
        S = Group(j); 
        M = Models(q);
        Folder = strcat('lap_Subj_',M,'_',S); 
        Subj = strcat('lap_Subj_',M,'_',S, '_'); 
        file = strcat('FittedFiles/LaplaceFit_BPD/lap_',M,'_',S,'.mat');

        parfor i = 1:length(DatUsei) %nested parfor loop for fitting
        %n=fitme(i);
        % 1st input: data
        data_subj = DatUsei(i);
        % 2nd input: function handle of model (i.e. @model_mf)
        % 3rd input: a prior struct.
        % 4th input: output file
        fname_subj = fullfile(MetaFol, Folder,strcat(Subj, num2str(n), '.mat'));
        cbm_lap(data_subj, ModelsFunc{q}, prior{q}, fname_subj);
        end
       
        CALC = cell(length(DatUsei),1);
        for n=1:length(CALC)
            CALC{n} = fullfile(MetaFol, Folder,strcat(Subj, num2str(n), '.mat'));
        end

        CALCBIND = CALC;
       
        fname_BF = file;
        cbm_lap_aggregate(CALCBIND,fname_BF);
 
        end
end

%% HBM CBM for bayes run g1

fcbm_maps = {'FittedFiles/LaplaceFit_BPD/lap_M1_g1.mat',...
             'FittedFiles/LaplaceFit_BPD/lap_Beta_g1.mat',...
             'FittedFiles/LaplaceFit_BPD/lap_M2_g1.mat',...
             'FittedFiles/LaplaceFit_BPD/lap_M3_g1.mat',...
             'FittedFiles/LaplaceFit_BPD/lap_M4_g1.mat'};

BFS_hbi = 'FittedFiles/hbi_g1_BPD.mat';

cbm_hbi(DFS{1},ModelsFunc,fcbm_maps,BFS_hbi);

%% HBM CBM for bayes run g2

fcbm_maps = {'FittedFiles/LaplaceFit_BPD/lap_M1_g2.mat',...
             'FittedFiles/LaplaceFit_BPD/lap_Beta_g2.mat',...
             'FittedFiles/LaplaceFit_BPD/lap_M2_g2.mat',...
             'FittedFiles/LaplaceFit_BPD/lap_M3_g2.mat',...
             'FittedFiles/LaplaceFit_BPD/lap_M4_g2.mat'};

BFS_hbi = 'FittedFiles/hbi_g2_BPD_UNCONLY.mat';

cbm_hbi(DFS{2},ModelsFunc,fcbm_maps,BFS_hbi);
