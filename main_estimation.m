%% Provide larger data by panel

close all
clear all
clc

addpath("matlab_aux/")

% Import dataset
% data_import(0) = SAR
% data_import(1) = LAR

implicit = [11 12 5];

% options
% SAR = 0, LAR = 1
% implicit networks specification
% Explicit only = 0, nested = 1, implicit only = 2
% Inversed outcome: no = 0, yes = 1
% Symmetrize the surveyed networks: no = 0, yes = 1

[outcome, covariate, network, union_network, population] = data_import_test2(0,implicit,0);

%% Group fixed effects

G = size(outcome,2);

for i=1:G
    temp=covariate{i};
    n=size(temp,1);
    group_fixed=zeros(n, G-1);
    if i~=1
        group_fixed(:,i-1)=1;
    end
    temp=[temp group_fixed];
    covariate{i}=temp;
end

clear temp

%% Environment

Q = size(network{1,1},3);
K = size(covariate{1,1},2);
G = size(outcome,2);

N = sum(population);

%% Finding the spectral radius

[radius_bound, radius] = RadiusFinder(network, G, Q);
esti_bound = 1/(radius_bound * 0.25);

%% Estimation

options = optimoptions('fminunc','Display','iter-detailed','SpecifyObjectiveGradient',true,'CheckGradients',false,'Algorithm','trust-region');

% 0) Baseline Estimation
disp("Estimating the benchmark logit model (Model 1)")

initial = .1*ones(1,K);
f = @(B)baseline(B, outcome, covariate, G);
wo_network = fminunc(f, initial, options);
wo_network = wo_network';

% Inference
[stars_vanilla, stderror_vanilla]=inference(wo_network', outcome, covariate, 0, 0, G);

% 1) Union network

disp("Estimating the SAR model with the union network (Model 2)")

initial = .1*ones(1,1+K);
initial = [0 wo_network'];
f = @(B)nfxp(B, outcome, covariate, union_network, 1, G);
NFXP_hat_union = fminunc(f, initial, options);
NFXP_hat_union = NFXP_hat_union';

spectral_rad_u = sr(NFXP_hat_union', union_network, covariate, 1, G);
smtp_u = smtp(NFXP_hat_union', outcome, covariate, union_network, 1, G);

% 1) Higher-order SAR

disp("Estimating the higher-order GSAR model (Model 5)")

initial = .01*ones(1,Q+K);
f = @(B)nfxp_flex(B, outcome, covariate, network, Q, G);
temp = zeros(1, Q) +.001;
initial = [temp 0 wo_network'];
% initial = [temp 0 zeros(1,K)];

NFXP_hat = fminunc(f, initial, options);
NFXP_hat=NFXP_hat'; 


% Inference

[stars, stderror]=inference_flex(NFXP_hat', outcome, covariate, network, Q, G);

%% Spectral Radius Condition

[spectral_rad, outcome_network] = sr_flex(NFXP_hat', network, covariate, Q, G);

if nnz(find(abs(spectral_rad)>4)) ~= 0
    error('Error: spectral radius condition violated!')
end

%% Summary statistics

% Social multiplier
smtp=smtp_flex(NFXP_hat', outcome, covariate, network, Q, G);

% Marginal effect
marginal_vanilla=marginal_effect_vanilla(wo_network',covariate,G)';
[direct_union, indirect_union, naive_union] = marginal_effect(NFXP_hat_union',outcome, covariate, union_network, 1, G);
[direct, indirect, naive] = marginal_effect_flex(NFXP_hat', outcome, covariate, network, Q, G);


%% Print Results

results=cell2table(cell(0,64));
results=renamevars(results,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8"],["Money","Advice","Kerosine","Medical","Non-relatives","Relatives","Temple","Visit"]);
results=renamevars(results,"Var9","conf. param.");
results=renamevars(results,["Var10","Var11","Var12","Var13","Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21","Var22"],["constant","roof: tile", "room no.", "bed no.", "leader", "caste_gen", "caste_min", "caste_sc", "caste_st", "islam", "lat_none", "house_notowned", "noelectricity"]);
writetable(results,'results.csv','Delimiter','comma')

writematrix(NFXP_hat','results.csv','Delimiter','comma','WriteMode','append');
writematrix(stars','results.csv','Delimiter','comma','WriteMode','append');
writematrix(stderror','results.csv','Delimiter','comma','WriteMode','append');

disp('Job done!')