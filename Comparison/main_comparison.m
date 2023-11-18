%% Provide larger data by panel

close all
clear all
clc

[outcome, covariate, network, union_network, population, covariate_bnj, covariate_leaders, covariate_bnj_leaders, outcome_leaders, population_leaders, covariate_woleader] = data_import_comp();

%% Environment

Q = size(network{1,1},3);
K = size(covariate{1,1},2);
N = sum(population);
G = size(outcome,2);

%% 1) Comparision with Banerjee et al.(2014) with the union network
options = optimoptions('fminunc','Display','iter-detailed','SpecifyObjectiveGradient',true,'CheckGradients',false,'Algorithm','trust-region');

% Estimating the peer effect simultaneously (SAR)

% BNJ covariates
initial = .1*ones(1,6);
f = @(B)nfxp(B, outcome, covariate_bnj, union_network, 1, G);
sar_bj = fminunc(f, initial, options);
[stars_bj, stderror_bj]=inference(sar_bj, outcome, covariate_bnj, union_network, 1, G);

% Full covariates
initial = .1*ones(1,14);
f = @(B)nfxp(B, outcome, covariate, union_network, 1, G);
sar_full_cov = fminunc(f, initial, options);
[stars_full, stderror_full]=inference(sar_full_cov, outcome, covariate, union_network, 1, G);

sar_bj = sar_bj';
sar_full_cov = sar_full_cov';

% In the following two-step estimation, the betas are estimated first only with the leaders, and
% the peer effect is estimated with the full sample.

% first-step (BNJ covariates)
initial = .1*ones(1,5);
f = @(B)baseline(B, outcome_leaders, covariate_bnj_leaders, G);
firststep_bj = fminunc(f, initial, options);
[stars_first_bj, stderror_first_bj]=inference(firststep_bj, outcome_leaders, covariate_bnj_leaders, 0, 0, G);

% second-step
initial = 0.01;
f = @(B)twostep_obj(B, firststep_bj, outcome, covariate_bnj, union_network, G);
secondstep_bj = fminunc(f, initial, options);
[stars_second_bj, stderror_second_bj]=inference_twostep(secondstep_bj, firststep_bj, outcome, covariate_bnj, union_network, G);


% first-step (full covariates)
initial = .1*ones(1,12);
f = @(B)baseline(B, outcome_leaders, covariate_leaders, G);
firststep_full = fminunc(f, initial, options);
[stars_first_full, stderror_first_full]=inference(firststep_full, outcome_leaders, covariate_leaders, 0, 0, G);

% second-step
initial = 0.01;
f = @(B)twostep_obj(B, firststep_full, outcome, covariate_woleader, union_network, G);
secondstep_full = fminunc(f, initial, options);
[stars_second_full, stderror_second_full]=inference_twostep(secondstep_full, firststep_full, outcome, covariate_woleader, union_network, G);


firststep_bj = firststep_bj';
firststep_full = firststep_full';

disp('Job done!')