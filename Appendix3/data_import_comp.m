function [outcome, covariate, network, union_network, population, covariate_bnj, covariate_leaders, covariate_bnj_leaders, outcome_leaders, population_leaders, covariate_woleader] = data_import_comp()
% This function imports the microfinance data from the csv files. (Banarjee
% et al. 2014)
%
% The outcome networks can be 1) adjacency matrices and 2) negative graph
% Laplacians.
%
% The original definition of the graph Laplacian is D-A.
% To unify with the SAR model, a negative of the Laplacian, A-D is used.


% Set a location for the dataset
datapath = 'C:\Users\jkshi\Documents\py_projects\DataBuilder_hh\output\';

vil_used = [1,2,3,4,6,9,12,15,19,20,21,23,24,25,29,31,32,33,36,39,42,43,45,46,47,48,50,51,52,55,57,59,60,62,64,65,67,68,70,71,72,73,75];
G = length(vil_used);
Q = 12;

network=[];
union_network=[];
population=zeros(1,G);
covariate=cell(1,G);
outcome=cell(1,G);
covariate_bnj=cell(1,G);

population_leaders=zeros(1,G);
covariate_leaders=cell(1,G);
outcome_leaders=cell(1,G);
covariate_bnj_leaders=cell(1,G);

covariate_woleader=cell(1,G);

opts1 = detectImportOptions([datapath 'vil_covariate29.csv']);

for ind=1:G
    % Assign the label of the village used
    village_label = vil_used(ind);
    % disp(['Importing village ', num2str(village_label), '...'])

    % Load a village folder with csv files
    myfiles = dir([datapath num2str(village_label)]);
    filenames={myfiles(:).name};
    csvfiles=filenames(endsWith(filenames,'.csv'));

    % Using the first network, find its population and record it
    temp = load([datapath num2str(village_label) '\' csvfiles{1}]);
    population(ind) = size(temp,1);

    % Assign a template matrix
    temp_net = zeros(population(ind),population(ind),Q);

    % Convert each csv file to a specified weight matrix

    for i = 1:12
        adjacency = load([datapath num2str(village_label) '\' csvfiles{i}]);
        temp_net(:,:,i) = adjacency + adjacency';
        temp_net(temp_net>0)=1;
    end

    combined = zeros(population(ind),population(ind),8);
    union_temp = zeros(population(ind),population(ind));

    borrowlend = temp_net(:,:,1) + temp_net(:,:,6);
    borrowlend(borrowlend>0)=1;
    combined(:,:,1) = borrowlend;

    advice = temp_net(:,:,2) + temp_net(:,:,3);
    advice(advice>0)=1;
    combined(:,:,2) = advice;

    kerocine = temp_net(:,:,4) + temp_net(:,:,5);
    kerocine(kerocine>0)=1;
    combined(:,:,3) = kerocine;

    combined(:,:,4) = temp_net(:,:,7);

    combined(:,:,5) = temp_net(:,:,8);

    combined(:,:,6) = temp_net(:,:,9);

    combined(:,:,7) = temp_net(:,:,10);

    visit = temp_net(:,:,11) + temp_net(:,:,12);
    visit(visit>0)=1;
    combined(:,:,8) = visit;

    for q=1:8
        union_temp = union_temp + combined(:,:,q);
    end
    union_temp(union_temp~=0)=1;


    % Load the Take-up data
    MF_takeup = load([datapath 'MF_outcome' num2str(village_label) '.csv']);

    outcome{ind} = MF_takeup;

    % Load the basic household covariates
    cov_import = readmatrix([datapath 'vil_covariate' num2str(village_label) '.csv'], opts1);
    const_term = ones(size(cov_import,1),1);
    leaders = cov_import(:,8);

    covariate{ind} = [const_term, cov_import(:, 5:end-1)];
    covariate_bnj{ind} = [const_term, cov_import(:, 6:7), cov_import(:,14), cov_import(:,16)]; % using only the covariates used in Banerjee et al.(2014)

    cov_temp = covariate{ind};
    cov_bnj_temp = covariate_bnj{ind};
    outcome_temp = outcome{ind};

    
    population_leaders(ind) = nnz(leaders);

    cov_leaders_temp = zeros(population_leaders(ind), size(cov_temp, 2));
    cov_bnj_leaders_temp = zeros(population_leaders(ind), size(cov_bnj_temp,2));
    outcome_leaders_temp = zeros(population_leaders(ind), 1);

    % Pick leaders only
    j=1;
    for i=1:population(ind)
        if leaders(i)==1
            cov_leaders_temp(j,:) = cov_temp(i,:);
            cov_bnj_leaders_temp(j,:) = cov_bnj_temp(i,:);

            outcome_leaders_temp(j,:) = outcome_temp(i,:);
            j=j+1;
        end
    end

    cov_leaders_temp(:,5)=[]; % the leader dummy variable is redundant
    
    covariate_leaders{ind} = cov_leaders_temp;
    covariate_bnj_leaders{ind} = cov_bnj_leaders_temp;
    cov_temp(:,5)=[];
    
    covariate_woleader{ind} = cov_temp;

    outcome_leaders{ind} = outcome_leaders_temp;

    network{end+1} = combined;
    union_network{end+1}=union_temp;

end
disp('Import successful')
end

