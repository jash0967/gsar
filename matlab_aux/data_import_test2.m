function [outcome, covariate, network, union_network, population] = data_import_test2(weight_matrix, implicit, model)
% This function imports the microfinance data from the csv files. (Banarjee
% et al. 2014)
%
% The outcome networks can be 1) adjacency matrices and 2) negative graph
% Laplacians.
%
% The original definition of the graph Laplacian is D-A.
% To unify with the SAR model, a negative of the Laplacian, A-D is used.

disp("Importing data...")
% Set a location for the dataset
datapath = './output/';

vil_used = [1,2,3,4,6,9,12,15,19,20,21,23,24,25,29,31,32,33,36,39,42,43,45,46,47,48,50,51,52,55,57,59,60,62,64,65,67,68,70,71,72,73,75];
G = length(vil_used);

if implicit == 0
        Q = 12;
        P = Q;
else
        Q = 12 + 3;
        P = 12;
end

network=[];
union_network=[];
population=zeros(1,G);
covariate=cell(1,G);
outcome=cell(1,G);

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
    if weight_matrix == 0
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


    elseif weight_matrix == 1
        borrowlend = temp_net(:,:,1) + temp_net(:,:,6);
        borrowlend(borrowlend>0)=1;
        degree = diag(sum(borrowlend,2));
        combined(:,:,1) = borrowlend-degree;

        advice = temp_net(:,:,2) + temp_net(:,:,3);
        advice(advice>0)=1;
        degree = diag(sum(advice,2));
        combined(:,:,2) = advice-degree;

        kerocine = temp_net(:,:,4) + temp_net(:,:,5);
        kerocine(kerocine>0)=1;
        degree = diag(sum(kerocine,2));
        combined(:,:,3) = kerocine-degree;
        
        medic = temp_net(:,:,7);
        degree = diag(sum(medic,2));
        combined(:,:,4) = medic-degree;
        
        nonrel = temp_net(:,:,8);
        degree = diag(sum(nonrel,2));
        combined(:,:,5) = nonrel-degree;

        rel = temp_net(:,:,9);
        degree = diag(sum(rel,2));
        combined(:,:,6) = rel-degree;

        temple = temp_net(:,:,10);
        degree = diag(sum(temple,2));
        combined(:,:,7) = temple-degree;

        visit = temp_net(:,:,11) + temp_net(:,:,12);
        visit(visit>0)=1;
        degree = diag(sum(visit,2));
        combined(:,:,8) = visit-degree;

    end


    % Load the basic household covariates
    cov_temp = readmatrix([datapath 'vil_covariate' num2str(village_label) '.csv'], opts1);
    const_term = ones(size(cov_temp,1),1);
    covariate{ind} = [const_term, cov_temp(:, 5:end-1)];

    if implicit ~= 0
        
        for i=1:length(implicit)
            implicit_temp = zeros(population(ind),population(ind));
            implicit_cov = cov_temp(:,4+implicit(i));
            for pop1=1:population(ind)
                for pop2=1:population(ind)
                    if implicit_cov(pop1) == implicit_cov(pop2) && implicit_cov(pop1) ~= 0 && implicit_cov(pop2) ~=0
                        implicit_temp(pop1, pop2) = 1;
                    end
                end
            end

            if weight_matrix == 0
                combined(:,:,8+i) = implicit_temp - diag(diag(implicit_temp));

            elseif weight_matrix == 1
                adjacency = implicit_temp - diag(diag(implicit_temp));
                degree = diag(sum(adjacency,2));
                laplacian = adjacency-degree;
                combined(:,:,8+i) = laplacian;
            end
        end
    end


    if model == 0
        network{end+1} = combined(:,:,1:8); %#ok<*AGROW>
    elseif model == 1
        network{end+1} = combined;
    elseif model == 2
        network{end+1} = combined(:,:,8+1:end);
    end

    union_network{end+1}=union_temp;

    % Load the Take-up data
    MF_takeup = load([datapath 'MF_outcome' num2str(village_label) '.csv']); 

    outcome{ind} = MF_takeup;

    
end
disp("Import successful")
end

