clear all
clc
closed = true;
medium = "High" ;%Low High

[num,txt,~] = xlsread("../Data/gene_fpkm.csv");
genes_in_dataset = txt(2:end,1);
conditions = txt(1,2:end);
gene_exp = num(:, 1:end);

load('..\Model\tinit_model.mat')
init_model.c(contains(init_model.rxns,'MAR13082')) = 1;
model = init_model;

%load('..\Model\Human-GEM.mat')
%model = ihuman;
model.c(contains(model.rxns,'MAR13082')) = 1;
%[model, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

% set medium constrain
if medium =="High"
    [model,basisMedium]  = change_medium_high_glucose(model);
else
    [model,basisMedium]  = change_medium_low_glucose(model);
end

model = setParam(model, 'obj', 'MAR03964', 1);
FBAsolution = solveLP(model);
max_atp = FBAsolution.f;
disp(max_atp)
% define parameters over which to iterate

gamma = [2];
threshold = [25];

[reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(model);

genes = model.genes;
Cmin{size(gene_exp, 2), length(model.rxns)} = [];
Cmax{size(gene_exp, 2), length(model.rxns)} = [];
GeneExpressionArray = ones(numel(genes),1); 

k = 1;
tic

for g = 1:numel(gamma)                
    changeCobraSolver('gurobi', 'all');
    %changeCobraSolverParams('QP', 'feasTol', 1e-3);
    %changeCobraSolverParams('LP', 'feasTol', 1e-3);
    %changeCobraSolverParams('QP', 'method', 1);
    %changeCobraSolver('ibm_cplex', 'all');
    
	gam = gamma(g);
    new_k = main(medium,threshold, gam, gene_exp,genes_in_dataset,conditions,model,genes,GeneExpressionArray,g,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length,k);
    k = new_k;
end
toc

save(['Results\CRC'])

function new_k = main(medium,threshold, gam, gene_exp,genes_in_dataset,conditions,t_model,genes,GeneExpressionArray,g,reaction_expression,pos_genes_in_react_expr,ixs_geni_sorted_by_length,k)
    gamma = gam;
    if isempty(gcp('nocreate'))
        parpool();   
    end
    
    pool = gcp('nocreate'); % check if pool was successfully created
    if isempty(pool) % if there is no active pool, then throw an error
        error('\nError! No parallel pool created!')
    end
    for tr = 1:numel(threshold)                                                       
        % cut data with respect to a threshold  
        %data = gene_exp./mean(gene_exp);                  % get the fold change (If A is a matrix, then mean(A) returns a row vector containing the mean of each column.)
        % data = data.*(data>prctile(data,prc_cut,1));    % if you want to binaruse data according to a percentile                                        
        data = gene_exp;
       % applying the bounds
        fprintf('Iteration (k): %d, Gamma: %d, Threshold: %d\n',k,gamma(g),threshold(tr));
        for t=1:size(data,2)          % in here we select a unique profile(one patient at the time)
           	disp(t)
            model = t_model;
            
            %set ATP as objective 
            model = setParam(model, 'obj', 'MAR03964', 1);
            FBAsolution = solveLP(model);
            max_atp = FBAsolution.f;
            disp(max_atp)
            if medium =="High"
                disp("High glucose medium")
                if contains(conditions(t), 'm4')
                % set ATP to 155.12%
                    disp(1.55*max_atp)
                    model = changeRxnBounds(model, 'MAR03964', 1.55*max_atp, 'l');
                else 
                % set ATP to 100%
                    model = changeRxnBounds(model, 'MAR03964', max_atp, 'l');
                end

            else
                disp("Low glucose medium")
                if contains(conditions(t), 'm4')
                % set ATP to 133.7%
                    disp(1.34*max_atp)
                    model = changeRxnBounds(model, 'MAR03964', 1.34*max_atp, 'l');
                else 
                % set ATP to 100%
                    model = changeRxnBounds(model, 'MAR03964', max_atp, 'l');
                end
            end
            

            % set biomass as objective function
            model = setParam(model, 'obj', 'MAR13082', 1);

            if mod(t,3) ==0
                delete(pool)
                if isempty(gcp('nocreate'))
                    parpool();   
                end
                pool = gcp('nocreate'); % check if pool was successfully created
                if isempty(pool) % if there is no active pool, then throw an error
                    error('\nError! No parallel pool created!')
                end
            end
            expr_profile = data(:,t);
            pos_genes_in_dataset = zeros(numel(genes),1);% gene in the model human 1
            for i=1:length(pos_genes_in_dataset)
                position = find(strcmp(genes{i},genes_in_dataset),1); 
                if ~isempty(position)                                   
                    pos_genes_in_dataset(i) = position(1);              
                    GeneExpressionArray(i) = expr_profile(pos_genes_in_dataset(i));         
                end
            end
            if or(sum(isinf(GeneExpressionArray)) >= 1, sum(isnan(GeneExpressionArray)) >= 1)
                fprintf('\nError in the gene expression data!');
            end
            [minfluxes, maxfluxes] = transcriptomic_bounds(conditions(t), gamma(g), GeneExpressionArray, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length);
            Cmin(k,:) = [conditions(k), num2cell(transpose(minfluxes))];
            Cmax(k,:) = [conditions(k), num2cell(transpose(maxfluxes))];
            k = k +1;
        end
    end
    new_k = k; 
    Tmin = cell2table(Cmin);
    Tmin.Properties.VariableNames = vertcat({'conditions'}, model.rxns);

    Tmax = cell2table(Cmax);
    Tmax.Properties.VariableNames = vertcat({'conditions'}, model.rxns);
    %T.Properties.VariableNames{1} = 'patient_id';  
    writetable(Tmin,'Results\tinit_FVA_min_high_glucose.csv');     
    writetable(Tmax,'Results\tinit_FVA_max_high_glucose.csv');     % set the folder where you want the data to be saved
end


