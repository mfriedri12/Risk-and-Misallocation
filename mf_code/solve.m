%% solve  
%-------------------------------------------------------------------------%
% This code finds the steady state solution of a heterogenous agent model 
% with idiosyncratic risk. Alternatively, it can be used to find the global 
% solution of a representative agent model with aggregate risk. Running this 
% function with an empty argument (solve()) will use baselines for the "setup" 
% functions. If this code is adapted to a new model, the "setup" functions 
% will have to be adjusted but the rest of the code can ideally stay the same 
%-------------------------------------------------------------------------%

function [model, method, solution] = solve(model_calibration, method_calibration, initial_guesses)
tic; 

%% setup 
% setup_model defines the model, setup_method defines the particular
% solution method (algorithms, etc), and setup_solution initiatizes the objects

[model] = setup_model ; % model has 9 fields: parameters, variables, grid, functions, utilties, conditions, exogenous, clearing, and other 
if exist('model_calibration','var'); model = substitute(model, model_calibration); end % this step subs out baseline settings entered through the calibration argument 

[method] = setup_method; % method has a field corresponding to each subfunction in this folder, plus "grid" each of which has at least 3 subfields: algorithm, threholds, and iterations
if exist('method_calibration','var'); method = substitute(method, method_calibration); end 

[solution] = setup_solution(model, method); % solution has 5 fields: states, controls, value_function, distribution, prices,  
if exist('initial_guesses','var'); solution = substitute(solution, initial_guesses); end 


%% solve 
% this is a loop over prices; update controls and update distribution are
% loops themselves, while update_prices just updates price(s). 

iteration = 1; 
solution.prices.accuracy = 1e5; 
while method.update_prices.threshold < solution.prices.accuracy && iteration < method.update_prices.iterations

    fprintf('\nITERATION %d', iteration); 
    for i=1:length(model.variables.prices); fprintf(' PRICE (%s): %f', model.variables.prices(i), solution.prices.values.(model.variables.prices(i))); end 
    fprintf('\n'); 

    tic; disp('1. update solution.controls given solution.prices');       [solution] = update_controls(model, method, solution); toc; 
    tic; disp('2. update solution.distribution given solution.controls'); [solution] = update_distribution(model, method, solution); toc; 
    tic; disp('3. update solution.prices given solution.distribution');   [solution] = update_prices(model, method, solution, iteration); toc; 

    iteration = iteration + 1; 
end

fprintf('\n SOLVED!'); 

toc; 
end


%% calibration function 
% this just subs out variables in model, method, or solution, entered in calibration 

function output = substitute(default, substitute)
    output = default; 
        firstnames = fieldnames(substitute);
        for i=1:length(firstnames)
            secondnames = fieldnames(substitute.(firstnames{i}));
            for j=1:length(secondnames)
                output.(firstnames{i}).(secondnames{j}) = substitute.(firstnames{i}).(secondnames{j});
            end
        end
end
