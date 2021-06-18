%% update controls 
%-------------------------------------------------------------------------%
% given prices, solve single agent problem. all algorithms follow these steps 
% 0. guess the VALUES of the value or the control (policy) function on the state space grid (first guess made in "setup_solution")
% 1. usually, *interpolate* the above to get a guess for the value or the control (policy) FUNCTION. (Can avoid interpolaiton if you integrate out exogenous states via a markov approximation of the exogenous state; butthis solution will also usually be less accurate)
% 2. update your VALUES using the bellman equation or some deriviative of it (FOC, EC, EE) that is also a contraction
%      2.1 usually you'll need to *integrate* over exogenous states to numerically calculate the expected value  
%      2.2 sometimes you'll need to nonlinearly *maximize* to get the updated VALUES, using a hill-climbing algorithm
%      2.2 other tiems you'll need to take a derivative and *root-find*
%      2.3 usually you'll need to enforce borrowing constraints 
% 3. check for convergence, that the distance between your guess & update is less than thresholds chosen in setup_method
% 4. if not converged, set guess = update and repeat until convergence 

% *interpolate* takes in values guesses and returns an object with 4 fields: values, basis, coefficients, and nodes, all defined over the state space   
% *integrate* takes in an object returned by interpolate and returns expected values integrating over exogenous states  
% *solve* takes in a function and returns its roots OR its maximization  

%-------------------------------------------------------------------------%

function [solution] = update_controls(model, method, solution)


%% reset 
if method.update_controls.reset 
    for i=1:length(model.variables.controls)
        solution.controls.values.(model.variables.controls(i)) = ones(solution.states.dimensions.states);  % guess controls as 1 
    end
    solution.value_function.values =  ones(solution.states.dimensions.states) ; % guess valus as 1 
end


%% solve for controls that are endogenous states by following steps 1-4 
iteration = 1; 
solution.controls.accuracy = 1e5;
last = solution.controls.accuracy; 
while method.update_controls.threshold < solution.controls.accuracy && iteration < method.update_controls.iterations     
    
    switch method.update_controls.algorithm
        
        case 'VFI'
            if method.update_controls.vtrick
                expected_value_function = model.utilities.vtrick_inv(solution.value_function.values, model.parameters);
            else 
                expected_value_function = solution.value_function.values;
            end
            expected_value_function = interpolate(method, solution, expected_value_function); % 1 
            expected_value_function = integrate(model, method, solution, expected_value_function, 'basis'); % 2.1           
            if length(model.variables.endogenous)==1 % no portfolio problem defining bellman on state space grid
                if method.update_controls.vtrick
                    bellman = @(c1, solution) model.conditions.b1(c1, solution.states.ndgrid, solution.prices.values, model.parameters) + model.parameters.beta*model.utilities.vtrick(expected_value_function(solution.states.ndgrid.z, c1), model.parameters); % 2.2  
                else 
                    bellman = @(c1, solution) model.conditions.b1(c1, solution.states.ndgrid, solution.prices.values, model.parameters) + model.parameters.beta*expected_value_function(solution.states.ndgrid.z, c1) ; % 2.2  
                end
            elseif length(model.variables.endogenous)==2 %portfolio problem defining bellman on state space plus control grid
                if method.update_controls.vtrick
                    bellman = @(c1, solution) model.conditions.b1(c1, solution.states.ndgridplus.(strcat(model.variables.endogenous(2),model.variables.endogenous(2))), solution.states.ndgridplus, solution.prices.values, model.parameters) + model.parameters.beta*model.utilities.vtrick(expected_value_function(solution.states.ndgridplus.z, c1, solution.states.ndgridplus.(strcat(model.variables.endogenous(2),model.variables.endogenous(2)))), model.parameters); % 2.2  
                else 
                    bellman = @(c1, solution) model.conditions.b1(c1, solution.states.ndgridplus.(strcat(model.variables.endogenous(2),model.variables.endogenous(2))), solution.states.ndgridplus, solution.prices.values, model.parameters) + model.parameters.beta*                      (expected_value_function(solution.states.ndgridplus.z, c1, solution.states.ndgridplus.(strcat(model.variables.endogenous(2),model.variables.endogenous(2))))); % 2.2  
                end
            end 
            search_max = model.utilities.max.(model.variables.endogenous(1))(solution, model.parameters); % maximium of search space, must not be so big that c is negative while also respecting borrowing constraint
            search_min = model.grid.min(1)*ones(size(search_max)); % minimum of search space, borrowing constraint
            search_max = max(search_max,search_min);
            [solution.controls.values.(model.variables.endogenous(1)), update] = maximize(method, bellman, search_min, search_max, solution); % 2.2
            if length(model.variables.endogenous)==2
                solution.controls.basis.(model.variables.endogenous(1)) = interpolate(method, solution, solution.controls.values.(model.variables.endogenous(1)), solution.states.cell_plus); 
                solution.controls.basis.(model.variables.endogenous(1)) = @(c2,states) solution.controls.basis.(model.variables.endogenous(1))(states.(model.variables.states(1)),states.(model.variables.states(2)),states.(model.variables.states(3)),c2); % hard coded for 3 states right now 
                if method.update_controls.vtrick
                    bellman = @(c2, solution) model.conditions.b1(solution.controls.basis.(model.variables.endogenous(1))(c2,solution.states.ndgrid), c2, solution.states.ndgrid, solution.prices.values, model.parameters) + model.parameters.beta*model.utilities.vtrick(expected_value_function(solution.states.ndgrid.z, solution.controls.basis.(model.variables.endogenous(1))(c2, solution.states.ndgrid), c2), model.parameters); % 2.2  
                else 
                    bellman = @(c2, solution) model.conditions.b1(solution.controls.basis.(model.variables.endogenous(1))(c2,solution.states.ndgrid), c2, solution.states.ndgrid, solution.prices.values, model.parameters) + model.parameters.beta*expected_value_function(solution.states.ndgrid.z, solution.controls.basis.(model.variables.endogenous(1))(c2, solution.states.ndgrid), c2); % 2.2  
                end
                search_max = model.utilities.max.(model.variables.endogenous(2))(solution, model.parameters); % maximium of search space, must not be so big that c is negative while also respecting borrowing constraint
                search_min = model.grid.min(2)*ones(size(search_max)); % minimum of search space, borrowing constraint
                search_max = max(search_max,search_min);
                [solution.controls.values.(model.variables.endogenous(2)), update] = maximize(method, bellman, search_min, search_max, solution); % 2.2    
                solution.controls.values.(model.variables.endogenous(1)) = solution.controls.basis.(model.variables.endogenous(1))(solution.controls.values.(model.variables.endogenous(2)), solution.states.ndgrid);   
            end
            solution.controls.accuracy = max(abs(update - solution.value_function.values), [], 'all'); % 3
            solution.value_function.values = update; % 4 
        
          
    end

    if mod(iteration,method.update_controls.print)==0
        fprintf('iteration: %i, precision: %i \n', iteration, solution.controls.accuracy);
        if last*10 < solution.controls.accuracy; fprintf('WARNING! algorithm exploding at iteration: %i. Paused. \n',iteration); pause; end
        last = solution.controls.accuracy;
    end
    iteration = iteration + 1; 
end


%% get "other" controls 
% this is usually consumption 

for i=1:length(model.variables.other)
    solution.controls.values.(model.variables.other(i)) = model.other.(model.variables.other(i))(solution, model.parameters); 
end


%% get control funciton basis 
% sometimes this is incidental in the algorithm but not always 
for i=1:length(model.variables.controls)
    solution.controls.basis.(model.variables.controls(i)) = interpolate(method, solution, solution.controls.values.(model.variables.controls(i)));
end

 
end
        