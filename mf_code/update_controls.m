%% update controls 
%-------------------------------------------------------------------------%
% given prices, solve single agent problem and get control functions, which is what you need for update_distribution 
% *interpolate* takes in values guesses and returns an interpolant 
% *integrate* takes in an object returned by interpolate and returns expected values integrating over exogenous states  
% *solve* takes in a function and returns its roots OR its maximization  
%-------------------------------------------------------------------------%

function [solution] = update_controls(model, method, solution)

%% solve single agent probem 
% variety of algorithms to do this
iteration = 1; 
solution.controls.accuracy = 1e15;
while method.update_controls.threshold < solution.controls.accuracy && iteration < method.update_controls.iterations     
    
  switch method.update_controls.algorithm      
    
    case 'VFI'
            
    % 1. get expected value function with value function guess
    if method.update_controls.vtrick;  solution.value_function.values = model.utilities.vtrick_inv(solution.value_function.values, model.parameters); end
    solution.value_function.function = interpolate(model, method, solution, solution.value_function.values, 'value_function'); % 1
    expected_value_function = integrate(  model, method, solution, solution.value_function.function); % 2.1 
    expected_value_function = interpolate(model, method, solution, expected_value_function,        'value_function');
    if method.update_controls.vtrick; expected_value_function = @(states) model.utilities.vtrick(expected_value_function(states), model.parameters); end

    % 2. maximize to update 
    if length(model.variables.endogenous)==1 % no portfolio problem defining bellman on state space grid
      bellman = @(c1,solution) model.conditions.b1(c1, solution.states.stack, solution.prices.values, model.parameters)...
                             + model.parameters.beta*expected_value_function([solution.states.stack.(model.variables.states)(:,1:length(model.variables.exogenous)) c1]) ; % 2.2  
      search_max = model.utilities.max.(model.variables.endogenous(1))(solution.states.stack, solution.prices.values, model.parameters); % maximium of search space, must not be so big that c is negative while also respecting borrowing constraint
    elseif length(model.variables.endogenous)==2 %portfolio problem defining bellman on state space plus control grid
      bellman = @(c1,solution) model.conditions.b1(c1, solution.states.stackplus.(strcat(model.variables.endogenous(2),model.variables.endogenous(2))), solution.states.stackplus, solution.prices.values, model.parameters) ...
                             + model.parameters.beta*expected_value_function([solution.states.stackplus.(model.variables.states)(:,1:length(model.variables.exogenous)) c1 solution.states.stackplus.(strcat(model.variables.endogenous(2),model.variables.endogenous(2)))]); 
      search_max = model.utilities.max.(model.variables.endogenous(1))(solution.states.stackplus, solution.prices.values, model.parameters); % maximium of search space, must not be so big that c is negative while also respecting borrowing constraint     
    end
    if method.update_controls.enforce_grid_max; search_max = min(search_max, model.grid.max.(model.variables.endogenous(1))); end
    search_min = model.grid.min.(model.variables.endogenous(1))*ones(size(search_max));
    [solution.controls.values.(model.variables.endogenous(1)), update] = maximize(method, bellman, search_min, search_max, solution); % 2.2

    % 2. second maximize loop if portfolio problem 
    if length(model.variables.endogenous)==2  
      solution.controls.function.(model.variables.endogenous(1)) = interpolate(model, method, solution, solution.controls.values.(model.variables.endogenous(1)), 'portfolio_problem'); 
      bellman = @(c2,solution) model.conditions.b1(solution.controls.function.(model.variables.endogenous(1))([solution.states.stack.(model.variables.states) c2]), c2, solution.states.stack, solution.prices.values, model.parameters)  ...
                             + model.parameters.beta*expected_value_function([solution.states.stack.(model.variables.states)(:,1:length(model.variables.exogenous))  solution.controls.function.(model.variables.endogenous(1))([solution.states.stack.(model.variables.states) c2]) c2]) ;           % 2.2    
      search_max = model.utilities.max.(model.variables.endogenous(2))(solution.states.stack, solution.prices.values, model.parameters); % maximium of search space, must not be so big that c is negative while also respecting borrowing constraint
      if method.update_controls.enforce_grid_max; search_max = min(search_max, model.grid.max.(model.variables.endogenous(2))); end
      search_min = model.grid.min.(model.variables.endogenous(2))*ones(size(search_max)); % minimum of search space, borrowing constraint
      [solution.controls.values.(model.variables.endogenous(2)), update] = maximize(method, bellman, search_min, search_max, solution); % 2.2    
      solution.controls.values.(model.variables.endogenous(1)) = solution.controls.function.(model.variables.endogenous(1))([solution.states.stack.(model.variables.states) solution.controls.values.(model.variables.endogenous(2))]);   
     end
    
    % 3. Check for convergence & update
    accuracy = max(abs(update - solution.value_function.values), [], 'all'); % 3
    solution.value_function.values = update;
  
  if mod(iteration,method.update_controls.print)==0        
    fprintf('iteration: %i, precision: %i \n', iteration, accuracy);
    if ~method.update_controls.enforce_grid_max 
      for i=1:length(model.variables.endogenous)
        fprintf('percent of control %s above top of grid: %f \n',model.variables.endogenous(i), 100*(sum(solution.controls.values.(model.variables.endogenous(i))>max(solution.states.grid.(model.variables.endogenous(i))),'all')/length(solution.controls.values.(model.variables.endogenous(i)))))
      end
    end
    if (accuracy - solution.controls.accuracy) > 5e-1; fprintf('WARNING! algorithm exploding at iteration: %i. Paused. \n',iteration); pause; end
    if abs(accuracy - solution.controls.accuracy) < 1e-8; fprintf('Convergence plateauing, halting iteration at %i. Paused. \n',iteration); accuracy = method.update_controls.threshold; end      
  end
  solution.controls.accuracy = accuracy;
  iteration = iteration + 1;
end


%% get "other" control values
% this is usually consumption 
for i=1:length(model.variables.other)
  solution.controls.values.(model.variables.other(i)) = model.other.(model.variables.other(i))(solution.controls.values, solution.states.stack, solution.prices.values, model.parameters); 
end


%% get control functions  
% sometimes this is incidental in the algorithm but not always 
for i=1:length(model.variables.controls)
  solution.controls.function.(model.variables.controls(i)) = interpolate(model, method, solution, solution.controls.values.(model.variables.controls(i)),'value_function');
end

 
end        