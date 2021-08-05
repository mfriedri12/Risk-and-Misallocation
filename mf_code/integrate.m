%% integrate
%-------------------------------------------------------------------------%
% Calculate expected value as a function of the controls that are also states 
% and the CURRENT exogenous state; integrate out exogenous states.
%-------------------------------------------------------------------------%

function [expected_values] = integrate(model, method, solution, func)

%% get expected values  
expected_values = func(solution.states.stack.prime);  % takes stack       
expected_values = transpose(solution.states.weights)*reshape(expected_values, method.integrate.nodes^(length(model.variables.exogenous)), prod(solution.states.dimensions.states)); 
expected_values = reshape(expected_values,solution.states.dimensions.states);


end

