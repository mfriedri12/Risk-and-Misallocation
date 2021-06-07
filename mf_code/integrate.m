%% integrate
%-------------------------------------------------------------------------%
% Calculate expected value as a function of the controls that are also states 
% current exogenous state; integrate out exogenous states. In other words, takes a guess
% for the next period value function, solution.values, (a function of states one period ahead)
% and then calculates the expected value function, returning a function of 
% endogenous states one period ahead (controls) and the current exogenous state. 
% Can do this either by simply premultiplying the value function guess by some
% exogenous state transistion matrix (ie an AR(1) approximagted by rowenhurst or tauchen
% method), or by taking a value function approximated by a basis and doing 
% gaussian quadrature. Just does one expected variable at the moment 
% NOTE this is hard for just one exogenous shock right now but could be adapted 
%-------------------------------------------------------------------------%

% I need to rewrite this so it returns both a function and values 

function [expected_values] = integrate(model, method, solution, basis, output)

%% get expected values  

switch method.interpolate.algorithm 
    case 'griddedInterpolant'
        expected_values = basis(solution.states.cell_prime);
        expected_values = transpose(solution.states.weights.(model.variables.exogenous(1)))*reshape(expected_values, method.integrate.nodes, prod(solution.states.dimensions.states)); 
        expected_values = reshape(expected_values,solution.states.dimensions.states);
    case 'basis'
end

%% get a basis if necessary 

switch output
    case 'basis'
        expected_values = interpolate(method, solution, expected_values); %not correct
    case 'values'
        expected_values = expected_values;
end


end