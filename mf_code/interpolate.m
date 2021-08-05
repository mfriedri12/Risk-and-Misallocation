%% interpolate 
%-------------------------------------------------------------------------%
% interpolate a state-space sized guess either for a value or a control (policy) function. 
% can do this either with a basis or with griddedInterpolent; if basis has optimal nodes, update the state space to 
% correspond to those optimal nodes; this is why this function returns "solution" as well as guess.  
% griddedInterpolant came first so my code is kind of subordinate to it's weird arguments 
%-------------------------------------------------------------------------%

function func = interpolate(model, method, solution, values, type) 

switch type 
    
  case 'value_function'    
    switch method.interpolate.toolkit.value_function  
      case 'griddedInterpolant'
        func = griddedInterpolant(solution.states.cell, reshape(values,solution.states.dimensions.states), method.interpolate.basis.value_function,method.interpolate.extrapolate.value_function); %takes cells, ndgrids, OR stacks 
      case 'compecon'   
        coefficients =     funbas(solution.value_function.basis, solution.states.stack.(model.variables.states), zeros(1,length(model.variables.states)))\reshape(values,[],1);
        func = @(states)   funbas(solution.value_function.basis, states,                                         zeros(1,length(model.variables.states)))*coefficients; %takes stacks
  
    end
  
 case 'portfolio_problem'  
    switch method.interpolate.toolkit.portfolio_problem   
      case 'griddedInterpolant'
        func = griddedInterpolant(solution.states.cell_plus, reshape(values,[solution.states.dimensions.states solution.states.dimensions.states(end)]), method.interpolate.basis.portfolio_problem,method.interpolate.extrapolate.portfolio_problem); 
      case 'compecon'  
        coefficients =     funbas(solution.portfolio_problem.basis, gridmake(funnode(solution.portfolio_problem.basis)), zeros(1,length(model.variables.states)+1))\reshape(values,[],1);
        func = @(states)   funbas(solution.portfolio_problem.basis, states,                                              zeros(1,length(model.variables.states)+1))*coefficients; %takes stack     
    end
        
end

end
