%% maximize  
%-------------------------------------------------------------------------%
% pretty slow right now because I haven't gotten aroudn to paralelling
% (dissapointed that golden from compecon seems to to paralellize incorrectly) 
%-------------------------------------------------------------------------%

function [argmax, maximum] = maximize(model, method, solution, objective)


switch method.maximize.algorithm 
    
    case {'golden search'} % hill climibing; uses bellman; really only works for models with one endogenous state variable 
       
       for i=1:length(model.variables.endogenous) 
            a = model.grid.min(i)*ones(solution.states.dimensions.states); 
            b = max(model.conditions.max.b(solution, model.parameters),a);  % model.grid.max(i)*ones(solution.states.dimensions.states); % must not be so big that c is negative 
            argmax.(model.variables.endogenous(i)) = zeros(solution.states.dimensions.states); 
            maximum = zeros(solution.states.dimensions.states); 
            for j=1:prod(solution.states.dimensions.states)
                states.z = solution.states.ndgrid.z(j);
                states.b = solution.states.ndgrid.b(j);
                [argmax.(model.variables.endogenous(i))(j), maximum(j)] = golden(objective, a(j), b(j), states);
            end 
       end
               
end


        