%% interpolate 
%-------------------------------------------------------------------------%
% interpolate a state-space sized guess either for a value or a control (policy) function. 
% can do this either with a basis or with griddedInterpolent; if basis has optimal nodes, update the state space to 
% correspond to those optimal nodes; this is why this function returns "solution" as well as guess.  
%-------------------------------------------------------------------------%

function basis = interpolate(method, solution, values, special_grid) 

switch method.interpolate.algorithm
    
    case 'griddedInterpolant'
       
        if exist('special_grid','var')
            basis = griddedInterpolant(special_grid,        values,method.interpolate.basis,method.interpolate.extrapolation_method ); 
        else 
            basis = griddedInterpolant(solution.states.cell,values,method.interpolate.basis,method.interpolate.extrapolation_method ) ;
        end
        
end

end
