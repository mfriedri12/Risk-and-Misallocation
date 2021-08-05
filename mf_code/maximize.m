%% maximize  
%-------------------------------------------------------------------------%
% Now it still only works over 2 variables 
%-------------------------------------------------------------------------%
function [argmax, maximum] = maximize(method, objective, search_min, search_max, additional_arguments)

%% check search interval 
search_max = max(search_max,search_min); 

%% find maximum 
switch method.maximize.algorithm 
    
    case {'golden search'} % hill climibing; uses bellman; really only works for models with one endogenous state variable right now 
       
      [argmax, maximum] = golden_mrf(objective, method.maximize.threshold, search_min, search_max, additional_arguments);  
    
      %return endpoint if larger 
        if method.maximize.checkendpoint
            fa = feval(objective,search_min,additional_arguments);
            selector = fa > maximum;
            argmax = selector.*search_min + (1-selector).*argmax;
            maximum = selector.*fa + (1-selector).*maximum;
                       
            fb = feval(objective,search_max,additional_arguments);
            selector = fb > maximum;
            argmax = selector.*search_max + (1-selector).*argmax;
            maximum = selector.*fb + (1-selector).*maximum;
        end 
end

