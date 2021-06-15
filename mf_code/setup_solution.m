%% setup solution 
%-------------------------------------------------------------------------%
% initialize all the solution parameters
%-------------------------------------------------------------------------%

function [solution] = setup_solution(model, method)


%% states 
% the state space is defined in several ways because different parts of the
% program need it different ways. If memory starts to be a problem should
% look at this and see what parts can cut out 

solution.states.dimensions.endogenous = method.grid.endogenous.nodes*ones(1,length(model.variables.endogenous));
solution.states.dimensions.exogenous = method.grid.exogenous.nodes*ones(1,length(model.variables.exogenous ));
solution.states.dimensions.states = [solution.states.dimensions.exogenous,  solution.states.dimensions.endogenous];

% get appropriate grid for endogenous states
for i=1:length(model.variables.endogenous) 
   switch method.grid.endogenous.spacing
     case 'even'
       solution.states.grid.(model.variables.endogenous(i)) = linspace(model.grid.min(i),model.grid.max(i), solution.states.dimensions.endogenous(i))'; % for looking at mainly 
     case 'log1'
       solution.states.grid.(model.variables.endogenous(i)) = (exp(linspace(0,log(model.grid.max(i) - model.grid.min(i)+1), solution.states.dimensions.endogenous(i)))-1+model.grid.min(i));  % set up quadruple exponential grid    
     case 'log2'
       solution.states.grid.(model.variables.endogenous(i)) = (exp(exp(linspace(0,log(log(model.grid.max(i) - model.grid.min(i)+1)+1), solution.states.dimensions.endogenous(i)))-1)-1+model.grid.min(i));  % set up quadruple exponential grid    
     case 'log3'
       solution.states.grid.(model.variables.endogenous(i)) = (exp(exp(exp(linspace(0,log(log(log(model.grid.max(i) - model.grid.min(i)+1)+1)+1), solution.states.dimensions.endogenous(i)))-1)-1)-1+model.grid.min(i));  % set up quadruple exponential grid    
     case 'log4'
       solution.states.grid.(model.variables.endogenous(i)) = (exp(exp(exp(exp(linspace(0,log(log(log(log(model.grid.max(i) - model.grid.min(i)+1)+1)+1)+1), solution.states.dimensions.endogenous(i)))-1)-1)-1)-1+model.grid.min(i));  % set up quadruple exponential grid    
   end
end 

% get appropriate grid for exogenous states 
for i=1:length(model.variables.exogenous)  
   switch method.grid.exogenous.spacing
     case 'custom'
       solution.states.grid.(model.variables.exogenous(i)) = model.exogenous.grid.(model.variables.exogenous(i))(method.grid.exogenous.nodes, model.parameters);
   end
end

% reshape this grid into other formats for use in other places 
for i=1:length(model.variables.states) 
    solution.states.cell{i} = solution.states.grid.(model.variables.states(i)) ; % for interpolant / griddedInterpolant 
    solution.states.cell_prime{i} = solution.states.cell{i}; % for integrate 
    twist = reshape(solution.states.grid.(model.variables.states(i)) ,[ones(1,i-1),solution.states.dimensions.states(i),ones(1,length(model.variables.states)-i)]);
    solution.states.ndgrid.(model.variables.states(i)) = repmat(twist, [solution.states.dimensions.states(1:(i-1)),1, solution.states.dimensions.states((i+1):length(model.variables.states))]); % for functions 
    solution.states.stack(:,i) = reshape(solution.states.ndgrid.(model.variables.states(i)), [],1);
end

% need that extra big grid for maximizing only works for 2 variables 
if length(model.variables.endogenous) > 1 && contains(method.update_controls.algorithm,'VFI') 
  solution.states.cell_plus = solution.states.cell; 
  solution.states.cell_plus{length(solution.states.cell)+1} = solution.states.grid.(model.variables.endogenous(2));   
  for i=1:length(model.variables.states) 
    solution.states.ndgridplus.(model.variables.states(i)) = repmat(solution.states.ndgrid.(model.variables.states(i)),[ones(1,length(model.variables.states)) solution.states.dimensions.endogenous(2)]);
  end
  twist = reshape(solution.states.grid.(model.variables.endogenous(2)),[ones(1,length(model.variables.states)) solution.states.dimensions.endogenous(2)]);
  solution.states.ndgridplus.(strcat(model.variables.endogenous(2),model.variables.endogenous(2))) = repmat(twist,[solution.states.dimensions.states 1]);
end


% approximate exogenous state shocks 
for i=1:length(model.variables.exogenous) 
    switch method.integrate.algorithm
       case 'gaussian quadrature' 
           switch method.integrate.gaussian  
                case 'qnwnorm'; [shocks, solution.states.weights.(model.variables.exogenous(i))] = qnwnorm(method.integrate.nodes, 0, model.parameters.sigma); 
           end
           shocks = shocks';    % shocks should be a horizontal vector
    end  
    solution.states.cell_prime{i} = reshape(transpose(model.exogenous.process.z(solution,shocks, model.parameters)),[],1);
end 
    

%% distribution  
% guess distribution, get exogenous transition matricies 

solution.distribution.values = ones(prod(solution.states.dimensions.states),1)/prod(solution.states.dimensions.states); % guess uniform distribution 
for i=1:length(model.variables.exogenous)   
    [~,solution.distribution.transition.exogenous] = model.exogenous.transition.z(method.grid.exogenous.nodes, model.parameters);
end


%% prices 

% guess prices & define interval
for i=1:length(model.variables.prices)
    solution.prices.values.(model.variables.prices(i)) = 1.25;  
    solution.prices.interval.(model.variables.prices(i))(1) = (1/method.update_prices.interval)*solution.prices.values.(model.variables.prices(i));
    solution.prices.interval.(model.variables.prices(i))(2) = (method.update_prices.interval)*solution.prices.values.(model.variables.prices(i));
    if method.update_prices.check_interval
       solution.prices.values.(model.variables.prices(i)) = solution.prices.interval.(model.variables.prices(i))(1); 
    end
end
%solution.prices.algorithm.allexcess = zeros(15,2);

%% controls & values 
% guess controls & value function. They are ndgrid sized 

for i=1:length(model.variables.endogenous)
    solution.controls.values.(model.variables.endogenous(i)) = (solution.states.ndgrid.(model.variables.endogenous(i))).^.8; % (ones(solution.states.dimensions.states);  % guess controls as 1 (solution.states.ndgrid.z.*solution.states.ndgrid.a.*solution.states.ndgrid.k).^.25
end

for i=1:length(model.variables.other)
    solution.controls.values.(model.variables.other(i)) = max(model.other.(model.variables.other(i))(solution, model.parameters),.01); % (ones(solution.states.dimensions.states);  % guess controls as 1 (solution.states.ndgrid.z.*solution.states.ndgrid.a.*solution.states.ndgrid.k).^.25
end

solution.value_function.values = ones(solution.states.dimensions.states);  %solution.controls.values.(model.variables.other(i))/(1-.9); %ones(solution.states.dimensions.states) ; % guess valus as 1 
 
 
end