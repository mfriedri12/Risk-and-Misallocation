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
solution.states.dimensions.exogenous =  method.grid.exogenous.nodes*ones(1,length(model.variables.exogenous ));
solution.states.dimensions.states =    [solution.states.dimensions.exogenous, solution.states.dimensions.endogenous];

% get appropriate grids for states 
switch method.update_distribution.exogenous_transition
  case 'gaussian quadrature'
    for i=1:length(model.variables.exogenous) 
      solution.states.grid.(model.variables.exogenous(i)) = linspace(model.grid.min.(model.variables.exogenous(i)),model.grid.max.(model.variables.exogenous(i)), solution.states.dimensions.exogenous(i))'; % for looking at mainly    
    end
 case 'tauchen'
   if length(model.variables.exogenous)>1
       disp('tauchen not ok for more than one exogenous variable'); pause;
   end
   solution.states.grid.(model.variables.exogenous(1)) = tauchen(0,model.parameters.rho.(model.variables.exogenous(1)), model.parameters.sigma.(model.variables.exogenous(1)),method.grid.exogenous.nodes);       
end
   

switch method.grid.endogenous.spacing
  case 'compecon'
    method.interpolate.basis.value_function = extractBetween(method.interpolate.basis.value_function ,1,4);
    method.interpolate.basis.value_function =  method.interpolate.basis.value_function{1};
    for i = 1:length(model.variables.states)
      minimums(i) = model.grid.min.(model.variables.states(i));
      maximums(i) = model.grid.max.(model.variables.states(i));
    end
    solution.value_function.basis = fundefn( method.interpolate.basis.value_function ,solution.states.dimensions.states, minimums, maximums,  method.interpolate.degree); %solution.states.dimensions.states
    nodes = funnode(solution.value_function.basis); 
    for i = 1:length(model.variables.endogneous)
      solution.states.grid.(model.variables.endogneous(i)) = nodes{i}; 
    end       
   otherwise
    for i=1:length(model.variables.endogenous) 
      switch method.grid.endogenous.spacing
        case 'even'
          solution.states.grid.(model.variables.endogenous(i)) = linspace(model.grid.min.(model.variables.endogenous(i)),model.grid.max.(model.variables.endogenous(i)), solution.states.dimensions.endogenous(i))'; % for looking at mainly    
        case 'custom'
          solution.states.grid.(model.variables.endogenous(i)) = model.makegrid.(model.variables.endogenous(i))(solution.states.dimensions.endogenous(i), model.parameters);
        case 'log1'
          solution.states.grid.(model.variables.endogenous(i)) = (exp(linspace(0,log(model.grid.max.(model.variables.endogenous(i)) - model.grid.min.(model.variables.endogenous(i))+1), solution.states.dimensions.endogenous(i)))-1+model.grid.min.(model.variables.endogenous(i)))';  % set up quadruple exponential grid    
        case 'log2'
          solution.states.grid.(model.variables.endogenous(i)) = (exp(exp(linspace(0,log(log(model.grid.max.(model.variables.endogenous(i)) - model.grid.min.(model.variables.endogenous(i))+1)+1), solution.states.dimensions.endogenous(i)))-1)-1+model.grid.min.(model.variables.endogenous(i)))';  % set up quadruple exponential grid    
        case 'log3'
         solution.states.grid.(model.variables.endogenous(i)) = (exp(exp(exp(linspace(0,log(log(log(model.grid.max.(model.variables.endogenous(i)) - model.grid.min.(model.variables.endogenous(i))+1)+1)+1), solution.states.dimensions.endogenous(i)))-1)-1)-1+model.grid.min.(model.variables.endogenous(i)))';  % set up quadruple exponential grid    
        case 'log4'
         solution.states.grid.(model.variables.endogenous(i)) = (exp(exp(exp(exp(linspace(0,log(log(log(log(model.grid.max.(model.variables.endogenous(i)) - model.grid.min.(model.variables.endogenous(i))+1)+1)+1)+1), solution.states.dimensions.endogenous(i)))-1)-1)-1)-1+model.grid.min.(model.variables.endogenous(i)))';  % set up quadruple exponential grid    
      end 
    end
end 


% reshape this grid into other formats for use in other places -- now stack only! 
for i=1:length(model.variables.states) 
    solution.states.cell{i} = solution.states.grid.(model.variables.states(i)) ; % for interpolant / griddedInterpolant 
    solution.states.cell_prime{i} = solution.states.cell{i}; % for integrate 
end 
solution.states.stack.(model.variables.states) = gridmake(solution.states.cell); % if memory a problem could replace all incidences of this with combos of below 
for i=1:length(model.variables.states) 
   solution.states.stack.(model.variables.states(i)) = solution.states.stack.(model.variables.states)(:,i) ; % for interpolant / griddedInterpolant 
end 
if length(model.variables.endogenous) > 1 && contains(method.update_controls.algorithm,'VFI') % need that extra big grid for maximizing only works for 2 variables
  solution.states.cell_plus = solution.states.cell; 
  solution.states.cell_plus{length(solution.states.cell)+1} = solution.states.grid.(model.variables.endogenous(2));    
  solution.states.stackplus.(model.variables.states) = gridmake(solution.states.cell_plus);
  for i=1:length(model.variables.states) 
    solution.states.stackplus.(model.variables.states(i)) = solution.states.stackplus.(model.variables.states)(:,i) ; % for interpolant / griddedInterpolant 
  end 
  solution.states.stackplus.(strcat(model.variables.endogenous(2),model.variables.endogenous(2))) = solution.states.stackplus.(model.variables.states)(:,i+1);
end 

% approximate exogenous state shocks 
switch method.integrate.algorithm
  case 'gaussian quadrature' 
    if length(model.variables.exogenous)==1
      model.parameters.covariance =  model.parameters.sigma.(model.variables.exogenous(1));
    elseif length(model.variables.exogenous)==2
       model.parameters.covariance = [model.parameters.sigma.(model.variables.exogenous(1)) model.parameters.correlation*model.parameters.sigma.(model.variables.exogenous(1))*model.parameters.sigma.(model.variables.exogenous(2));...
                                      model.parameters.correlation*model.parameters.sigma.(model.variables.exogenous(1))*model.parameters.sigma.(model.variables.exogenous(2)) model.parameters.sigma.(model.variables.exogenous(2)) ];
    end
    switch method.integrate.gaussian  
      case 'qnwnorm'
        [shocks, solution.states.weights] = qnwnorm(method.integrate.nodes*ones(1,length(model.variables.exogenous)), zeros(1,length(model.variables.exogenous)), model.parameters.covariance); 
    end
    shocks = transpose(shocks);    % shocks should be a horizontal vector
    for i=1:length(model.variables.exogenous)
      solution.states.stack.prime(:,i) =  reshape(transpose(model.parameters.rho.(model.variables.exogenous(i))*solution.states.stack.(model.variables.exogenous(i)) + shocks(i,:)),[],1); %reshape(transpose(model.exogenous.process.(model.variables.exogenous(i))(solution.states.stack.(model.variables.states)(:,i), shocks(i,:), model.parameters)),[],1);
    end 
    for i=1:length(model.variables.endogenous)
      solution.states.stack.prime(:,length(model.variables.exogenous)+i) = ...
          reshape(repmat(transpose(solution.states.stack.(model.variables.states)(:,length(model.variables.exogenous)+i)),size(shocks,2),1),[],1); 
    end
end  


%% distribution  
% guess distribution, get exogenous transition matricies 
solution.distribution.values = ones(prod(solution.states.dimensions.states),1)/prod(solution.states.dimensions.states); % guess uniform distribution 
switch method.update_distribution.exogenous_transition
  case 'gaussian quadrature'
    solution.distribution.transition.exogenous=sparse(prod(solution.states.dimensions.exogenous),prod(solution.states.dimensions.exogenous));
    
    for i=1:length(model.variables.exogenous) 
      solution.states.stack.exogenous{i} = model.parameters.rho.(model.variables.exogenous(i))*solution.states.cell{i};
    end
    solution.states.stack.exogenous = gridmake(solution.states.stack.exogenous);
     
    for i = 1:length(model.variables.exogenous)
      minimums(i) = model.grid.min.(model.variables.exogenous(i));
      maximums(i) = model.grid.max.(model.variables.exogenous(i));
    end
    exogenous_basis = fundefn('spli' ,solution.states.dimensions.exogenous, minimums, maximums,  1); %solution.states.dimensions.states
    
    for i=1:size(shocks,2)
      next = solution.states.stack.exogenous + repmat(shocks(:,i)',size(solution.states.stack.exogenous,1),1);
      for j=1:length(model.variables.exogenous)
        next(:,j) = max(next(:,j), min(solution.states.grid.(model.variables.exogenous(j))));
        next(:,j) = min(next(:,j), max(solution.states.grid.(model.variables.exogenous(j))));
      end
      solution.distribution.transition.exogenous = solution.distribution.transition.exogenous + solution.states.weights(i)*funbas(exogenous_basis, next); 
    end
    solution.distribution.transition.exogenous = full(solution.distribution.transition.exogenous); % can figure out how to put sparse matricies into update_distribution later 
  
 case 'tauchen'
   if length(model.variables.exogenous)>1
       disp('tauchen not ok for more than one exogenous variable'); pause;
   end
   [~,solution.distribution.transition.exogenous] = tauchen(0,model.parameters.rho.(model.variables.exogenous(1)), model.parameters.sigma.(model.variables.exogenous(1)),method.grid.exogenous.nodes );       
end


%% prices 
% guess prices & define interval
for i=1:length(model.variables.prices)
    solution.prices.values.(model.variables.prices(i)) = 1.00;  
    solution.prices.interval.(model.variables.prices(i))(1) = (1/method.update_prices.interval)*solution.prices.values.(model.variables.prices(i));
    solution.prices.interval.(model.variables.prices(i))(2) = (method.update_prices.interval)*solution.prices.values.(model.variables.prices(i));
    if method.update_prices.check_interval
       solution.prices.values.(model.variables.prices(i)) = solution.prices.interval.(model.variables.prices(i))(1); 
    end
    solution.prices.algorithm.excess.(model.variables.prices(i)) = 1e5; 
    solution.prices.algorithm.prices.(model.variables.prices(i)) = solution.prices.values.(model.variables.prices(i)); 
    solution.prices.algorithm.update.(model.variables.prices(i)){1} = 'bisection';
end


%% controls
% guess controls They are stack sized 
for i=1:length(model.variables.endogenous)
    solution.controls.values.(model.variables.endogenous(i)) = (max(solution.states.stack.(model.variables.endogenous(i)),0)).^.8; % zeros(size(solution.states.stack.(model.variables.endogenous(1))));  
end

switch method.interpolate.toolkit.control_function
  case 'compecon'
    for i = 1:length(model.variables.endogneous)
      solution.controls.basis.(model.variables.endogenous(i)) =  solution.value_function.basis;  
    end       
    solution.portfolio_problem.basis = fundefn( method.interpolate.basis.value_function , ...
                                             [solution.states.dimensions.states solution.states.dimensions.states(end)],...
                                             [minimums minimums(end)], ...
                                             [maximums maximums(end)], ...
                                             method.interpolate.degree); %solution.states.dimensions.states
end                                      

for i=1:length(model.variables.other)
    solution.controls.values.(model.variables.other(i)) = max(model.other.(model.variables.other(i))(solution.controls.values, solution.states.stack, solution.prices.values, model.parameters),.01);
end


%% value function 
solution.value_function.values = model.functions.objective(solution.controls.values.(model.variables.other(i)),model.parameters.gamma); %zeros(size(solution.states.stack.(model.variables.endogenous(1)))); 

end