%% update distribution 
%-------------------------------------------------------------------------%
% given controls, update distribuiton
% more useful view that I might rewrite as --
% reshape(solution.distribution.values, solution.states.dimensions.states)
%-------------------------------------------------------------------------%

function [solution] = update_distribution (model, method, solution)


%% get endogenous transition matricies

solution.distribution.transition.endogenous = zeros(prod(solution.states.dimensions.states));
spacer = eye(prod(solution.states.dimensions.exogenous));
s = 1; 
for i=1:length(solution.distribution.values)
    
     for j=1:length(model.variables.endogenous)
        
        % find transition probability from state point i to grid space j  
        transition_i_to_j = zeros(solution.states.dimensions.endogenous(j),1) ;
        lookup = solution.controls.basis.(model.variables.endogenous(j))(solution.states.stack(i,:));
        h = find(lookup<=solution.states.grid.(model.variables.endogenous(j)),1); % first index greater 
        if h==1  % lookup value below grid 
            transition_i_to_j(1) = 1; 
        elseif isempty(h) % lookup value above grid 
            transition_i_to_j(end) = 1; 
        else 
            transition_i_to_j(h-1) = 1-((lookup-solution.states.grid.(model.variables.endogenous(j))(h-1))/(solution.states.grid.(model.variables.endogenous(j))(h) - solution.states.grid.(model.variables.endogenous(j))(h-1)));
            transition_i_to_j(h)   = 1-((solution.states.grid.(model.variables.endogenous(j))(h)-lookup)  /(solution.states.grid.(model.variables.endogenous(j))(h) - solution.states.grid.(model.variables.endogenous(j))(h-1)));
        end
       
        % add these probabilities to the previous ones  
        if j==1 
            transition_i_to_i = transition_i_to_j; 
        else
            transition_i_to_i = reshape(transition_i_to_i*transition_i_to_j',[],1);
        end
          
    end   
    
    % add exogenous state zeros or duplicate for full matrix 
    
    switch method.update_distribution.algorithm
       case 'tan' 
          solution.distribution.transition.endogenous(i,:) = reshape(transpose(transition_i_to_i*spacer(s,:)),[],1); 
        case 'default'
          solution.distribution.transition.endogenous(i,:) = reshape(transpose(transition_i_to_i*solution.distribution.transition.exogenous(s,:)),[],1); 
    end
    if s < prod(solution.states.dimensions.exogenous) 
      s = s + 1;
    else 
      s = 1;
    end
    
end



%% solve for distribution 
% start with some guess, and then iterate forward to solve for distribution 

if method.update_distribution.reset || ~all(solution.distribution.values > method.update_distribution.threshold )
    solution.distribution.values = ones(prod(solution.states.dimensions.states),1)/prod(solution.states.dimensions.states); % guess uniform distribution 
end

iteration = 0; 
solution.distribution.accuracy = 1e5;
while method.update_distribution.threshold < solution.distribution.accuracy && iteration < method.update_distribution.iterations   
    
    switch method.update_distribution.algorithm
        case 'tan'
            intermediate = reshape(transpose(solution.distribution.transition.endogenous)*solution.distribution.values, prod(solution.states.dimensions.exogenous), prod(solution.states.dimensions.endogenous));
            distribution = reshape(transpose(solution.distribution.transition.exogenous )*intermediate, prod(solution.states.dimensions.states),1);
        case 'default'
            distribution = transpose(solution.distribution.transition.endogenous)*solution.distribution.values;
    end
    
solution.distribution.accuracy = abs(max(distribution - solution.distribution.values));
solution.distribution.values = distribution;  
 
iteration = iteration + 1;
if mod(iteration,method.update_distribution.print)==0; fprintf('iteration: %i, precision: %i \n', iteration, solution.distribution.accuracy); end  
end


 
end