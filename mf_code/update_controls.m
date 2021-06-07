%% update controls 
%-------------------------------------------------------------------------%
% given prices, solve single agent problem. all algorithms follow these steps 
% 0. guess the VALUES of the value or the control (policy) function on the state space grid (first guess made in "setup_solution")
% 1. usually, *interpolate* the above to get a guess for the value or the control (policy) FUNCTION. (Can avoid interpolaiton if you integrate out exogenous states via a markov approximation of the exogenous state; butthis solution will also usually be less accurate)
% 2. update your VALUES using the bellman equation or some deriviative of it (FOC, EC, EE) that is also a contraction
%      2.1 usually you'll need to *integrate* over exogenous states to numerically calculate the expected value  
%      2.2 sometimes you'll need to nonlinearly *maximize* to get the updated VALUES, using a hill-climbing algorithm
%      2.2 other tiems you'll need to take a derivative and *root-find*
%      2.3 usually you'll need to enforce borrowing constraints 
% 3. check for convergence, that the distance between your guess & update is less than thresholds chosen in setup_method
% 4. if not converged, set guess = update and repeat until convergence 

% *interpolate* takes in values guesses and returns an object with 4 fields: values, basis, coefficients, and nodes, all defined over the state space   
% *integrate* takes in an object returned by interpolate and returns expected values integrating over exogenous states  
% *solve* takes in a function and returns its roots OR its maximization  

%-------------------------------------------------------------------------%

function [solution] = update_controls(model, method, solution)

%% reset 
if method.update_controls.reset 
    for i=1:length(model.variables.controls)
        solution.controls.values.(model.variables.controls(i)) = ones(solution.states.dimensions.states);  % guess controls as 1 
    end
    solution.value_function.values =  ones(solution.states.dimensions.states) ; % guess valus as 1 
end


%% solve for controls that are endogenous states by following steps 1-4 
iteration = 1; 
solution.accuracy.controls = 1e5;
last = solution.accuracy.controls;
while method.update_controls.threshold < solution.accuracy.controls && iteration < method.update_controls.iterations     
    
    switch method.update_controls.algorithm
        
        case 'VFI'
            expected_value_function = interpolate(method, solution, solution.value_function.values); % 1 
            expected_value_function = integrate(model, method, solution, expected_value_function, 'basis'); % 2.1 
            bellman = @(controls, states) model.conditions.b1(states, controls, solution.prices.values, model.parameters) + model.parameters.beta*expected_value_function(states.z, controls) ; % 2.2  
            [solution.controls.values, update] = maximize(model, method, solution, bellman); % 2.2
            solution.accuracy.controls = max(abs(update - solution.value_function.values), [], 'all'); % 3
            solution.value_function.values = update; % 4     
        
        case 'BLPT' % I still hate this algorithm
            % first update c on a' k', z grid 
            solution.controls.basis.c = interpolate(method, solution, solution.controls.values.c); % 1 
            rhs = model.conditions.rhs(solution.states.ndgrid.z, solution.states.ndgrid.a, solution.states.ndgrid.k, solution.controls.basis.c, model.parameters,solution.prices.values); 
            rhs = interpolate(method, solution, rhs);
            rhs = integrate(model, method, solution, rhs,'values');
            lhs = model.functions.objective_foc_inv(rhs,model.parameters.gamma); %still normal grid size 
            lhs = interpolate(method, solution, lhs); % easier than subseting? 
            for i=1:2 
                qtype = strcat('q',string(i));
                
                % second find a' on k' z grid given c                                                                                                    endogenous state/control associated with second endogenous state/control, exogenous state, given c 
                eulercombo = model.conditions.eulercombo.(qtype)(solution.states.ndgrid.z, solution.states.ndgrid.a, solution.states.ndgrid.k,solution.controls.basis.c, model.parameters,solution.prices.values); % model.functions.objective_foc(c(z,a,k))*(model.functions.profit(k,z)
                eulercombo = interpolate(method, solution, eulercombo);
                eulercombo = integrate(model, method, solution, eulercombo,'values'); 
                [~, a_index.(qtype)] = min(abs(eulercombo),[],2); % feel like I should do a root find? but that is so slow? and would it help? 
                %[~, k_index.(qtype)] = min(abs(eulercombo),[],3); % alt pick k 
                
                % third calculate resources 
                a.(qtype) = reshape(squeeze(solution.states.grid.a(a_index.(qtype))),[],1); 
                %a.(qtype) = reshape(squeeze(solution.states.ndgrid.a(:,:,1)),[],1); % alt pick k 
                k.(qtype) = reshape(squeeze(solution.states.ndgrid.k(:,1,:)),[],1);
                %k.(qtype) = reshape(squeeze(solution.states.grid.k(k_index.(qtype))),[],1);   % alt pick k 
                z.(qtype) = reshape(squeeze(solution.states.ndgrid.z(:,1,:)),[],1);
                c.(qtype) = reshape(squeeze(lhs(z.(qtype), a.(qtype), k.(qtype))),[],1);
                r.(qtype) = model.functions.resources.a.(qtype)(c.(qtype),a.(qtype),k.(qtype),model.parameters); % z by k' grid      
                
                % fourth subset to just the rz_grid a & interpolate     
                [r.(qtype), r_index.(qtype)] = sort(r.(qtype));
                c.(qtype)  =     interpolate(method, solution, c.(qtype)(r_index.(qtype)),  r.(qtype)); 
                k.(qtype)  =     interpolate(method, solution, k.(qtype)(r_index.(qtype)),  r.(qtype)); %kprime 
                a.(qtype)  =     interpolate(method, solution, a.(qtype)(r_index.(qtype)),  r.(qtype)); %aprime 
                
                %fifth update policy functions based on interpolation
                c.(qtype) = max(c.(qtype)(model.conditions.resources.b.(qtype)(solution,model.parameters)),0.01);
                k.(qtype) = max(k.(qtype)(model.conditions.resources.b.(qtype)(solution,model.parameters)),model.grid.min(1));
                a.(qtype) = max(a.(qtype)(model.conditions.resources.b.(qtype)(solution,model.parameters)),model.grid.min(2));
            end 
            q2 =  k.q2 < solution.states.ndgrid.k;
            q1 =  1 -  q2; %this is default now for small numerical erros 
            update =                     c.q1.*q1 + c.q2.*q2;
            solution.controls.values.k = k.q1.*q1 + k.q2.*q2;
            solution.controls.values.a = a.q1.*q1 + a.q2.*q2;
            solution.accuracy.controls = max(abs(update - solution.controls.values.c), [], 'all'); % 3
            solution.controls.values.c = update; % 4  
            
        case 'EE' 
            % interpolate control functions 
            for i=1:length(model.variables.endogenous)
                 solution.controls.basis.(model.variables.controls(i)) = interpolate(method, solution, solution.controls.values.(model.variables.controls(i))); % 1 - interpolate policy function
            end
            % update control functions via euler equations 
            for i=1:length(model.variables.endogenous)   
                expected_ee_RHS = interpolate(method, solution, model.conditions.eeRHS.(model.variables.controls(i))(solution, model.parameters)); % 1 
                expected_ee_RHS = integrate(model, method, solution, expected_ee_RHS,'values'); % 2.1
                update.(model.variables.controls(i)) = model.conditions.ee.(model.variables.controls(i))(solution,model.parameters,expected_ee_RHS);  % 2
                update.(model.variables.controls(i))(update.(model.variables.controls(i))< model.grid.min(i)) = model.grid.min(i); % 2.3 
                solution.accuracy.controls = max(abs(update.(model.variables.controls(i)) - solution.controls.values.(model.variables.controls(i))), [], 'all'); % 3
                solution.controls.values.(model.variables.controls(i)) = update.(model.variables.controls(i)); % 4
            end
            
    end

    if mod(iteration,method.update_controls.print)==0
        fprintf('iteration: %i, precision: %i \n', iteration, solution.accuracy.controls);
        if last < solution.accuracy.controls; fprintf('WARNING! algorithm exploding at iteration: %i. Paused. \n',iteration); pause; end
        last = solution.accuracy.controls;
    end
    iteration = iteration + 1; 
end


%% get "other" controls 
% this is usually consumption 

for i=1:length(model.variables.other)
    solution.controls.values.(model.variables.other(i)) = model.other.(model.variables.other(i))(solution, model.parameters); 
end


%% get control funciton basis 
% sometimes this is incidental in the algorithm but not always 
for i=1:length(model.variables.controls)
    solution.controls.basis.(model.variables.controls(i)) = interpolate(method, solution, solution.controls.values.(model.variables.controls(i)));
end

 
end
        