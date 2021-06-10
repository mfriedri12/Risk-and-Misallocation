%% setup solution 
%-------------------------------------------------------------------------%
% define all model parameters & model functions; ideally should be no
% model-specific info in the rest of the code 
%-------------------------------------------------------------------------%

function [model] = setup_model


%% parameters
% these are for use in the functions below 

model.parameters.alpha = 1/3; % capital share 
model.parameters.beta = .95; % discount  
model.parameters.gamma = 2; % risk aversion / IES 
model.parameters.delta = .05; % capital depreciation
model.parameters.eta = .85; % DRS at individual level 
model.parameters.lambda = .25; % haircut 
model.parameters.theta = 1; % borrowing constraint 
model.parameters.rho = .9; % productivity persistance 
model.parameters.sigma = .25; % productivity innovation variance


%% variables  
% this just records the basic dimensions of the model

model.variables.exogenous = ['z']; %  ["z: productivity"]; exogenous states; assume these are AR(1)s for now 
model.variables.endogenous = ['a','k']; % ["a: liquid assets" "k: capital"];  % engogenous states = intersection of states & controls. the key variables 
model.variables.other = ['c']; %[ "c: consumption"]
model.variables.prices = ['w','r']; % ["w: wages" "r: liquid asset return"];
model.variables.states = ['z', 'a','k']; % [model.variables.exogenous, model.variables.endogenous];
model.variables.controls = ['a','k','c']; %[model.variables.endogenous, model.variables.other]; 


%% grid
% this defines endogneous states grid. States are listed in order
% noted above. Agents are not allowed to choose a state variable outside of this
% grid so that's why it's part of model 

model.grid.min = [0.01 0.01];  % borrowing constraints
model.grid.max = [500 500]; % this is not so important but we're still gonna call it model specific 


%% exogenous state processes
%used by setup solution to get grid for these and also by *integrate* to calculate expected values

model.exogenous.process.z = @(solution, shocks, parameters) (solution.states.grid.z.^parameters.rho).*exp(shocks); 
model.exogenous.grid.z = @(gridsize, parameters) exp(tauchen(0,parameters.rho,parameters.sigma,gridsize)); %   *parameters.sigma/sqrt(1-parameters.rho^2));
model.exogenous.transition.z = @(gridsize, parameters) tauchen(0,parameters.rho, parameters.sigma, gridsize); % second argument 


%% objective & constraint functions 
% first define subfunctions naively -- as a function of model-specific states, parameters

% objective 
model.functions.objective = @(c, gamma) (c.^(1-gamma))./(1-gamma); 
model.functions.objective_foc = @(c, gamma) (c.^(-gamma)); 
model.functions.objective_foc_inv = @(argument, gamma) argument.^(1/(-gamma));

% the price function                                     
model.functions.q = @(k, kk, delta, lambda) (kk >= (1-delta)*k) * 1  + (kk <  (1-delta)*k) * (1-lambda);

% the optimized labor functions 
model.functions.labor = @(w, k, z, alpha, eta) w./ (((z.^(1-(1-alpha)*eta)).* (k.^(alpha*eta)).*  ((1-alpha)*eta)).^    (1/((1-alpha)*eta - 1)));  

% the optimized profit functions 
model.functions.profit =     @(w, k, z, alpha, eta)                                       z.* (k.^((alpha*eta)/(1-(1-alpha)*eta))) .*  (1/w).^(((1-alpha)*eta)/(1-(1-alpha)*eta)) .*  ((1-alpha)*eta).^(((1-alpha)*eta )/(1-(1-alpha)*eta)) .*  (1-  ((1-alpha)*eta) );
model.functions.profit_foc = @(w, k, z, alpha, eta)     ((alpha*eta)/(1-(1-alpha)*eta)).* z.* (k.^((eta-1)/(1-(1-alpha)*eta)))     .*  (1/w).^(((1-alpha)*eta)/(1-(1-alpha)*eta)) .*  ((1-alpha)*eta).^(((1-alpha)*eta )/(1-(1-alpha)*eta)) .*  (1-  ((1-alpha)*eta) );

% the constraint solved for c
model.functions.c = @(w, r, a, k, z, aa, kk, alpha, delta, eta, lambda, theta) w  + model.functions.profit(w, k, z, alpha, eta)  + r*(a - theta*lambda*(1-delta)*k)  - (aa - theta*lambda*(1-delta)*kk) - model.functions.q(k, kk, delta, lambda).*(kk - (1-delta)*k) ; 

% resources 
model.functions.resources.a.q1 = @(c,a,k,lambda) c + a + 1*k; 
model.functions.resources.a.q2 = @(c,a,k,lambda) c + a + lambda*k; 
model.functions.resources.b.q1 = @(w,r,a,k,z,alpha,delta,eta,lambda) w + model.functions.profit(w, k, z, alpha, eta) + r.*a + 1.*(1-delta).*k ; 
model.functions.resources.b.q2 = @(w,r,a,k,z,alpha,delta,eta,lambda) w + model.functions.profit(w, k, z, alpha, eta) + r.*a + lambda.*(1-delta).*k; 

model.functions.resources.a.q1 = @(c,a,k,parameters) c + a + 1*k; 
model.functions.resources.a.q2 = @(c,a,k,parameters) c + a + parameters.lambda*k; 
model.conditions.resources.b.q1 = @(solution,parameters) model.functions.resources.b.q1(solution.prices.values.w,solution.prices.values.r,solution.states.ndgrid.a,solution.states.ndgrid.k,solution.states.ndgrid.z,parameters.alpha,parameters.delta,parameters.eta,parameters.lambda) ;
model.conditions.resources.b.q2 = @(solution,parameters) model.functions.resources.b.q2(solution.prices.values.w,solution.prices.values.r,solution.states.ndgrid.a,solution.states.ndgrid.k,solution.states.ndgrid.z,parameters.alpha,parameters.delta,parameters.eta,parameters.lambda) ;




%% optimal conditions    
% next define optimal conditions, which we will need for update_controls, as functions of states, controls, prices, parameters OR solution, parameters 
%     1. first roll up objective & constraint functions to be functions of states, controls, prices, paramaters (s c p par) 
%     2. second combine these to be functions of states, controls, prices, paramaters, or functions of solutions, parameters

%qq = @(s,c,par) 1; 
qq = @(s,c,par)  model.functions.q(c.k(s.a,s.k,s.z), c.k(c.a(s.a,s.k,s.z), c.k(s.a,s.k,s.z), s.z), par.delta, par.lambda); % not technically correct                                   
pipi_foc = @(s,c,p,par) model.functions.profit_foc(p.w,c.k(s.a, s.k, s.z), s.z, par.alpha, par.eta);
c  =  @(s,c,p,par) model.functions.c(p.w, p.r, s.a, s.k, s.z, c.a,              c.k,              par.alpha, par.delta, par.eta, par.lambda, par.theta);                                                                             
cf  = @(s,c,p,par) model.functions.c(p.w, p.r, s.a, s.k, s.z, c.a(s.a,s.k,s.z), c.k(s.a,s.k,s.z), par.alpha, par.delta, par.eta, par.lambda, par.theta);                                                                             
ccf = @(s,c,p,par) model.functions.c(p.w, p.r, c.a(s.a,s.k,s.z), c.k(s.a,s.k,s.z), s.z, c.a(c.a(s.a,s.k,s.z),c.k(s.a,s.k,s.z),s.z), c.k(c.a(s.a,s.k,s.z),c.k(s.a,s.k,s.z),s.z), par.alpha, par.delta, par.eta, par.lambda, par.theta);

model.conditions.b1 = @(states, controls, prices, parameters) model.functions.objective(c(states, controls, prices, parameters),paramters.gamma); % bellman term 1
model.conditions.b1f = @(states, controls, prices, parameters) model.functions.objective(cf(states, controls, prices)); % bellman term 1 

model.conditions.eulercombo.q1   = @(z,a,k,c,parameters,prices) model.functions.objective_foc(c(z,a,k),parameters.gamma).*(model.functions.profit_foc(prices.w, k, z, parameters.alpha, parameters.eta)+                     (1-parameters.delta)-prices.r);
model.conditions.eulercombo.q2   = @(z,a,k,c,parameters,prices) model.functions.objective_foc(c(z,a,k),parameters.gamma).*(model.functions.profit_foc(prices.w, k, z, parameters.alpha, parameters.eta)./parameters.lambda + (1-parameters.delta)-prices.r);
model.conditions.eulercombo.k.q1 = @(z,a,k,c,parameters,prices) model.functions.objective_foc(c(z,a,k),parameters.gamma).*(model.functions.profit_foc(prices.w, k, z, parameters.alpha, parameters.eta)+                     (1-parameters.delta));
model.conditions.eulercombo.a.q1 = @(z,a,k,c,parameters,prices) model.functions.objective_foc(c(z,a,k),parameters.gamma).*(prices.r);

model.conditions.rhs = @(z,a,k,c,parameters,prices)           model.functions.objective_foc(c(z,a,k),parameters.gamma).*parameters.beta*prices.r;

model.conditions.eeRHS.k = @(solution, parameters) parameters.beta*model.functions.objective_foc(ccf(solution.states.ndgrid, solution.controls.basis, solution.prices.values, parameters), parameters.gamma).*(qq(solution.states.ndgrid, solution.controls.basis, parameters) + pipi_foc(solution.states.ndgrid, solution.controls.basis, solution.prices.values, parameters));  % euler equation RHS  
model.conditions.eeRHS.a = @(solution, parameters) parameters.beta*model.functions.objective_foc(ccf(solution.states.ndgrid, solution.controls.basis, solution.prices.values, parameters), parameters.gamma).*solution.prices.values.r;   % euler equation RHS for 
model.conditions.ee.k = @(solution, parameters, expected_ee_RHS) (cf(solution.states.ndgrid, solution.controls.basis, solution.prices.values, parameters) + (parameters.delta.*parameters.lambda.*(1-parameters.delta)).*solution.states.ndgrid.k - model.functions.objective_foc_inv(expected_ee_RHS, parameters.gamma))./(parameters.delta.*parameters.lambda.*(1-parameters.delta));
model.conditions.ee.a = @(solution, parameters, expected_ee_RHS)  cf(solution.states.ndgrid, solution.controls.basis, solution.prices.values, parameters) + solution.states.ndgrid.a - model.functions.objective_foc_inv(expected_ee_RHS, parameters.gamma);


%% market clearing conditions
% next define market clearing conditions, to be used by *update_prices* to calculate excess demand, as functions of solution, parameters 
%     1. first roll up objective & constraint functions to be functions of solutions, paramaters (s c p par) 
%     2. second combine these to be functions of solutions, parameters

labor = @(solution, parameters) model.functions.labor(solution.prices.values.w, solution.states.ndgrid.k, solution.states.ndgrid.z, model.parameters.alpha, model.parameters.eta);

model.clearing.r = @(solution, parameters) dot(solution.states.stack(:,1) - parameters.theta*parameters.lambda*(1-parameters.delta)*solution.states.stack(:,2) - parameters.bmin, solution.distribution.values); 
model.clearing.w = @(solution, parameters) dot(reshape(labor(solution),[],1), solution.distribution.values) -1; 


%% controls
% final define functions to get equations for controls that aren't also endogenous states
%     1. first roll up objective & constraint functions to be functions of solutions, paramaters (s c p par) 
%     2. second combine these to be functions of solutions, parametersk

model.other.c = @(solution, parameters) c(solution.states.ndgrid, solution.controls.values, solution.prices.values, parameters);


end