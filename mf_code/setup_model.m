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
model.parameters.theta = 0; % borrowing constraint capital doens't help lending
model.parameters.rho = .9; % productivity persistance 
model.parameters.sigma = .15; % productivity innovation variance
model.parameters.bmin = 5; % productivity innovation variance  % should be smaller than like 20 or so since that's roughly lifetime wage income 


%% variables  
% this just records the basic dimensions of the model

model.variables.exogenous = ['z']; %  ["z: productivity"]; exogenous states; assume these are AR(1)s for now 
model.variables.endogenous = ['a','k']; % ["a: liquid assets" "k: capital"];  % engogenous states = intersection of states & controls. the key variables 
model.variables.other = ['c']; %[ "c: consumption"]
model.variables.prices = ['r','w']; % ["w: wages" "r: liquid asset return"];
model.variables.states = ['z', 'a','k']; % [model.variables.exogenous, model.variables.endogenous];
model.variables.controls = ['a','k','c']; %[model.variables.endogenous, model.variables.other]; 


%% grid
% this defines endogneous states grid. States are listed in order
% noted above. Agents are not allowed to choose a state variable outside of this
% grid so that's why it's part of model 

model.grid.min = [0.0001 0.0001];  % borrowing constraints -- not zero because sometimes there are functions where the marginal return is infinite
model.grid.max = [100 100]; % this is not so important but we're still gonna call it model specific 


%% exogenous state processes
%used by setup solution to get grid for these and also by *integrate* to calculate expected value

model.exogenous.process.z = @(solution, shocks, parameters) (solution.states.grid.z.^parameters.rho).*exp(shocks); %lognormal 
model.exogenous.grid.z = @(gridsize, parameters) exp(tauchen(0,parameters.rho,parameters.sigma,gridsize)*(parameters.sigma/sqrt(1-parameters.rho^2)));  % took this from other code, it's "correcting for human capital" exp(tauchen(0,parameters.rho,parameters.sigma,gridsize)*(parameters.sigma/sqrt(1-parameters.rho^2))); %   
model.exogenous.transition.z = @(gridsize, parameters) tauchen(0,parameters.rho, parameters.sigma, gridsize); % second argument .. might be wrong 


%% objective & constraint functions 
% first define subfunctions naively -- as a function of model-specific states, parameters

% objective 
model.functions.objective = @(c, gamma) (c>0).*(c.^(1-gamma))./(1-gamma) + (c<0).*-10e10; %hate CRRA  
model.functions.objective_foc = @(c, gamma) (c.^(-gamma)); 
model.functions.objective_foc_inv = @(argument, gamma) argument.^(1/(-gamma));

% the price function                                     
model.functions.q = @(k, kk, delta, lambda) (kk >= (1-delta)*k) * 1  + (kk <  (1-delta)*k) * (1-lambda);

% the optimized labor functions 
model.functions.labor = @(w, k, z, alpha, eta) (w./ ((z.^(1-(1-alpha)*eta)).* (k.^(alpha*eta)).*  ((1-alpha)*eta))).^    (1/((1-alpha)*eta - 1));  

% the optimized profit functions 
model.functions.profit =     @(w, k, z, alpha, eta)                                       z.* (k.^((alpha*eta)/(1-(1-alpha)*eta))) .*  (1/w).^(((1-alpha)*eta)/(1-(1-alpha)*eta)) .*  ((1-alpha)*eta).^(((1-alpha)*eta )/(1-(1-alpha)*eta)) .*  (1-  ((1-alpha)*eta) );
model.functions.profit_foc = @(w, k, z, alpha, eta)     ((alpha*eta)/(1-(1-alpha)*eta)).* z.* (k.^((eta-1)/(1-(1-alpha)*eta)))     .*  (1/w).^(((1-alpha)*eta)/(1-(1-alpha)*eta)) .*  ((1-alpha)*eta).^(((1-alpha)*eta )/(1-(1-alpha)*eta)) .*  (1-  ((1-alpha)*eta) );

% the constraint solved for c
model.functions.c = @(w, r, a, k, z, aa, kk, alpha, delta, eta, lambda, theta, bmin) w  + model.functions.profit(w, k, z, alpha, eta)  + r*(a - theta*lambda*(1-delta)*k - bmin)  - (aa - theta*lambda*(1-delta)*kk - bmin) - model.functions.q(k, kk, delta, lambda).*(kk - (1-delta)*k) ; 

% max a and k allowed for VFI algorithm -- we're assuming in teh max k that it's high q and a = 0 and relying on the piecewise CRRA to not get us into trouble 
model.functions.max.a = @(w, r, a, k, z, kk, alpha, delta, eta, lambda, theta, bmin)     w  + model.functions.profit(w, k, z, alpha, eta)  + r*(a - theta*lambda*(1-delta)*k - bmin) + theta*lambda*(1-delta)*kk + bmin - model.functions.q(k, kk, delta, lambda).*(kk - (1-delta)*k) ; 
model.functions.max.k = @(w, r, a, k, z,     alpha, delta, eta, lambda, theta, bmin)    (w  + model.functions.profit(w, k, z, alpha, eta)  + r*(a - theta*lambda*(1-delta)*k - bmin)                                                           + 1.*(1-delta)*k)                  ./(1 - r*theta*lambda*(1-delta)) ; % a bit too loose exact but shouldn't matter as long as consumption isn't negative 


%% utilities for specific alogirithms

% for VFI, flip to make interpolation more linear for vtrick
model.utilities.vtrick_inv  = @(v, parameters) v.^(1-parameters.gamma); 
model.utilities.vtrick =      @(v, parameters) v.^(1/(1-parameters.gamma));

% for VFI  maximium a', k' allowable so that c not negative, 
maxa = @(s, p, par, kk) model.functions.max.a(p.w, p.r, s.a, s.k, s.z, s.kk, par.alpha, par.delta, par.eta, par.lambda, par.theta, par.bmin) ; 
maxk = @(s, p, par) model.functions.max.k(p.w, p.r, s.a, s.k, s.z, par.alpha, par.delta, par.eta, par.lambda, par.theta, par.bmin) ; 
model.utilities.max.a = @(solution, parameters) maxa(solution.states.ndgridplus, solution.prices.values, parameters) ; 
model.utilities.max.k = @(solution, parameters) maxk(solution.states.ndgrid, solution.prices.values, parameters) ; 


%% optimal conditions    
% next define optimal conditions, which we will need for update_controls, as functions of states, controls, prices, parameters OR solution, parameters 
%     1. first roll up objective & constraint functions to be functions of states, controls, prices, paramaters (s c p par) 
%     2. second combine these to be functions of states, controls, prices, paramaters, or functions of solutions, parameters

% for VFI 
cak = @(s,a,k,p,par) model.functions.c(p.w, p.r, s.a, s.k, s.z, a,                k,              par.alpha, par.delta, par.eta, par.lambda, par.theta, par.bmin);                                                                             
model.conditions.b1 =  @(c1, c2, states, prices, parameters) model.functions.objective(cak(states, c1, c2,  prices, parameters),parameters.gamma); % bellman term 1


%% market clearing conditions
% next define market clearing conditions, to be used by *update_prices* to calculate excess demand, as functions of solution, parameters 
%     1. first roll up objective & constraint functions to be functions of solutions, paramaters (s c p par) 
%     2. second combine these to be functions of solutions, parameters

labor = @(solution, parameters) model.functions.labor(solution.prices.values.w, solution.states.ndgrid.k, solution.states.ndgrid.z, parameters.alpha, parameters.eta);
model.clearing.w = @(solution, parameters) 1-dot(reshape(labor(solution, parameters),[],1), solution.distribution.values); % negative means too much labor demand-- wage got to be higher. positive means too little labor demand -wage got to be lower. 
model.clearing.r = @(solution, parameters) dot(solution.states.stack(:,2) - parameters.theta*parameters.lambda*(1-parameters.delta)*solution.states.stack(:,3) - parameters.bmin, solution.distribution.values); 


%% controls
% final define functions to get equations for controls that aren't also endogenous states
%     1. first roll up objective & constraint functions to be functions of solutions, paramaters (s c p par) 
%     2. second combine these to be functions of solutions, parametersk

c  =  @(s,c,p,par) model.functions.c(p.w, p.r, s.a, s.k, s.z, c.a,              c.k,              par.alpha, par.delta, par.eta, par.lambda, par.theta, par.bmin);                                                                             
model.other.c = @(solution, parameters) c(solution.states.ndgrid, solution.controls.values, solution.prices.values, parameters);


end