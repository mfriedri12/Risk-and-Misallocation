%% setup solution 
%-------------------------------------------------------------------------%
% define all model parameters & model functions; ideally should be no
% model-specific info in the rest of the code 
%-------------------------------------------------------------------------%

function [model] = setup_model


%% parameters
% these are for use in the functions below 
model.parameters.A = 1;
model.parameters.alpha = 1/3; % capital share 
model.parameters.beta = .95; % discount  
model.parameters.gamma = 1.5; % risk aversion / IES 
model.parameters.delta = .05; % capital depreciation
model.parameters.eta = .85; % DRS at individual level 
model.parameters.lambda = .25; % what you get to keep  
model.parameters.theta = 0; % borrowing constraint capital doens't help lending
model.parameters.rho.z = .9; % productivity persistance 
model.parameters.sigma.z = .2; % productivity innovation variance
model.parameters.bmin = 10; % productivity innovation variance  % should be smaller than like 20 or so since that's roughly lifetime wage income 
model.parameters.L = 1; 

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
model.grid.min.z = -3;
model.grid.min.a = .001;
model.grid.min.k = .001;
model.grid.max.z = 3;
model.grid.max.a = 200;
model.grid.max.k = 250;


%% objective & constraint functions -- functions of primatives 
% objective 
model.functions.objective = @(c, gamma)  (c>0).*(c.^(1-gamma))./(1-gamma) + (c<=0).*-1e3; %hate CRRA  
model.functions.objective_foc = @(c, gamma) (c.^(-gamma)); 
model.functions.objective_foc_inv = @(argument, gamma) argument.^(1/(-gamma));

% the price function                                     
model.functions.q = @(k, kk, delta, lambda) (kk >= (1-delta)*k) * 1  + (kk <  (1-delta)*k) * lambda;  %

% the optimized labor functions 
model.functions.labor = @(w, k, z, alpha, eta) (w./ (( model.parameters.A .* exp(z).^(1-(1-alpha)*eta)).* (k.^(alpha*eta)).*  ((1-alpha)*eta))).^    (1/((1-alpha)*eta - 1));  

% the optimized profit functions 
model.functions.profit =     @(w, k, z, alpha, eta)                                 model.parameters.A .* exp(z).* (k.^((alpha*eta)/(1-(1-alpha)*eta))) .*  (1/w).^(((1-alpha)*eta)/(1-(1-alpha)*eta)) .*  ((1-alpha)*eta).^(((1-alpha)*eta )/(1-(1-alpha)*eta)) .*  (1-  ((1-alpha)*eta) );
model.functions.profit_foc = @(w, k, z, alpha, eta)     ((alpha*eta)/(1-(1-alpha)*eta)).* exp(z).* (k.^((eta-1)/(1-(1-alpha)*eta)))     .*  (1/w).^(((1-alpha)*eta)/(1-(1-alpha)*eta)) .*  ((1-alpha)*eta).^(((1-alpha)*eta )/(1-(1-alpha)*eta)) .*  (1-  ((1-alpha)*eta) );

% the constraint solved for c
model.functions.c = @(w, r, a, k, z, aa, kk, alpha, delta, eta, lambda, theta, bmin) exp(z).*w  + model.functions.profit(w, k, z, alpha, eta)  + r*(a - theta*lambda*(1-delta)*k - bmin)  - (aa - theta*lambda*(1-delta)*kk - bmin) - model.functions.q(k, kk, delta, lambda).*(kk - (1-delta)*k) ; 

% max a and k allowed for VFI algorithm -- we're assuming in teh max k that it's high q and a = 0 and relying on the piecewise CRRA to not get us into trouble 
model.functions.max.a = @(w, r, a, k, z, kk, alpha, delta, eta, lambda, theta, bmin)      exp(z).*w  + model.functions.profit(w, k, z, alpha, eta)  + r*(a - theta*lambda*(1-delta)*k - bmin) + theta*lambda*(1-delta)*kk + bmin - model.functions.q(k, kk, delta, lambda).*(kk - (1-delta)*k) ; 
model.functions.max.k = @(w, r, a, k, z,     alpha, delta, eta, lambda, theta, bmin)    ( exp(z).*w  + model.functions.profit(w, k, z, alpha, eta)  + r*(a - theta*lambda*(1-delta)*k - bmin)                                                           + 1.*(1-delta)*k)                  ./(1 - r*theta*lambda*(1-delta)) ; % a bit too loose exact but shouldn't matter as long as consumption isn't negative 


%% utilities for specific alogirithms -- functions of states, prices, parameters 
% for VFI, flip to make interpolation more linear for vtrick
model.utilities.vtrick_inv  = @(v, parameters) v.^(1-parameters.gamma); 
model.utilities.vtrick =      @(v, parameters) v.^(1/(1-parameters.gamma));

% for VFI  maximium a', k' allowable so that c not negative, 
model.utilities.max.a = @(states, prices, parameters) model.functions.max.a(prices.w, prices.r, states.a, states.k, states.z, states.kk, parameters.alpha, parameters.delta, parameters.eta, parameters.lambda, parameters.theta, parameters.bmin) ; 
model.utilities.max.k = @(states, prices, parameters) model.functions.max.k(prices.w, prices.r, states.a, states.k, states.z, parameters.alpha, parameters.delta, parameters.eta, parameters.lambda, parameters.theta, parameters.bmin) ; 

% optimal conditions  
model.conditions.b1 =  @(c1, c2, states, prices, parameters) model.functions.objective(model.functions.c(prices.w, prices.r, states.a, states.k, states.z, c1,c2, parameters.alpha, parameters.delta, parameters.eta, parameters.lambda, parameters.theta, parameters.bmin),parameters.gamma); % bellman term 1

% optimal conditions  
model.other.c = @(controls, states, prices, parameters) model.functions.c(prices.w, prices.r, states.a, states.k, states.z, controls.a,controls.k, parameters.alpha, parameters.delta, parameters.eta, parameters.lambda, parameters.theta, parameters.bmin); % bellman term 1


%% market clearing conditions
% next define market clearing conditions, to be used by *update_prices* to calculate excess demand, as functions of solution, parameters 
model.clearing.w = @(states, prices, distribution, parameters) parameters.L-dot(reshape(model.functions.labor(prices.w, states.k, states.z, parameters.alpha, parameters.eta),[],1), distribution); % negative means too much labor demand-- wage got to be higher. positive means too little labor demand -wage got to be lower. 
model.clearing.r = @(states, prices, distribution, parameters) dot(states.a - parameters.theta*parameters.lambda*(1-parameters.delta)*states.k - parameters.bmin, distribution); 


end