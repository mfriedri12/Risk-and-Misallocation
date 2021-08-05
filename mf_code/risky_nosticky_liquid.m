%-------------------------------------------------------------------------%
% This is like the BLPT model 
%-------------------------------------------------------------------------%
function [model] = risky_nosticky_liquid()

%% parameters
model.parameters.rb =1.02; 
model.parameters.qk =1; 
model.parameters.rho.z = 0; % productivity persistance 
model.parameters.sigma.z = .15; % productivity innovation variance
model.parameters.tau = .3; 

%% variables
model.variables.exogenous = ['z']; %  ["z: productivity"]; exogenous states; assume these are AR(1)s for now 
model.variables.endogenous = ['b','k']; % ["a: liquid assets" "k: capital"];  % engogenous states = intersection of states & controls. the key variables 
model.variables.other = ['c']; %[ "c: consumption"]
model.variables.prices = ['r']; % ["w: wages" "r: liquid asset return"];
model.variables.states = ['z', 'b','k']; % [model.variables.exogenous, model.variables.endogenous];
model.variables.controls = ['b','k','c']; %[model.variables.endogenous, model.variables.other]; 

%% grid 
model.grid.min.z = -5; %0.01;
model.grid.min.b = 0.001;  %0.01;
model.grid.min.k = 0.001; %0.01;
model.grid.max.z = 5;
model.grid.max.b = 300;
model.grid.max.k = 200;

%% objective & constraint problem  -- functions of primatives 
model.functions.objective = @(c, gamma)  (c>0).*(c.^(1-gamma))./(1-gamma) + (c<=0).*-10e5; %hate CRRA  

model.functions.r =     @(K, L, alpha) alpha.*(K.^(alpha-1)).*(L.^(1-alpha)); 
model.functions.K =     @(r, L, alpha) L./((r./alpha).^(1./(1-alpha))); % (r./(alpha.*(L.^(1-alpha)))).^(1./(alpha-1)); 
model.functions.w =     @(r, L, alpha) (1-alpha).*(model.functions.K(r,L,alpha).^(alpha)).*(L.^(-alpha)); 

model.functions.c =     @(rb, w, r, qk, b, k, z, bb, kk, delta, tau) (1-tau).*exp(z).*w + exp(z).*k.*(r + qk.*(1-delta))  + rb.*b  - bb - qk.*kk;

model.functions.max.b = @(rb, w, r, qk, b, k, z, kk, delta, tau) (1-tau).*exp(z).*w + exp(z).*k.*(r + qk.*(1-delta))  + rb.*b   - qk.*kk; 
model.functions.max.k = @(rb, w, r, qk, b, k, z, delta, tau)    ((1-tau).*exp(z).*w + exp(z).*k.*(r + qk.*(1-delta))  + rb.*b)./ qk; 

%% update_controls functions  -- functions of states, prices, parameters 
model.utilities.max.b = @(states, prices, parameters) model.functions.max.b(parameters.rb, model.functions.w(prices.r,parameters.L,parameters.alpha), prices.r, parameters.qk, states.b, states.k, states.z, states.kk, parameters.delta, parameters.tau);
model.utilities.max.k = @(states, prices, parameters) model.functions.max.k(parameters.rb, model.functions.w(prices.r,parameters.L,parameters.alpha), prices.r, parameters.qk, states.b, states.k, states.z, parameters.delta, parameters.tau);

model.conditions.b1 =   @(c1, c2, states, prices, parameters) model.functions.objective(model.functions.c(parameters.rb, model.functions.w(prices.r,parameters.L,parameters.alpha), prices.r, parameters.qk, states.b, states.k, states.z, c1, c2, parameters.delta, parameters.tau),parameters.gamma); 

model.other.c = @(controls, states, prices, parameters) model.functions.c(parameters.rb, model.functions.w(prices.r,parameters.L,parameters.alpha), prices.r, parameters.qk, states.b, states.k, states.z, controls.b, controls.k, parameters.delta, parameters.tau); 

%% market clearing conditions -- functions of states, prices, distribution, parameters 
model.clearing.r = @(states, prices, distribution, parameters)  dot(states.k.*exp(states.z),distribution) - model.functions.K(prices.r, parameters.L, parameters.alpha); 

end