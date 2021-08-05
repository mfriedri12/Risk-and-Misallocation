%% huggett
%-------------------------------------------------------------------------%
% This is a huggett calibration 
%-------------------------------------------------------------------------%

function [model, method] = huggett()

%% parameters 
model.parameters.sigma.z = .1; 

%% variables 
model.variables.exogenous = ['z']; %  ["z: productivity"]; exogenous states 
model.variables.endogenous = ['b']; % ["a: liquid assets" "k: capital"];  % engogenous states = intersection of states & controls. the key variables 
model.variables.other =['c'];
model.variables.prices = ['r']; % ["w: wages" "r: liquid asset return"];
model.variables.states = ['z', 'b'];
model.variables.controls = ['b', 'c'];

%% grid 
%limits for endogenous states 
model.grid.min.z = -.6882; 
model.grid.min.b = -3;  % borrowing constraints
model.grid.max.z = .6882; 
model.grid.max.b = 150;

%% objective & constraint problem  
model.functions.objective = @(c, gamma) (c.^(1-gamma))./(1-gamma); 
model.functions.c = @(bb, z, b, r) exp(z)  + r*b  - bb; 

%% setup for value function maximization 
model.utilities.max.b = @(states, prices, parameters) states.z  + prices.r*states.b; % maximum b possible so that c is not negative    
model.conditions.b1 =   @(control, states, prices, parameters) model.functions.objective(model.functions.c(control, states.z, states.b, prices.r), parameters.gamma); % bellman term 1 

%% market clearing conditions
model.clearing.r = @(states, prices, distribution, parameters) dot(states.b,distribution); 

%% controls
model.other.c = @(controls, states, prices, parameters) model.functions.c(controls.b, states.z, states.b, prices.r);

%% other 
method.update_prices.interval = 1.04;

end