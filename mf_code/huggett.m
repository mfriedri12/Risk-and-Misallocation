%% huggett
%-------------------------------------------------------------------------%
% This is a huggett calibration 
%-------------------------------------------------------------------------%

function [model] = huggett()


%% variables 
model.variables.exogenous = ['z']; %  ["z: productivity"]; exogenous states 
model.variables.endogenous = ['b']; % ["a: liquid assets" "k: capital"];  % engogenous states = intersection of states & controls. the key variables 
model.variables.other =['c'];
model.variables.prices = ['r']; % ["w: wages" "r: liquid asset return"];
model.variables.states = ['z', 'b'];
model.variables.controls = ['b', 'c'];


%% grid 
%limits for endogenous states 
model.grid.min = [-5];  % borrowing constraints
model.grid.max = [100]; % this is not so important but we're still gonna call it model specific 


%% objective & constraint problem  
model.functions.objective = @(c, gamma) (c.^(1-gamma))./(1-gamma); 
model.functions.c = @(r, b, z, bb) z  + r*b  - bb; 


%% utilities
model.utilities.max.b = @(solution, parameters) solution.states.ndgrid.z  + solution.prices.values.r*solution.states.ndgrid.b; % maximum b possible so that c is not negative  


%% optimal conditions   
c  = @(s,b,p,par) model.functions.c(p.r, s.b, s.z, b);   
model.conditions.b1 = @(c1, states, prices, parameters) model.functions.objective(c(states, c1, prices, parameters), parameters.gamma); % bellman term 1 


%% market clearing conditions
model.clearing.r = @(solution, parameters) dot(solution.states.stack(:,2),solution.distribution.values); 


%% controls
model.other.c = @(solution, parameters) c(solution.states.ndgrid, solution.controls.values.b, solution.prices.values, parameters);


end