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
model.functions.c = @(r, b, z, bb) 2*z  + r*b  - bb; 


%% optimal conditions   
c  = @(s,c,p,par) model.functions.c(p.r, s.b, s.z, c);   
cf  = @(s,c,p,par) model.functions.c(p.r, s.b, s.z, c.b(s.a,s.z));                                                                             
ccf = @(s,c,p,par) model.functions.c(p.r, c.b(s.b,s.z), s.z, c.b(c.b(s.a,s.z),s.z));
 
model.conditions.max.b = @(solution, parameters) solution.states.ndgrid.z  + solution.prices.values.r*solution.states.ndgrid.b  ; 
model.conditions.b1 = @(states, controls, prices, parameters) model.functions.objective(c(states, controls, prices, parameters), parameters.gamma); % bellman term 1 

%% market clearing conditions
model.clearing.r = @(solution, parameters) dot(solution.states.stack(:,2),solution.distribution.values); 


%% controls
model.other.c = @(solution, parameters) c(solution.states.ndgrid, solution.controls.values.b, solution.prices.values);


end