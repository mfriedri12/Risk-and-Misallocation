%% maximize  
%-------------------------------------------------------------------------%
% Now it still only works over 2 variables 
%-------------------------------------------------------------------------%

function [argmax, maximum] = maximize(method, objective, search_min, search_max, additional_arguments)

switch method.maximize.algorithm 
    
    
    case {'golden search'} % hill climibing; uses bellman; really only works for models with one endogenous state variable right now 
       
      [argmax, maximum] = golden_mrf(objective, method.maximize.threshold, search_min, search_max, additional_arguments);  
             
end

%need ot check c not negative bc my max isn't doing it properly 




% objective should be defined as V(aa,kk,z,a,k) 

% I need to make model.conditions.max.a and  model.conditions.max.k a function of ndgridplus -- aka
% ndgrid plus another dimension for k', and I need to add ndgridplus to my
% list of states (will that be too much memory?) 
% There must in general be a way to reduce the amount of memory used I
% don't know it though I don't understand how to ...

% could I do it as a dot, at least for the state variables? possibly. I
% will at the end of the day need zakk for the min and max and to do the
% search over that space so I can knock on down to a zak 

% objective needs to a be a function of 5 variables -- the three states
% plus the two other variables. It's ok, we can do that without too much
% effort since it's the addition of the value function which is an interpolant but a function
% of only three and the utility function which is explict and a function of
% all 5
 

%let's say I remake the objective every time so that it's only a function
%of the things we're maximizing over... that makes sense eh 
%then the only states argument I have to enter is for the other variable if
%I'm maximizing over several things 
% oh no that won't work becasue of the grid thing 

% idea is that we first find a' or model.variables.endogenous(1) for every
% k' meaning that we're searching on a zxaxkxk grid 

% I can sub out that for loop for a proper paralellizaiton 

%but how do I flexibly change the grid ? oooh I dot it I get it. 
 
