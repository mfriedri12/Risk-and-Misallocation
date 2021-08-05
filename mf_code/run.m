%% Prelims 
clear all; clc;
cepath='/Applications/MATLAB_R2021a.app/toolbox/matlab/'; path([cepath 'cetools'],path); %for compecon toolkit


%% Huggett model 
%[huggettmodel, huggettmethod] =  huggett();
%[model, method, solution]  = solve(huggettmodel, huggettmethod);


%% Risky nosticky liquid model 
% guesses.prices.interval.r = [.05 .2];
% guesses.prices.values.r = guesses.prices.interval.r(1); 
% [model, method, solution]  = solve(risky_nosticky_liquid(),struct(), guesses);


%% Twostates  model 
guesses.prices.interval.r = [.05 .2];
guesses.prices.values.r = guesses.prices.interval.r(1); 
[model, method, solution]  = solve(twostates(),struct(), guesses);


%% Risk-and-Misallocation model 
%guesses.prices.interval.w  = [.5 1.4]; %[0.775000 1]; %[.6 1.4];
%guesses.prices.interval.r = [.98 1.02]; %[1.011250 1.2];  %[.98 1.02];
%guesses.prices.values.w = guesses.prices.interval.w(1); 
%guesses.prices.values.r = guesses.prices.interval.r(1);  
%[model, method, solution]  = solve(struct(),struct(), guesses);

