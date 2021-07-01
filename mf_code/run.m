clear all; clc;

%cepath='/Applications/MATLAB_R2021a.app/toolbox/matlab/'; path([cepath 'cetools'],path);
 
guesses.prices.interval.w = [1 5];
guesses.prices.interval.r = [1 1.5];
guesses.prices.values.w = guesses.prices.interval.w(1); 
guesses.prices.values.r = guesses.prices.interval.r(1);  

[model, method, solution]  = solve(struct(),struct(), guesses);


%[model, method, solution]  = solve(huggett());



