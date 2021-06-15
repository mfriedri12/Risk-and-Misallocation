clear all; clc;
 
guesses.prices.values.w = .32; 
guesses.prices.values.r = 1.005; 
guesses.prices.interval.w = [.32 .33];
guesses.prices.interval.r = [1.005 1.02];
 
[model, method, solution]  = solve(struct(),struct(), guesses);

%[model, method, solution]  = solve(huggett());


