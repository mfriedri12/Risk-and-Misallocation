%% update prices 
%-------------------------------------------------------------------------%
% given distribution, update prices
%-------------------------------------------------------------------------%

function [solution] = update_prices(model, method, solution, price, solve_iteration)


%% caluculate excess demand 
% see how much guessed prices vary from solution prices 

excess = model.clearing.(price)(solution, model.parameters);
fprintf('EXCESS: %f', excess); 


%% check interval  

if solve_iteration==0
   
   solution.prices.values.(price) = solution.prices.interval.(price)(2);
   
   fprintf('\n ITERATION %f PRICE: %f \n', solve_iteration, solution.prices.values.(price));
   
   tic; disp('1. update solution.controls given solution.prices');       [solution] = update_controls(model, method, solution); toc; 
   tic; disp('2. update solution.distribution given solution.controls'); [solution] = update_distribution(model, method, solution); toc; 
   
   excess2 = model.clearing.(price)(solution, model.parameters);
   fprintf('EXCESS: %f \n', excess2); 
   
   if sign(excess/excess2) == 1 
        disp('WARNING: interval not valid! Paused.'); 
        pause
   else
       solution.prices.values.(price) = solution.prices.interval.(price)(2)/method.update_prices.interval;
   end
       
else 
    

%% update price

switch method.update_prices.algorithm 
   
    case 'brent-dekker'
        
        % if fist iteration 
        if solve_iteration == 1
            solution.prices.algorithm.lastexcess.(price) = excess; % ensures that bisection will be used on the first go round
            solution.prices.algorithm.lastlastprice.(price) = 0; % won't be relevant until refilled 
            solution.prices.algorithm.lastupdate.(price) = 'bisection'; % won't be relevant until refilled 
            if excess > 0
                solution.prices.algorithm.lastprice.(price) = max(solution.prices.interval.(price));  
            elseif excess < 0      
               solution.prices.algorithm.lastprice.(price) = min(solution.prices.interval.(price));  
            end
        end  
        
        % using excess demand, update contrapoint and interval 
        if excess > 0
            contrapoint = min(solution.prices.interval.(price));  % demand too large so need to lower price 
            solution.prices.interval.(price) = [contrapoint solution.prices.values.(price)];
        elseif excess < 0      
            contrapoint = max(solution.prices.interval.(price));  % demand too small so need to increase price 
            solution.prices.interval.(price) = [solution.prices.values.(price) contrapoint];
        end   
        
        % calculate update options 
        secant = solution.prices.values.(price)-excess*((solution.prices.values.(price)-solution.prices.algorithm.lastprice.(price))/(excess-solution.prices.algorithm.lastexcess.(price))); % secant_guess;
        bisection = (contrapoint + solution.prices.values.(price))/2 ;  % bisection_guess
        
        % calculate update conditions  
        dekker_condition_1 = (excess > 0 && secant < bisection && secant > solution.prices.values.(price)) || (excess < 0 && secant > bisection && secant < solution.prices.values.(price)); % secant_guess betwen midpoint and guess 
        dekker_condition_2 = abs(excess - solution.prices.algorithm.lastexcess.(price))>0;
        switch solution.prices.algorithm.lastupdate.(price) 
            case 'bisection'
                brent_condition_1 = abs(solution.prices.values.(price) - solution.prices.algorithm.lastprice.(price)) > method.update_prices.brent_threshold;
                brent_condition_2 = abs(secant - solution.prices.values.(price)) < .5*abs(solution.prices.values.(price) - solution.prices.algorithm.lastprice.(price));
            case 'secant'
                brent_condition_1 = abs(solution.prices.algorithm.lastprice.(price) - solution.prices.algorithm.lastlastprice.(price)) > method.update_prices.brent_threshold;
                brent_condition_2 = abs(secant - solution.prices.values.(price)) < .5*abs(solution.prices.algorithm.lastprice.(price) - solution.prices.algorithm.lastlastprice.(price));
        end 
           
        % save info for next iteration 
        solution.prices.algorithm.lastlastprice.(price) = solution.prices.algorithm.lastprice.(price);
        solution.prices.algorithm.lastprice.(price) = solution.prices.values.(price);
        solution.prices.algorithm.lastexcess.(price) = excess; 
        
        % update price
        if dekker_condition_1 && dekker_condition_2 && brent_condition_1 && brent_condition_2
            solution.prices.values.(price) = secant;
            solution.prices.algorithm.lastupdate.(price) = 'secant';
            pause;
        else
            solution.prices.values.(price) = bisection;
            solution.prices.algorithm.lastupdate.(price) = 'bisection';
        end
        fprintf(' UPDATE: %s \n',  solution.prices.algorithm.lastupdate.(price)); 
       
    case 'bisection'
    case 'newton'
end

solution.accuracy.prices.(price) = abs(solution.prices.values.(price)- solution.prices.algorithm.lastprice.(price));


end


end