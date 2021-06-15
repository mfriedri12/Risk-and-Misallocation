%% update prices 
%-------------------------------------------------------------------------%
% given distribution, update prices
%-------------------------------------------------------------------------%

function [solution] = update_prices(model, method, solution, i, solve_iteration)


%% caluculate excess demand 
% see how much guessed prices vary from solution prices 

excess = model.clearing.(model.variables.prices(i))(solution, model.parameters);
fprintf('EXCESS (%s): %f',model.variables.prices(i), excess); 

if length(model.variables.prices)==2
   other_excess = model.clearing.(model.variables.prices(2/i))(solution, model.parameters);
   fprintf(', EXCESS (%s): %f',model.variables.prices(2/i), other_excess); 
end


%% check interval  

if solve_iteration==1 
    
    solution.prices.algorithm.lastlastprice.(model.variables.prices(i)) = 0; % won't be relevant until refilled 
    solution.prices.algorithm.lastupdate.(model.variables.prices(i)) = 'bisection'; % won't be relevant until refilled   
      
    
   if method.update_prices.check_interval

       solution.prices.values.(model.variables.prices(i)) = solution.prices.interval.(model.variables.prices(i))(2);

       fprintf('\n\nITERATION %f PRICE: %f \n', solve_iteration, solution.prices.values.(model.variables.prices(i)));

       tic; disp('1. update solution.controls given solution.prices');       [solution] = update_controls(model, method, solution); toc; 
       tic; disp('2. update solution.distribution given solution.controls'); [solution] = update_distribution(model, method, solution); toc; 

       excess2 = model.clearing.(model.variables.prices(i))(solution, model.parameters);
       fprintf('EXCESS (%s): %f',model.variables.prices(i), excess2); 
       
       if length(model.variables.prices)==2
           fprintf(', EXCESS (%s): %f',model.variables.prices(2/i), model.clearing.(model.variables.prices(2/i))(solution, model.parameters)); 
       end


       if sign(excess/excess2) == 1 
            fprintf('\nWARNING: interval not valid! Paused.'); 
            pause
       else
            solution.prices.algorithm.lastexcess.(model.variables.prices(i)) = excess2; % ensures that bisection will be used on the first go round
            solution.prices.algorithm.lastprice.(model.variables.prices(i)) = solution.prices.interval.(model.variables.prices(i))(2);
            solution.prices.values.(model.variables.prices(i)) = solution.prices.interval.(model.variables.prices(i))(1); % reset so rest of algorithm works 
       end
       
   else
       
       solution.prices.algorithm.lastexcess.(model.variables.prices(i)) = excess; % ensures that bisection will be used on the first go round
       if excess > 0
           solution.prices.algorithm.lastprice.(model.variables.prices(i)) = max(solution.prices.interval.(model.variables.prices(i)));  
       elseif excess < 0      
           solution.prices.algorithm.lastprice.(model.variables.prices(i)) = min(solution.prices.interval.(model.variables.prices(i)));  
       end
    
   end    
       
end 
    

%% update price

switch method.update_prices.algorithm 
   
    case 'brent-dekker'
         
        % using excess demand, update contrapoint and interval 
        if excess > 0
            contrapoint = min(solution.prices.interval.(model.variables.prices(i)));  % demand too large so need to lower price 
            solution.prices.interval.(model.variables.prices(i)) = [contrapoint solution.prices.values.(model.variables.prices(i))];
        elseif excess < 0      
            contrapoint = max(solution.prices.interval.(model.variables.prices(i)));  % demand too small so need to increase price 
            solution.prices.interval.(model.variables.prices(i)) = [solution.prices.values.(model.variables.prices(i)) contrapoint];
        end   
        
        % calculate update options 
        secant = solution.prices.values.(model.variables.prices(i))-excess*((solution.prices.values.(model.variables.prices(i))-solution.prices.algorithm.lastprice.(model.variables.prices(i)))/(excess-solution.prices.algorithm.lastexcess.(model.variables.prices(i)))); % secant_guess;
        bisection = (contrapoint + solution.prices.values.(model.variables.prices(i)))/2 ;  % bisection_guess
        
        % calculate update conditions  
        dekker_condition_1 = (excess > 0 && secant < bisection && secant > solution.prices.values.(model.variables.prices(i))) || (excess < 0 && secant > bisection && secant < solution.prices.values.(model.variables.prices(i))); % secant_guess betwen midpoint and guess 
        dekker_condition_2 = abs(excess - solution.prices.algorithm.lastexcess.(model.variables.prices(i)))>0;
        switch solution.prices.algorithm.lastupdate.(model.variables.prices(i)) 
            case 'bisection'
                brent_condition_1 = abs(solution.prices.values.(model.variables.prices(i)) - solution.prices.algorithm.lastprice.(model.variables.prices(i))) > method.update_prices.brent_threshold;
                brent_condition_2 = abs(secant - solution.prices.values.(model.variables.prices(i))) < .5*abs(solution.prices.values.(model.variables.prices(i)) - solution.prices.algorithm.lastprice.(model.variables.prices(i)));
            case 'secant'
                brent_condition_1 = abs(solution.prices.algorithm.lastprice.(model.variables.prices(i)) - solution.prices.algorithm.lastlastprice.(model.variables.prices(i))) > method.update_prices.brent_threshold;
                brent_condition_2 = abs(secant - solution.prices.values.(model.variables.prices(i))) < .5*abs(solution.prices.algorithm.lastprice.(model.variables.prices(i)) - solution.prices.algorithm.lastlastprice.(model.variables.prices(i)));
        end 
           
        % save info for next iteration 
        solution.prices.algorithm.lastlastprice.(model.variables.prices(i)) = solution.prices.algorithm.lastprice.(model.variables.prices(i));
        solution.prices.algorithm.lastprice.(model.variables.prices(i)) = solution.prices.values.(model.variables.prices(i));
        solution.prices.algorithm.lastexcess.(model.variables.prices(i)) = excess; 
        
        % update price
        if dekker_condition_1 && dekker_condition_2 && brent_condition_1 && brent_condition_2
            solution.prices.values.(model.variables.prices(i)) = secant;
            solution.prices.algorithm.lastupdate.(model.variables.prices(i)) = 'secant';
            pause;
        else
            solution.prices.values.(model.variables.prices(i)) = bisection;
            solution.prices.algorithm.lastupdate.(model.variables.prices(i)) = 'bisection';
        end
        fprintf(', UPDATE: %s \n\n',  solution.prices.algorithm.lastupdate.(model.variables.prices(i))); 
       
    case 'bisection'
    case 'newton'
end

solution.accuracy.prices.(model.variables.prices(i)) = abs(solution.prices.values.(model.variables.prices(i))- solution.prices.algorithm.lastprice.(model.variables.prices(i)));
    

solution.prices.algorithm.allexcess(solve_iteration,:) = [excess other_excess];
disp(solution.prices.algorithm.allexcess);



end


