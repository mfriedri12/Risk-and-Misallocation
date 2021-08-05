%% update prices 
%-------------------------------------------------------------------------%
% given distribution, update prices
%-------------------------------------------------------------------------%

function [solution] = update_prices(model, method, solution, solve_iteration)


%% caluculate excess demand 
% see how much guessed prices vary from solution prices 

for i=1:length(model.variables.prices)
    solution.prices.algorithm.prices.(model.variables.prices(i))(solve_iteration) = solution.prices.values.(model.variables.prices(i));         
    solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration) = model.clearing.(model.variables.prices(i))(solution.states.stack, solution.prices.values, solution.distribution.values, model.parameters);
    fprintf('EXCESS (%s): %f ',model.variables.prices(i), solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration)); 
end
fprintf('\n');


%% check interval  

if solve_iteration==1 && method.update_prices.check_interval         
      
       % calculate excess value at all points in the interval (or grid, if it's 2 prices) abcd
       % note this code assumes that you're first run was at min-min (which is how it's coded in setup_solution if this option is true)
       for i=1:(length(model.variables.prices)*2-1)  
           % for i=0 aka first run should be min min (min first price & min second price) excess should be - -
           if i==1 % max max (max first price & max second price) excess should be + + 
               solution.prices.values.(model.variables.prices(1)) = solution.prices.interval.(model.variables.prices(1))(2);
               if length(model.variables.prices)>1 
                    solution.prices.values.(model.variables.prices(2)) = solution.prices.interval.(model.variables.prices(2))(2);
               end
           elseif i==2 % min max (min first price & max second price) excess should be - + 
               solution.prices.values.(model.variables.prices(1)) = solution.prices.interval.(model.variables.prices(1))(1);
           elseif i==3 % max min (max first price & min second price) excess should be + -
               solution.prices.values.(model.variables.prices(1)) = solution.prices.interval.(model.variables.prices(1))(2);
               solution.prices.values.(model.variables.prices(2)) = solution.prices.interval.(model.variables.prices(2))(1);
           end
   
           fprintf('\nITERATION %d (CHECK GRID %d) ', solve_iteration, i); 
           for j=1:length(model.variables.prices); fprintf('PRICE (%s): %f ', model.variables.prices(j), solution.prices.values.(model.variables.prices(j))); end 
           fprintf('\n');
            
           disp('1. update solution.controls given solution.prices');       [solution] = update_controls(model, method, solution); toc; 
           disp('2. update solution.distribution given solution.controls'); [solution] = update_distribution(model, method, solution); toc; 

           for j=1:length(model.variables.prices)    
               solution.prices.algorithm.prices.(model.variables.prices(j))(solve_iteration+i) = solution.prices.values.(model.variables.prices(j));
               solution.prices.algorithm.excess.(model.variables.prices(j))(solve_iteration+i) = model.clearing.(model.variables.prices(j))(solution.states.stack, solution.prices.values, solution.distribution.values, model.parameters);
               fprintf('EXCESS (%s): %f ',model.variables.prices(j), solution.prices.algorithm.excess.(model.variables.prices(j))(solve_iteration+i)); 
           end     
           fprintf('\n');
       end
       
       % check that excesses are appropriate signs 
       if (length(model.variables.prices)==1 && ...
          all((sign(solution.prices.algorithm.excess.(model.variables.prices(1)))==[-1 1])))  ...
          ||...
          (length(model.variables.prices)==2 && ...
          all((sign(solution.prices.algorithm.excess.(model.variables.prices(1)))==[-1 1 -1 1])) && ...
          all((sign(solution.prices.algorithm.excess.(model.variables.prices(2)))==[-1 1 1 -1])))
          fprintf('Interval OK! \n');
          solve_iteration = 2; %solve iterations 3 & 4 will be overwritten if it's 2 prices, but that's ok 
          if length(model.variables.prices)==2
            solution.prices.values.(model.variables.prices(2)) = solution.prices.algorithm.prices.(model.variables.prices(2))(2);
          end
       else 
          fprintf('\nWARNING: interval not valid! Paused.'); 
          pause 
       end
       
end 
    

%% update price

for i=1:length(model.variables.prices)

    switch method.update_prices.algorithm 

        case 'brent-dekker'

            % using excess demand, update contrapoint and interval 
            if solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration) > 0
                contrapoint = min(solution.prices.interval.(model.variables.prices(i)));  % demand too large so need to lower price 
                solution.prices.interval.(model.variables.prices(i)) = [contrapoint solution.prices.values.(model.variables.prices(i))];
            else  
                contrapoint = max(solution.prices.interval.(model.variables.prices(i)));  % demand too small so need to increase price 
                solution.prices.interval.(model.variables.prices(i)) = [solution.prices.values.(model.variables.prices(i)) contrapoint];
            end   

            % calculate update options 
            secant = solution.prices.values.(model.variables.prices(i))-solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration)*...
                    (solution.prices.values.(model.variables.prices(i))-solution.prices.algorithm.prices.(model.variables.prices(i))(solve_iteration-1))...
                   /(solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration)-solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration-1)); % secant_guess;
            bisection = (contrapoint + solution.prices.values.(model.variables.prices(i)))/2 ;  % bisection_guess

            % calculate update conditions  
            dekker_condition_1 = (solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration) > 0 && secant < bisection && secant > solution.prices.values.(model.variables.prices(i))) || (solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration) < 0 && secant > bisection && secant < solution.prices.values.(model.variables.prices(i))); % secant_guess betwen midpoint and guess 
            dekker_condition_2 = abs(solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration) - solution.prices.algorithm.excess.(model.variables.prices(i))(solve_iteration-1))>0;
            switch solution.prices.algorithm.update.(model.variables.prices(i)){solve_iteration-1} 
                case 'bisection'
                    brent_condition_1 = abs(solution.prices.values.(model.variables.prices(i)) - solution.prices.algorithm.prices.(model.variables.prices(i))(solve_iteration-1)) > method.update_prices.brent_threshold;
                    brent_condition_2 = abs(secant - solution.prices.values.(model.variables.prices(i))) < .5*abs(solution.prices.values.(model.variables.prices(i)) - solution.prices.algorithm.prices.(model.variables.prices(i))(solve_iteration-1));
                case 'secant'
                    brent_condition_1 = abs(solution.prices.algorithm.prices.(model.variables.prices(i))(solve_iteration-1) - solution.prices.algorithm.(model.variables.prices(i))(max(solve_iteration-2,1))) > method.update_prices.brent_threshold;
                    brent_condition_2 = abs(secant - solution.prices.values.(model.variables.prices(i))) < .5*abs(solution.prices.algorithm.prices.(model.variables.prices(i))(solve_iteration-1) - solution.prices.algorithm.prices.(model.variables.prices(i))(max(solve_iteration-2,1)));
            end 

            % update price
            if dekker_condition_1 && dekker_condition_2 && brent_condition_1 && brent_condition_2
                solution.prices.values.(model.variables.prices(i)) = secant;
                solution.prices.algorithm.update.(model.variables.prices(i)){solve_iteration} = 'secant';
                pause;
            else
                solution.prices.values.(model.variables.prices(i)) = bisection;
                solution.prices.algorithm.update.(model.variables.prices(i)){solve_iteration} = 'bisection';
            end
            fprintf('UPDATE (%s): %s ',model.variables.prices(i),  solution.prices.algorithm.update.(model.variables.prices(i)){solve_iteration}); 

        case 'bisection'
        case 'newton'
    end

    if i==1
        solution.prices.accuracy = abs(solution.prices.values.(model.variables.prices(i))- solution.prices.algorithm.prices.(model.variables.prices(i))(solve_iteration-1));
    elseif i==2
        solution.prices.accuracy = max(solution.prices.accuracy, abs(solution.prices.values.(model.variables.prices(i))- solution.prices.algorithm.prices.(model.variables.prices(i))(solve_iteration-1)));
    end
    
end


end


