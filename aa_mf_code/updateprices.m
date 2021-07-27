%% update prices 
%-------------------------------------------------------------------------%
% given distribution, update prices
%-------------------------------------------------------------------------%

function [prices,saveprices, interval] = updateprices(param,prices,saveprices,interval,distribution,iter)

%% prelims 
pricenames(1) = 'w';
pricenames(2) = 'r';
method = 'brent-dekker';


%% caluculate excess demand 
ugrid = linspace(0,1,param.ngrid);
bgrid = param.bmin + (param.bmax - param.bmin).*ugrid.^2; bgrid = bgrid';
kgrid = param.kmax.*ugrid.^2; kgrid = kgrid';
%zgrid = linspace(param.zmin,param.zmax,param.ngrid); zgrid = zgrid';
zgrid = tauchen(0,param.rho, param.sgm, param.ngrid); % only slightly different than yours and let's me do the distribution
[~,X2,X3] = ndgrid(bgrid,kgrid,zgrid);
bkz = [repmat(bgrid,param.ngrid*param.ngrid,1),...
      repmat(repelem(kgrid,param.ngrid,1),param.ngrid,1),...
      repelem(zgrid,param.ngrid*param.ngrid,1)];

[~,n] = staticchoices(X2,X3,param,prices);

excess.w = 1-dot(reshape(n,[],1), distribution); % negative means too much labor demand-- wage got to be higher. positive means too little labor demand -wage got to be lower. 
excess.r = dot(bkz(:,1), distribution); 

for i=1:length(pricenames)
    saveprices.prices.(pricenames(i))(iter) = prices.(pricenames(i));         
    saveprices.excess.(pricenames(i))(iter) = excess.(pricenames(i));
    fprintf('EXCESS (%s): %f ',pricenames(i), saveprices.excess.(pricenames(i))(iter)); 
end
fprintf('\n');


%% check interval  

if iter==1 
      
       % calculate excess value at all points in the interval (or grid, if it's 2 prices) abcd
       % note this code assumes that you're first run was at min-min (which is how it's coded in setup_solution if this option is true)
       for i=1:(length(pricenames)*2-1)  
           % for i=0 aka first run should be min min (min first price & min second price) excess should be - -
           if i==1 % max max (max first price & max second price) excess should be + + 
               prices.(pricenames(1)) = interval.(pricenames(1))(2);
               if length(pricenames)>1 
                    prices.(pricenames(2)) = interval.(pricenames(2))(2);
               end
           elseif i==2 % min max (min first price & max second price) excess should be - + 
               prices.(pricenames(1)) = interval.(pricenames(1))(1);
           elseif i==3 % max min (max first price & min second price) excess should be + -
               prices.(pricenames(1)) = interval.(pricenames(1))(2);
               prices.(pricenames(2)) = interval.(pricenames(2))(1);
           end
   
           fprintf('\nITERATION %d (CHECK GRID %d) ', iter, i); 
           for j=1:length(pricenames); fprintf('PRICE (%s): %f ', pricenames(j), prices.(pricenames(j))); end 
           fprintf('\n');
            
            disp('1. update solution.controls given solution.prices');       [V,KPOL,BPOL]  = VFI(param,prices); toc; 
            disp('2. update solution.distribution given solution.controls'); [distribution] = updatedistribution(param,KPOL,BPOL); toc; 
            
            excess.w = 1-dot(reshape(n,[],1), distribution); % negative means too much labor demand-- wage got to be higher. positive means too little labor demand -wage got to be lower. 
            excess.r = dot(bkz(:,1), distribution);
 
            for j=1:length(pricenames)
                saveprices.prices.(pricenames(j))(iter+i) = prices.(pricenames(j));         
                saveprices.excess.(pricenames(j))(iter+i) = excess.(pricenames(j));
                fprintf('EXCESS (%s): %f ',pricenames(j), saveprices.excess.(pricenames(j))(iter+i)); 
            end
            fprintf('\n');
       end
       
       % check that excesses are appropriate signs 
       if (length(pricenames)==1 && ...
          all((sign(saveprices.excess.(pricenames(1)))==[-1 1])))  ...
          ||...
          (length(pricenames)==2 && ...
          all((sign(saveprices.excess.(pricenames(1)))==[-1 1 -1 1])) && ...
          all((sign(saveprices.excess.(pricenames(2)))==[-1 1 1 -1])))
          fprintf('Interval OK! \n');
          iter = 2; %solve iterations 3 & 4 will be overwritten if it's 2 prices, but that's ok 
          if length(pricenames)==2
            prices.(pricenames(2)) = saveprices.prices.(pricenames(2))(2);
          end
       else 
          fprintf('\nWARNING: interval not valid! Paused.'); 
          pause 
       end
       
end 
    

%% update price DON'T RUN PAST HERE YET 

for i=1:length(pricenames)

    switch method 

        case 'brent-dekker'

            % using excess demand, update contrapoint and interval 
            if saveprices.excess.(pricenames(i))(iter) > 0
                contrapoint = min(interval.(pricenames(i)));  % demand too large so need to lower price 
                interval.(pricenames(i)) = [contrapoint prices.(pricenames(i))];
            else  
                contrapoint = max(interval.(pricenames(i)));  % demand too small so need to increase price 
                interval.(pricenames(i)) = [prices.(pricenames(i)) contrapoint];
            end   

            % calculate update options 
            secant = prices.(pricenames(i))-saveprices.excess.(pricenames(i))(iter)*...
                    (prices.(pricenames(i))-saveprices.prices.(pricenames(i))(iter-1))...
                   /(saveprices.excess.(pricenames(i))(iter)-saveprices.excess.(pricenames(i))(iter-1)); % secant_guess;
            bisection = (contrapoint + prices.(pricenames(i)))/2 ;  % bisection_guess

            % calculate update conditions  
            dekker_condition_1 = (saveprices.excess.(pricenames(i))(iter) > 0 && secant < bisection && secant > prices.(pricenames(i))) || (saveprices.excess.(pricenames(i))(iter) < 0 && secant > bisection && secant < prices.(pricenames(i))); % secant_guess betwen midpoint and guess 
            dekker_condition_2 = abs(saveprices.excess.(pricenames(i))(iter) - saveprices.excess.(pricenames(i))(iter-1))>0;
            switch saveprices.update.(pricenames(i)){iter-1} 
                case 'bisection'
                    brent_condition_1 = abs(prices.(pricenames(i)) - saveprices.prices.(pricenames(i))(iter-1)) > method.update_prices.brent_threshold;
                    brent_condition_2 = abs(secant - prices.(pricenames(i))) < .5*abs(prices.(pricenames(i)) - saveprices.prices.(pricenames(i))(iter-1));
                case 'secant'
                    brent_condition_1 = abs(saveprices.prices.(pricenames(i))(iter-1) - saveprices.(pricenames(i))(max(iter-2,1))) > method.update_prices.brent_threshold;
                    brent_condition_2 = abs(secant - prices.(pricenames(i))) < .5*abs(saveprices.prices.(pricenames(i))(iter-1) - saveprices.prices.(pricenames(i))(max(iter-2,1)));
            end 

            % update price
            if dekker_condition_1 && dekker_condition_2 && brent_condition_1 && brent_condition_2
                prices.(pricenames(i)) = secant;
                saveprices.update.(pricenames(i)){iter} = 'secant';
                pause;
            else
                prices.(pricenames(i)) = bisection;
                saveprices.update.(pricenames(i)){iter} = 'bisection';
            end
            fprintf('UPDATE (%s): %s ',pricenames(i),  saveprices.update.(pricenames(i)){iter}); 

        case 'bisection'
        case 'newton'
    end

    if i==1
        solution.prices.accuracy = abs(prices.(pricenames(i))- saveprices.prices.(pricenames(i))(iter-1));
    elseif i==2
        solution.prices.accuracy = max(solution.prices.accuracy, abs(prices.(pricenames(i))- saveprices.prices.(pricenames(i))(iter-1)));
    end
    
end


end


