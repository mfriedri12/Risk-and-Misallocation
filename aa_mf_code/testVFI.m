clear all; 
clc; 

%% parameters
param.th = 0.9; param.lmd = 0.7; param.dt = 0.06; 
param.rho = 0.9; param.sgm = 1; param.gm = 0.6; param.bt = 0.93;
param.bmax = 100; param.bmin = -20; param.bbar = 10;
param.kmax = 100;
param.zmin = -5*param.sgm; param.zmax = 5*param.sgm;
param.ngrid = 10;
param.al = 0.33; param.eta = 0.15; 

%% prices 
interval.w = [.05, 1.05];
interval.r = [-.02, .02]; 
prices.w = interval.w(1); %0.8;
prices.r = interval.r(1); %0.02;
saveprices.excess.r  = 1e5;  saveprices.excess.w = 1e5; 
saveprices.prices = prices; 
saveprices.update.r{1}  = 'bisection'; saveprices.update.w{1} = 'bisection';

%% solution
tic; 
iter = 1; 
maxiter = 10;
dif = 1e5; 
tol = 1e-5;
while tol < dif && iter < maxiter
    fprintf('\nITERATION %d :', iter); 

    disp('1. update solution.controls given solution.prices');       [V,KPOL,BPOL]  = VFI(param,prices); toc; 
    disp('2. update solution.distribution given solution.controls'); [distribution] = updatedistribution(param,KPOL,BPOL); toc; 
    disp('3. update solution.prices given solution.distribution');   [prices,saveprices, interval] = updateprices(param,prices,saveprices,interval,distribution,iter); toc; 

    iter = iter + 1; 
end

fprintf('\n SOLVED!\n'); 

toc; 
