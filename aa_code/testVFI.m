prices.w = 0.8;
prices.r = 0.02;
param.th = 0.9; param.lmd = 0.7; param.dt = 0.06; 
param.rho = 0.9; param.sgm = 1; param.gm = 0.5; param.bt = 0.96;
param.bmax = 100; param.bmin = -50; param.bbar = 10;
param.kmax = 100;
param.zmin = 1; param.zmax = 100;
param.ngrid = 10;
param.al = 0.33; param.eta = 0.15; 
test = VFI(param,prices);

