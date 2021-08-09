prices.w = 0.8;
prices.r = 0.02;
param.th = 0.9; param.lmd = 0.7; param.dt = 0.06; 
param.rho = 0.9; param.sgm = 1; param.gm = 0.6; param.bt = 0.93;
param.bmax = 100; param.bmin = -20; param.bbar = 10;
param.kmax = 100;
param.zmin = -5*param.sgm; param.zmax = 5*param.sgm;
param.ngrid = 6;
param.al = 0.33; param.eta = 0.15; 
param.infeaspen = 0;
L = 1; %labor supply
[V,KPOL,BPOL] = VFI(param,prices);

nb = 5; nk = 5; nz = 5;
% one issue: everyone is choosing bmin
% another issue: is everyone just endowed with 1 unit of labor they supply
% inelastically? Is this reflected in income?
bgrid = param.bmin + (param.bmax - param.bmin).*linspace(0,1,nb).^2; bgrid = bgrid';
kgrid = param.kmax.*linspace(0,1,nk).^2; kgrid = kgrid';
zgrid = linspace(param.zmin,param.zmax,nz); zgrid = zgrid';
n = ergdist(KPOL,BPOL,bgrid,kgrid,zgrid,param,prices,1);
bkz = [repmat(bgrid,nk*nz,1),...
      repmat(repelem(kgrid,nb,1),nz,1),...
      repelem(zgrid,nb*nk,1)];
[~,ld] = staticchoices(bkz(:,2),bkz(:,3),param,prices);
bd = BPOL(bkz);
kd = KPOL(bkz);

excl = ld'*n - L
excb = bd'*n






    