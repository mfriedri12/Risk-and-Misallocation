function [V0,KPOL,BPOL] = VFI(param,prices)
% Input: prices and parameters
% Output: optimal policy functions
w = prices.w;
r = prices.r;
th = param.th; lmd = param.lmd; dt = param.dt; bbar = param.bbar; bt = param.bt;
rho = param.rho; sgm = param.sgm; gm = param.gm;

ugrid = linspace(0,1,param.ngrid);
bgrid = param.bmin + (param.bmax - param.bmin).*ugrid.^2; bgrid = bgrid';
kgrid = param.kmax.*ugrid.^2; kgrid = kgrid';
%zgrid = linspace(param.zmin,param.zmax,param.ngrid); zgrid = zgrid';
zgrid = tauchen(0,param.rho, param.sgm, param.ngrid); % only slightly different than yours and let's me do the distribution
bkz = [repmat(bgrid,param.ngrid*param.ngrid,1),...
      repmat(repelem(kgrid,param.ngrid,1),param.ngrid,1),...
      repelem(zgrid,param.ngrid*param.ngrid,1)];
nbkz = size(bkz,1);      

[X1,X2,X3] = ndgrid(bgrid,kgrid,zgrid);

V0 = griddedInterpolant(X1,X2,X3,zeros(size(X1)),'makima','nearest');
dif = 1; tol = 1e-2;
iter = 1;
while dif >=tol
    
    [Vupeval,kpolup,bpolup] = computeVup(V0,param,prices,bkz);
    [Vdowneval,kpoldown,bpoldown] = computeVdown(V0,param,prices,bkz);
    upind = Vupeval >= Vdowneval;
    ufinal = Vdowneval; ufinal(upind) = Vupeval(upind);
    V1eval = reshape(ufinal,param.ngrid,param.ngrid,param.ngrid);
    V1 = griddedInterpolant(X1,X2,X3,V1eval,'makima','nearest');
    dif0 = (V1(bkz)-V0(bkz))./((V1(bkz)+V0(bkz))./2);
    dif = max(abs(dif0))
    iter = iter + 1;
    V0 = V1;
end

% policy functions
kpolfinal = kpoldown; kpolfinal(upind) =  kpolup(upind);
bpolfinal = bpoldown; bpolfinal(upind) = bpolup(upind);
kpoleval = reshape(kpolfinal,param.ngrid,param.ngrid,param.ngrid);
bpoleval = reshape(bpolfinal,param.ngrid,param.ngrid,param.ngrid);
KPOL = griddedInterpolant(X1,X2,X3,kpoleval,'makima','nearest');
BPOL = griddedInterpolant(X1,X2,X3,bpoleval,'makima','nearest');

end

