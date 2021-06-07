function [kbar] = kconst(aztgrid,param,prices)
% Compute credit constraint given function space
% param = [bt,sgm,T,gm,al,eta,kp,dt,phi,amin,amax,zmin,zmax,zta,rho,splinep,nT];
bt = param(1); sgm = param(2); T = param(3); gm = param(4);
al = param(5); eta = param(6); kp = param(7); dt = param(8); phi = param(9);
amin = param(10); amax = param(11);
zmin = param(12); zmax = param(13); zta = param(14); rho = param(15);
splinep = param(16);nT = param(17);
r = prices(1); w = prices(2);

% compute actual values of kbar at nodes

% value of k which minimizes lhs of incentive constraint
khat = (phi*(al*(1-eta))/(1 - phi + phi*dt + r))^((eta + al*(1-eta))/eta) ...
    *(((1-al)*(1-eta))/w)^((1-al)*(1-eta)/eta) .* aztgrid(:,2);
% evaluate incentive compatibility at khat
khatic = kbaric(khat,aztgrid,param,prices);
% firms for which kbar = 0
kbarzero = khatic(:,1) <= 0;

% newton method to find highest root
dif = 1;
tol = 1e-5;
kbar0 = 1e7.*ones(size(aztgrid,1),1);
while dif >= tol
    f = kbaric(kbar0,aztgrid,param,prices);
    kbar1 = kbar0 - (f(:,1)./f(:,2));
    kbar1 = max(kbar1,0);
    diff1 = kbar1-kbar0;
    diff1(kbarzero ==1) = 0;
    dif = max(abs(diff1));
    kbar0 = kbar1;
end
kbar = kbar0;
% question: Need to think about whether I still need this given I am
% now doing this bounded newton search thing
% update 1/3/2021: in the case that for all elements in vector, asset is
% equal to 0, diff gets evaluated as 0 and kbar gets assigned the first
% kbar1 which is not the correct value. For this reason, it is important to
% assign the proper value below.
kbar(kbarzero == 1) = 0;
end

