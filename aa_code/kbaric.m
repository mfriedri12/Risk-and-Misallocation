function [c] = kbaric(k,azt,param,prices)
% evaluate incentive compatibility constraint
% c(1) returns value of incentive compatibility constraint
% c(2) returns derivative to help with finding root
% Compute credit constraint given function space
% param = [bt,sgm,T,gm,al,eta,kp,dt,phi,amin,amax,zmin,zmax,zta,rho,splinep,nT];
bt = param(1); sgm = param(2); T = param(3); gm = param(4);
al = param(5); eta = param(6); kp = param(7); dt = param(8); phi = param(9);
amin = param(10); amax = param(11);
zmin = param(12); zmax = param(13); zta = param(14); rho = param(15);
splinep = param(16);nT = param(17);
r = prices(1); w = prices(2);


c(:,1) =  (1+r).*azt(:,1) - phi*kp - (1-phi+phi*dt+r).*k ...
    + phi.*azt(:,2).^(eta/(eta + al*(1-eta))).*k.^((al*(1-eta))/(eta + al*(1-eta)))...
    .*((1-al)*(1-eta)/w).^((1-al)*(1-eta)/(eta + al*(1-eta))).*(1-((1-al)*(1-eta)));


c(:,2) =  - (1-phi+phi*dt+r) ...
    + phi.*azt(:,2).^(eta/(eta + al*(1-eta)))...
    .*((al*(1-eta))/(eta + al*(1-eta))).*k.^((al*(1-eta))/(eta + al*(1-eta)) - 1)...
    .*((1-al)*(1-eta)/w).^((1-al)*(1-eta)/(eta + al*(1-eta))).*(1-((1-al)*(1-eta)));
end 

