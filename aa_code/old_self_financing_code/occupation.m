function [occ,y,k,l,yout] = occupation(azt,kbar,param,prices)
% Makes occupational decision given asset and talent
% and returns relavent variables
% param = [bt,sgm,T,gm,al,eta,kp,dt,phi,amin,amax,zmin,zmax,zta,rho,splinep,nT];
bt = param(1); sgm = param(2); T = param(3); gm = param(4);
al = param(5); eta = param(6); kp = param(7); dt = param(8); phi = param(9);
amin = param(10); amax = param(11);
zmin = param(12); zmax = param(13); zta = param(14); rho = param(15);
splinep = param(16);nT = param(17);
r = prices(1); w = prices(2);


k = zeros(length(azt),1);
l = zeros(length(azt),1);
% modern sector 
% unconstrained optimal input demand
lm_u = azt(:,2).*((al*(1-eta)/(r+dt))*((1-al)*(r+dt)/(al*w))^(1-al*(1-eta)))^(1/eta);
km_u = (al*w)/((1-al)*(r + dt)).*lm_u;
% constrained optimal input demand
lm_c = (((1-al)*(1-eta).*azt(:,2).^eta.*kbar.^(al*(1-eta)))./w).^(1/(eta + al*(1-eta)));
% check if constrained and assign optimal input demands
% WHY ARE SOME PEOPLE STILL CONSTRAINED EVEN THOUGH PHI = 1
c_ind = km_u > kbar; 
lm = (1-c_ind).*lm_u + c_ind.*lm_c;
km = (1-c_ind).*km_u + c_ind.*kbar;
% compute income in modern sector
yo(:,1) = azt(:,2).^eta.*(km.^al.*lm.^(1-al)).^(1-eta) - (r + dt).*km  - w.*lm - kp; 
% ymtest = azt(:,2).*Mm - kp;
% sum(yo(:,1)-ymtest)
% traditional sector
ltau = azt(:,2).*((1-eta)/w)^(1/eta);
yo(:,2) = azt(:,2).^eta.*ltau.^(1-eta) - w.*ltau;
% yttest = azt(:,2).*Mt;
% sum(yo(:,2)-yttest)

% sector choice
occ(:,1) = yo(:,1) >= yo(:,2) & yo(:,1) >= w;
occ(:,2) = yo(:,2) > yo(:,1) & yo(:,2) >= w;
occ(:,3) = w > max(yo(:,1),yo(:,2));

% check to see if all choose exactly one occupation
sumcheck = sum(occ,2);
sumcheck = sumcheck ~= 1;
sumcheck = sum(sumcheck);
if sumcheck >0 
    error('some agents not choosing exactly one occupation')
end

y = occ(:,1).*yo(:,1) + occ(:,2).*yo(:,2) + occ(:,3).*w;
k = k + occ(:,1).*km;
l = l + lm.*occ(:,1) + ltau.*occ(:,2);
ym = azt(:,2).^eta.*(km.^al.*lm.^(1-al)).^(1-eta) .* occ(:,1);
ytau = azt(:,2).^eta.*ltau.^(1-eta) .* occ(:,2);
yout = ym + ytau;
end

