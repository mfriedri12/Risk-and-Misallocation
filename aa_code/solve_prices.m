clear;clc; 

% In BKS they say sgm = 1.5 but this doesnt make sense to me
% remember B(1+r)< 1 in order for asset space to be compact

bt = 0.989453125000000; sgm = 1.5; T = 50; gm = 0.4534;
al = 0.33; eta = 0.15; kp = 0.6132; dt = 0.06;
% when magnitude of rho is large, value function doesnt converge. 
% in this case, when look at value of test cdf, it is clear that the 
% gaussian quadrature is not getting close to approximating function well
% as for rho = 5, cdf is 0.5 for t=0 and for rho = 2, cdf is 0.9 "---".
% this might be the reason it isnt converging
phi = 0.1; zta = 2.0900; rho = 0.4;
% I think this will cause a problem because don't know how to 
% implement tranversality condition numerically, so for now 
% limiting to case of no borrowing
% % amin = -w/r;
amin = 0; amax = 1000;
zmin = (1-0.001)^(-1/zta); zmax = (1-0.999)^(-1/zta);
splinep = 50; splinek = 3;  nT = T;
param = [bt,sgm,T,gm,al,eta,kp,dt,phi,amin,amax,zmin,zmax,zta,rho,splinep,nT]
options = optimoptions(@fsolve,'Display','iter');

solvemethod = 1
switch solvemethod
case 1 %--------------------------- BISECTION ----------------------------%
    akztol = 0.01;
    lztol = 0.01;
    akz = 10;
    rl = -dt; rh = 1/bt-1; rh = 0.2;
    while abs(akz)>akztol && abs(rh-rl)>1e-5
    r = (rl + rh)/2;
    lz = 10;
    wl = 0; wh = 2;
    while abs(lz)>lztol && abs(wh-wl)>1e-5
        w = (wl + wh)/2;
        term1 = al*w/((1-al)*(r+dt));
        term2 = (al*(1-eta)/(r+dt)*((1-al)*(r+dt)/(al*w))^(1-al*(1-eta)))^(1/eta);
        term3 = ((1-eta)/w)^(1/eta);
        Mm = (term1^al * term2)^(1-eta) - (r+dt)*term1*term2 - w*term2;
        Mt = term3^(1-eta) - w*term3;
        Zt = w/Mt;
        Zm = kp/(Mm-Mt);

        excess = excess_demand_sp(param,[r,w],4,0);
        lz = excess(1); akz = excess(2); cz = excess(3);
        display(strcat('r = ',num2str(r),', w = ',num2str(w)))
        display(strcat('Zt = ',num2str(Zt),' Zm = ',num2str(Zm)))
        display(strcat('lz = ',num2str(lz),' akz = ',num2str(akz)))
        disp('--------------------------------')
        if lz > 0
            wl = w;
        else
            wh = w;
        end
    end
    if akz >0
        rl = r;
    else
        rh = r;
    end
    end
    excess = excess_demand_sp(param,[r,w],4,1);
    [p,fval] = fsolve(@(X)excess_demand_sp(param,X,3,0),[r,w],options)
    excess = excess_demand_sp(param,p,4,1);

    
case 2 %----------------------------- FSOLVE -----------------------------%
    r =0.012232, w =0.94659
    [p,fval] = fsolve(@(X)excess_demand_sp2(param,X),[r,w],options)
    [cz,lz,akz] = excess_demand_spout(param,p);
    
case 3 %----------------------------- DEBUG ------------------------------%
    r = 0.048343, w = 0.83057
    [cz,lz,akz] = excess_demand_sp(param,[r,w]);
end

%% COMMENTS

% 7/11/2020
% I had an issue that even though labor market and goods market were
% clearing it seemed as if the asset market was not clearing.  
% What happened was that that even though the absolute
% difference between a_demand, a_supply is very small, becasuse they are
% both such small magnitudes (almost no capital in economy), the percentage
% difference looks large.
% Note, by walras law, the fact that the other two markets were clearing
% should imply that the last market clears, but since we are using a
% numerical approximation, the other markets clearing meant that the excess
% demand was "very small", but the sum of the "very small" excess demands
% is what is carried over to the asset market and this sum still creates a
% large excess demand in percentage terms becasue the magnitude of assets
% traded is so small. 
% (7/28/2020) Another reason that the asset excess demand can be large even
% though the other markets clear is because the walras residual is divided by r.
% so even though the residual <1, it is multiplied 20x by dividing by r.
% Thus need residualto be a lot smaller

% I also made the the interpolation as makima rather than spline (spline 
% seemed to be adding too much curvature) and also 
% changed zta so that the support of z would be large enough to evaluate 
% all of the glx points. I think the reason the value function was not
% I think need to have a good approximation of value function for 
% convergence, which is why I did these things. Also switched from
% frank copula to simple persistence initally because for rho values far
% from zero, the cdf was very properly evaluated so assumed that value
% function as well was probably being evaluated poorly. after changing to
% this simple persistence, i further was convinced to use it because it
% seems more simple and more in line with BKS, MX, Moll.

% 8/1/2020
% When Zm is negative it implies that modern sector will never be preferred
% to traditional sector
% I am still having trouble generating enough entry into the modern sector.
% I think this is because if the specification is z^eta*(k^(al)*l^(1-al))^1-eta for
% the modern sector and z^eta(l^1-eta) for the traditional sector, the
% modern sector really doesnt have an advantage over traditional sector. In
% MX, they add a exogenous productivity boost to the modern sector. I think
% they do this to generate entry into the modern sector. it doesnt sit well with me. 

% 7/27/2020, Plan going forward:
% 1. Need to show TFP is decreasing in phi, but in order to do this, need to
% have a large chunk of people working in modern sector. it seems like i
% cant generate this at the moment. need to double check whether something
% going wrong in my credit constraint function. going to debug phi = 1 case
% 2. After that, need to show that gamma and rho can either amplify or
% attenuate these losses
% 3. Finally need to calibrate and test numerically.

