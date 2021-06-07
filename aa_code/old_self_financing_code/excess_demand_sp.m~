function [out] = excess_demand_sp(param,prices,outputind,saveind)
%Computes excess demands given parameters and prices

% param = [bt,sgm,T,gm,al,eta,kp,dt,phi,amin,amax,zmin,zmax,zta,rho,splinep,nT];
bt = param(1); sgm = param(2); T = param(3); gm = param(4);
al = param(5); eta = param(6); kp = param(7); dt = param(8); phi = param(9);
amin = param(10); amax = param(11);
zmin = param(12); zmax = param(13); zta = param(14); rho = param(15);
splinep = param(16);nT = param(17);
r = prices(1); w = prices(2);

% initialize talent grid
zgrid = logspace(log10(zmin),log10(zmax),splinep);

% initialize wealth grid
ugrid = linspace(0,1,splinep); 
ugridT = linspace(0,1,nT); 
agrid = amin + (amax - amin).*ugrid.^2;
agrid = agrid'; zgrid = zgrid';

% initialize age grid
tgrid = 0 + ((T-1) - 0).*ugridT; tgrid = tgrid';

% initialize joint talent, wealth, age grid
aztgrid(:,1) = repmat(agrid,splinep*nT,1);
aztgrid(:,2) = repmat(repelem(zgrid,splinep,1),nT,1);
aztgrid(:,3) = repelem(tgrid,splinep^2,1);
azgrid(:,1) = repmat(agrid,splinep,1);
azgrid(:,2) = repelem(zgrid,splinep,1);
naz = size(azgrid,1);


%% COMPUTE CREDIT CONSTRAINTS

kbar = kconst(aztgrid,param,prices);

%% OCCUPATIONAL DECISION
% occ gives the occupational choices of each point on grid
% y gives the occupational choices of each point on grid
[occ,y,~,~,~] = occupation(aztgrid,kbar,param,prices);

%% VALUE FUNCTION ITERATION

% compute nodes and weights for gauss-laguerre quadrature
[glx,glw] = GaussLaguerre(30,0);

oldind = aztgrid(:,3) == tgrid(end);

dif = 1;
diflast = 10;
aprimelast = ones(size(oldind));
clast = aprimelast;
ulast = aprimelast;
aprimetol = 1e-3;

tol = 1e-2;

[X1,X2] = ndgrid(agrid,zgrid);

V{1} = griddedInterpolant(X1,X2,zeros(size(X1)),'makima','nearest');
V{T} = griddedInterpolant(X1,X2,zeros(size(X1)),'makima','nearest');

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
iter = 1;
politer = 1;
ctr = 0;
while dif >= tol
    if iter > 20
    % case when value function not converging
    if dif < 5e-2
        % let this go
        break;
    else
        % consider this an error
        l_excess = NaN; c_excess = NaN; ak_excess = NaN;    
    end
    return;
    end
    for iT = T:-1:1

    % use golden search to find optimal a'
    aprimea = ones(naz,1).*amin;
    aprimeb = max((1+r).*azgrid(:,1) + y((iT-1)*naz+1:iT*naz,:),0);

    d = aprimeb - aprimea;
    x(:,1) = aprimea + alpha1*d; 
    x(:,2) = aprimea + alpha2*d;
    dif2 = 1;
    while dif2 >= aprimetol
        for j = 1:2
            if iT == T
               % compute continuation value    
               % gaussian quadrature for old cohort
               aeval = repelem(x(:,j),1,size(glx,1));
               zpreval = repmat(glx',naz,1);
               Veval = V{1}(aeval,zpreval+1);
               pdfexp = zta.*(zpreval+1).^(-zta-1).*exp(zpreval);
               test = pdfexp*glw;
               feval = Veval.*pdfexp;
               CV = bt.*gm.*(rho.*V{1}(x(:,j),azgrid(:,2)) + (1-rho).*feval*glw);
               c = (1+r).*azgrid(:,1) + y((iT-1)*naz+1:iT*naz,:) - x(:,j);
               u(:,j) = (c.^(1-sgm))./(1-sgm) + CV;
           else
               % just update states for other cohorts
               CV = bt.*V{iT+1}(x(:,j),azgrid(:,2));
               c = (1+r).*azgrid(:,1) + y((iT-1)*naz+1:iT*naz,:) - x(:,j);
               u(:,j) = (c.^(1-sgm))./(1-sgm) + CV;
           end
        end
                
        aprimea = aprimea.*(u(:,1)>u(:,2)) + x(:,1).*(u(:,1)<=u(:,2));
        aprimeb = aprimeb.*(u(:,1)<=u(:,2)) + x(:,2).*(u(:,1)>u(:,2));
        d = aprimeb - aprimea;
        x(:,1) = aprimea + alpha1*d; 
        x(:,2) = aprimea + alpha2*d;
        dif2 = max(aprimeb - aprimea);
    end
    
    

    VVeval = reshape(u(:,1),splinep,splinep);
    aprimea = reshape(aprimea,splinep,splinep);

    c = reshape(c,splinep,splinep);
    VV{iT} = griddedInterpolant(X1,X2,VVeval,'makima','nearest');
    Apol{iT} = griddedInterpolant(X1,X2,aprimea,'makima','nearest');
    Apold{iT} = griddedInterpolant(X1,X2,aprimea,'nearest','nearest');
    Cpol{iT} = griddedInterpolant(X1,X2,c,'makima','nearest');
    
    if iT == T
        dif0 = (VV{iT}(azgrid(:,1),azgrid(:,2)) - V{iT}(azgrid(:,1),azgrid(:,2)))...
            ./ ((VV{iT}(azgrid(:,1),azgrid(:,2)) + V{iT}(azgrid(:,1),azgrid(:,2)))./2);
        dif = max(abs(dif0));
        iter = iter + 1;
    end

    V{iT} = VV{iT};    
    end
end 


for iT = 1:T
    apoleval((iT-1)*naz+1:iT*naz,:) = Apol{iT}(azgrid);
    tveval((iT-1)*naz+1:iT*naz,:) = V{iT}(azgrid);
end

% Check Euler Equation
elhs = 1 + r;
for iT = 1:T
    
if iT < T    
    ctoday{iT} = (Cpol{iT}(azgrid)).^(-sgm);
    cprime{iT} = (Cpol{iT+1}(Apol{iT}(azgrid),azgrid(:,2))).^(-sgm);
    erhs{iT} = ctoday{iT}./cprime{iT}./bt;

else
    ctoday{iT} = (Cpol{iT}(azgrid)).^(-sgm);
    % gaussian quadrature for old cohort
    aeval0 = Apol{iT}(azgrid);
    aeval = repelem(aeval0,1,size(glx,1));
    zpreval = repmat(glx',naz,1);    
    ceval = (Cpol{1}(aeval,zpreval+1)).^(-sgm);
    feval = ceval.*zta.*(zpreval+1).^(-zta-1).*exp(zpreval);
    cprime{iT} = rho.*(Cpol{1}(aeval0,azgrid(:,2))).^(-sgm) + (1-rho).*feval*glw;
    erhs{iT} = ctoday{iT}./cprime{iT}./bt./gm;
end
end

%% Ergodic Distribution

ergda = 50; ergdz = 150;
ergdp = ergda;
ugrid2 = linspace(0,1,ergda);
ugrid3 = linspace(0.001,0.999,ergdz);
agrid2 = amin + (amax - amin).*ugrid2.^2; agrid2 = agrid2';
% initialize talent grid so that gridpoints are at evenly spaced percentiles
zgrid2 = (1-ugrid3).^(-1/zta); zgrid2 = zgrid2'; 
tgrid2 = 0:T-1;tgrid2 = tgrid2';
[n] = ergdist_sp(Apol,agrid2,zgrid2,param,prices,ergda,ergdz,2);

%% Check Market Clearing
aztgrid2(:,1) = repmat(agrid2,ergdz*T,1);
aztgrid2(:,2) = repmat(repelem(zgrid2,ergda,1),T,1);
aztgrid2(:,3) = repelem(tgrid2,ergda*ergdz,1);
azgrid2(:,1) = repmat(agrid2,ergdz,1);
azgrid2(:,2) = repelem(zgrid2,ergda,1);
naz2 = size(azgrid2,1);
oldind2 = aztgrid2(:,3) == T-1;
youngind2 = aztgrid2(:,3) == 0;

kbar2 = kconst(aztgrid2,param,prices);
[occ2,y2,k2,l2,yout2] = occupation(aztgrid2,kbar2,param,prices);
for iT = 1:T
    aprime2((iT-1)*naz2+1:iT*naz2,:) = Apol{iT}(azgrid2);
end

% to prevent impossible asset choices. Because policy function not
% explicitly bounded, still might give asset choices outside of bounds for
% points that were not in inital value function grid. 
aprime2 = min(max(aprime2,amin),amax);
c2 = (1+r).*aztgrid2(:,1) + y2 - aprime2;

% % Walras Check
% in = (1+r).*aztgrid2(:,1) + yout2 + w.*occ2(:,3); 
% out = c2 + aprime2 + (r+dt).*k2 + w.*l2 + kp.*occ2(:,1); 
% t = [in,out,in-out];
% max(abs(t(:,3)))
% budgetcheck = (in-out)'*n
% % take out labor payments
% in = (1+r).*aztgrid2(:,1) + yout2; 
% out = c2 + aprime2 + (r+dt).*k2 + kp.*occ2(:,1);
% % ensure what we take out sums to zero
% lminus = w*(occ2(:,3) - l2)'*n
% lcheck = (in-out)'*n
% % check if assets saved = assets available in next period
% asum = [aztgrid2(:,1)'*n,aprime2'*n,(aztgrid2(:,1)-aprime2)'*n];
% % take out capital saved since asum shows they are about equal
% in = (r).*aztgrid2(:,1) + yout2; 
% out = c2 + (r+dt).*k2 + kp.*occ2(:,1);
% asumminus = asum(:,3)
% asumcheck = (in-out)'*n
% % take out good payments
% in = (r).*aztgrid2(:,1); 
% out = (r+dt).*k2 + kp.*occ2(:,1);
% cminus = (c2 - yout2)'*n
% ccheck = (in-out)'*n
% %take out depreciation and fixed cost
% in = (r).*aztgrid2(:,1); 
% out = (r).*k2;
% minus = (dt.*k2 + kp.*occ2(:,1))'*n
% kcheck = (in-out)'*n

c_demand = c2'*n; c_supply = yout2'*n;
ak_demand = k2'*n; ak_supply = aprime2'*n;
l_demand = l2'*n; l_supply = occ2(:,3)'*n;

if ((c_demand + c_supply)/2) == 0
    c_excess = (c_demand - c_supply);
else
    c_excess = (c_demand - c_supply)./((c_demand + c_supply)./2);
end

if ((l_demand + l_supply)/2) == 0
    l_excess = (l_demand - l_supply);
else
    l_excess = (l_demand - l_supply)./((l_demand + l_supply)./2);
end

if ((ak_demand + ak_supply)/2) == 0
    ak_excess = (ak_demand - ak_supply);
else
    ak_excess = (ak_demand - ak_supply)./((ak_demand + ak_supply)./2);
end

%% Compute Moments

% Elasticity of Earnings
indvec = zeros(ergda,1); indvec(1) = 1;
indvec = logical(repmat(indvec,ergdz,1));
Yp = log(y2(indvec)); Yp = repelem(Yp,ergdz,1); Yp = [ones(size(Yp)),Yp];
Yc = log(y2(indvec)); Yc = repmat(Yc,ergdz,1);

zmid = (zgrid2(1:end-1)' + zgrid2(2:end)')/2;
zmid = [1,zmid,inf];
zmargcdf = 1 - zmid.^(-zta);
zmargprob = zmargcdf(2:end) - zmargcdf(1:end-1);
Qzz = (1-rho).*repmat(zmargprob,ergdz,1);
Qzz = rho.*eye(ergdz) + Qzz;
pi = reshape(Qzz',ergdz^2,1);
nn = repelem(zmargprob',ergdz,1);
W = diag(nn.*pi);
nu = (Yp'*W*Yp)^(-1) * Yp'*W*Yc;
nu = nu(2);

% Slope of talent-rank
% Rp = repelem(1-zgrid2.^(-zta),ergdz,1); Rp = [ones(size(Rp)),Rp];
% Rc = repmat(1-zgrid2.^(-zta),ergdz,1);
% zmid = (zgrid2(1:end-1)' + zgrid2(2:end)')/2;
% zmid = [1,zmid,inf];
% zmargcdf = 1 - zmid.^(-zta);
% zmargprob = zmargcdf(2:end) - zmargcdf(1:end-1);
% Qzz = (1-rho).*repmat(zmargprob,ergdz,1);
% Qzz = rho.*eye(ergdz) + Qzz;
% pi = reshape(Qzz',ergdz^2,1);
% nn = repelem(zmargprob',ergdz,1);
% W = diag(nn.*pi);
% b = (Rp'*W*Rp)^(-1) * Rp'*W*Rc;
% nu2 = b(2);

% Consumption Ratio
temp = 1:T;
csum = (bt.*(1+r)).^((temp-1)./(sgm));
csum = sum(csum);
psum = (bt.*(1+r)).^((temp-T)./(sgm));
psum = sum(psum);
% parent consumption in final time period
Cp = Cpol{T}(azgrid2);
% parent lifetime average consumption
Cp = Cp.*psum./T;
% gaussian quadrature for old cohort
aeval0 = Apol{T}(azgrid2);
aeval = repelem(aeval0,1,size(glx,1));
zpreval = repmat(glx',naz2,1);    
ceval = Cpol{1}(aeval,zpreval+1);
feval = ceval.*zta.*(zpreval+1).^(-zta-1).*exp(zpreval);
% expectation of chilid consumption in first period
Cc = rho.*Cpol{1}(aeval0,azgrid2(:,2)) + (1-rho).*feval*glw;
% child lifetime lifetime average consumption
Cc = Cc.*csum./T;
xi = (Cc./Cp)'*n(oldind2)./sum(n(oldind2));

% Modern Sector Share
lmd = occ2(:,1)'*n/T;

% Average Modern Sector firm size
if (occ2(:,1)'*n) > 0
    lbar = (occ2(:,1).*n)'*l2/(occ2(:,1)'*n);
else
    lbar = 0;
end

% wealth rank correlation
Teval = T-1;
Tind = aztgrid2(:,3) == Teval;
% calculate ranks for each wealth state
sumind = eye(ergda); sumind = repmat(sumind,ergdz*nT,1);
apdf = n'*sumind; apdf = apdf./nT;
acdf = cumsum(apdf);
R = griddedInterpolant(agrid2,acdf,'nearest','nearest');

% compute transition probabilities
zmid = (zgrid2(1:end-1)' + zgrid2(2:end)')/2;
zmid = [1,zmid,inf];
zmargcdf = 1 - zmid.^(-zta);
zmargprob = zmargcdf(2:end) - zmargcdf(1:end-1);
Qzz = (1-rho).*repmat(zmargprob,ergdz,1);
Qzz = rho.*eye(ergdz) + Qzz;
% discrete optimal aprime
for iT = 1:T
    aprime0((iT-1)*naz2+1:iT*naz2,:) = Apol{iT}(azgrid2);
end
% still need to do this because function might give weird answers for
% points not on initial value function grid due to curvature
aprime0 = min(max(aprime0,amin),amax);
% lottery method
lottprob = nan(size(aprime0));
agridrep = repmat(agrid',length(aprime0),1);
aprime0rep = repmat(aprime0,1,length(agrid));
lb = aprime0rep >= agridrep;
% find the interval in which aprime is
lbind = sum(lb,2);
% have to do this in case some agents choose amax
ubind = min(lbind+1,ergda);
lb = agrid(lbind); 
ub = agrid(ubind);
%agents who choose amax
bdd = lb == ub;
lottprob(~bdd) = (aprime0(~bdd) - lb(~bdd))./(ub(~bdd)-lb(~bdd));
lottprob(bdd) = 1;
aprimeind = zeros(length(aztgrid2),ergda);
% convert to linear indices
ind0 = 1:length(aztgrid2); ind0 = ind0';
linind = sub2ind(size(aprimeind),ind0,lbind);
aprimeind(linind) = 1-lottprob;
linind = sub2ind(size(aprimeind),ind0,ubind);
aprimeind(linind) = lottprob;
aprimeind = sparse(aprimeind);
% generate Ac matrix to figure out transition probs from age 1 to age 30 wealth for children
erase = sparse(kron(eye(ergdz),ones(ergda,ergda)));
Ac = sparse(eye(length(erase)));
for t = 0:Teval-2
    tempmat = repmat(aprimeind(aztgrid2(:,3) == t,:)',ergdz,1).*erase;
    Ac = tempmat*Ac;
end
% generate Ap matrix to figure out transition prob from age 30  to age T wealth for parents
Ap = sparse(eye(length(erase)));
for t = Teval-1:T-1
    tempmat = repmat(aprimeind(aztgrid2(:,3) == t,:)',ergdz,1).*erase;
    Ap = tempmat*Ap;
end
clearvars -except Q Ac ergda ergdz Qzz Ap W n Tind aztgrid2 R T...
    l_excess ak_excess c_excess nu xi lmd lbar Apol Apold...
    outputind saveind bt rho gm kp zta phi agrid2 ...
    c2 k2 l2 occ occ2 yout2 aprime2 prices kbar2 param;
%out of memory error here need to fix need to sum over talent levels
Q = (Ac*kronm({eye(ergda),Qzz'},Ap))';
% sum over tommorrows asset states
summat = repmat(eye(ergda),ergdz,1);
Q = Q*summat;
% percentage of agents who start at (a,z) [row] today and end at (a)[column] tomorrow
W = repmat(n(Tind),1,ergda).*Q;
% percentage of agents who start at (a) [row] today and end at (a)[column] tomorrow
summat = repmat(eye(ergda),1,ergdz);
W = summat*W;
clear Q Ac Ap;
W = reshape(W',ergda^2,1);
W = sparse(1:length(W),1:length(W),W);
% rank of parents
Rp = R(agrid2); Rp = [ones(size(Rp)),Rp];
Rp = repelem(Rp,ergda,1);
% rank of children
Rc = R(agrid2);
Rc = repmat(Rc,ergda,1);

xi = (Rp'*W*Rp)^(-1) * Rp'*W*Rc;
xi = xi(2);




%% Output
if outputind == 1
    out = l_excess; 
elseif outputind == 2
    out = ak_excess;
elseif outputind == 3
    out = [l_excess,ak_excess];
elseif outputind == 4
    out = [l_excess,ak_excess,c_excess];
elseif outputind == 5
    out = [l_excess,ak_excess,c_excess,nu,xi,lmd,lbar];
end

if saveind == 1
%     filename = strcat('/scratch/aa6553/LCM/output/','bt',num2str(bt),'rho',num2str(rho),'gm',num2str(gm),'kp',num2str(kp),'zta',num2str(zta),'T',num2str(T),'phi',num2str(phi),'.mat');
    filename = strcat('C:\Users\anuae\OneDrive\Documents\MATLAB\AnuModel2019\LCM\output\','bt',num2str(bt),'rho',num2str(rho),'gm',num2str(gm),'kp',num2str(kp),'zta',num2str(zta),'phi',num2str(phi),'.mat');
    save(filename);
end

end

