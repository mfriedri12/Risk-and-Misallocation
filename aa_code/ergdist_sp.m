 function [n] = ergdist_sp(Apol,agrid,zgrid,param,prices,ergda,ergdz,method)
%   Compute the ergodic distribution
% param = [bt,sgm,T,gm,al,eta,kp,dt,phi,amin,amax,zmin,zmax,zta,rho,splinep,nT];
bt = param(1); sgm = param(2); T = param(3); gm = param(4);
al = param(5); eta = param(6); kp = param(7); dt = param(8); phi = param(9);
amin = param(10); amax = param(11);
zmin = param(12); zmax = param(13); zta = param(14); rho = param(15);
splinep = param(16);nT = param(17);
r = prices(1); w = prices(2);

tgrid = 0:(T-1); tgrid = tgrid';

azt(:,1) = repmat(agrid,ergdz*T,1);
azt(:,2) = repmat(repelem(zgrid,ergda,1),T,1);
azt(:,3) = repelem(tgrid,ergda*ergdz,1);
azgrid(:,1) = repmat(agrid,ergdz,1);
azgrid(:,2) = repelem(zgrid,ergda,1);
naz = size(azgrid,1);

% compute transition probabilities
% is it a bad idea to use midpoints with a nonlinearly spaced grid
zmid = (zgrid(1:end-1)' + zgrid(2:end)')/2;
zmid = [1,zmid,inf];
zmargcdf = 1 - zmid.^(-zta);
zmargprob = zmargcdf(2:end) - zmargcdf(1:end-1);
Qzz = (1-rho).*repmat(zmargprob,ergdz,1);
Qzz = rho.*eye(ergdz) + Qzz;

% discrete optimal aprime
for iT = 1:T
    aprime0((iT-1)*naz+1:iT*naz,:) = Apol{iT}(azgrid);
end

% still need to do this because function might give weird answers for
% points not on initial value function grid due to curvature
aprime0 = min(max(aprime0,amin),amax);


% % nearest node method
% agridrep = repmat(agrid',length(aprime0),1);
% aprime0rep = repmat(aprime0,1,length(agrid));
% [minval, ind1] = min(abs(aprime0rep - agridrep),[],2);
% aprime = agrid(ind1);
% aprimeind = zeros(length(azt),ergda);
% % convert to linear indices
% ind0 = 1:length(azt); ind0 = ind0';
% linind = sub2ind(size(aprimeind),ind0,ind1);
% aprimeind(linind) = ones(length(azt),1);


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
aprimeind = zeros(length(azt),ergda);
% convert to linear indices
ind0 = 1:length(azt); ind0 = ind0';
linind = sub2ind(size(aprimeind),ind0,lbind);
aprimeind(linind) = 1-lottprob;
linind = sub2ind(size(aprimeind),ind0,ubind);
aprimeind(linind) = lottprob;
aprimeind = sparse(aprimeind);
% generate A matrix
erase = sparse(kron(eye(ergdz),ones(ergda,ergda)));
A = sparse(eye(length(erase)));
for t = 0:T-1
    tempmat = repmat(aprimeind(azt(:,3) == t,:)',ergdz,1).*erase;
    A = tempmat*A;
end

switch method
case 1
    % initialize uniformly
    f0 = ones(ergda*ergdz,1);
    f0 = f0./length(f0);
    f1 = zeros(ergda*ergdz,1);

    ndif = 1;
    ctr = 0; 
    while ndif > 1e-15
        ctr = ctr+1;
        ftilde = A*f0;
        f1 = kronm({eye(ergda),Qzz'},ftilde);
        ndif = norm(abs(f1-f0));
        f0 = f1;
    end
    f = f1;

case 2
    Q2 = kronm({eye(ergda),Qzz'},A);
    Qsparse = sparse(Q2);
    [f,d] = eigs(Qsparse,1,1);
    f = f./sum(f);
    f = max(f,0);
    f = f./(sum(f));
    % ergdif0 = sum(abs(f2 -  Qsparse*f2),1);
end


n = zeros(T*ergda*ergdz,1);
n(1:ergda*ergdz) = f;
A = sparse(eye(length(erase)));
for t=0:T-2
    tempmat = repmat(aprimeind(azt(:,3) == t,:)',ergdz,1).*erase;
    A = tempmat*A;
    n(ergda*ergdz*(t+1)+1:ergda*ergdz*(t+2)) = A*f;
end

end

