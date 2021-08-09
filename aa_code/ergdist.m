 function [n] = ergdist(KPOL,BPOL,bgrid,kgrid,zgrid,param,prices,method)
%   Compute the ergodic distribution
w = prices.w;
r = prices.r;
th = param.th; lmd = param.lmd; dt = param.dt; bbar = param.bbar; bt = param.bt;
rho = param.rho; sgm = param.sgm; gm = param.gm;

nb = length(bgrid); nk = length(kgrid); nz = length(zgrid);

bkz = [repmat(bgrid,nk*nz,1),...
      repmat(repelem(kgrid,nb,1),nz,1),...
      repelem(zgrid,nb*nk,1)];
nbkz = size(bkz,1);  

infeas = bkz(:,1) < -th.*lmd.*(1-dt).*bkz(:,2)-bbar;

% compute transition probabilities (z are rows, z' are columns)
zmid = (zgrid(1:end-1)' + zgrid(2:end)')/2;
zmid = [-inf,zmid,inf];
ztoday = repmat(zgrid,1,length(zmid));
ztomorrow = repmat(zmid,length(zgrid),1);
zmidcdf = normcdf(ztomorrow-rho.*ztoday,0,sgm);
Qzz = zmidcdf(:,2:end) - zmidcdf(:,1:end-1);

% optimal bprime and kprime
bprime0 = BPOL(bkz);
kprime0 = KPOL(bkz);

bprime0 = min(max(bprime0,param.bmin),param.bmax);
kprime0 = min(max(kprime0,0),param.kmax);

% lottery method
blottprob = nan(size(bprime0));
bgridrep = repmat(bgrid',length(bprime0),1);
bprime0rep = repmat(bprime0,1,length(bgrid));
blb = bprime0rep >= bgridrep;
% find the interval in which bprime is
blbind = sum(blb,2);
% have to do this in case some agents choose max
bubind = min(blbind+1,nb);
blb = bgrid(blbind); 
bub = bgrid(bubind);

klottprob = nan(size(kprime0));
kgridrep = repmat(kgrid',length(kprime0),1);
kprime0rep = repmat(kprime0,1,length(kgrid));
klb = kprime0rep >= kgridrep;
% find the interval in which bprime is
klbind = sum(klb,2);
% have to do this in case some agents choose max
kubind = min(klbind+1,nk);
klb = kgrid(klbind); 
kub = kgrid(kubind);

% blb and klb
llind = nb*(klbind-1)+blbind;
% blb and kub
luind = nb*(kubind-1)+blbind;
% bub and klb
ulind = nb*(klbind-1)+bubind;
% bub and kub
uuind = nb*(kubind-1)+bubind;

%agents who choose bmax
bbdd = blb == bub;
blottprob(~bbdd) = (bprime0(~bbdd) - blb(~bbdd))./(bub(~bbdd)-blb(~bbdd));
blottprob(bbdd) = 1;

%agents who choose kmax
kbdd = klb == kub;
klottprob(~kbdd) = (kprime0(~kbdd) - klb(~kbdd))./(kub(~kbdd)-klb(~kbdd));
klottprob(kbdd) = 1;

% bkprime = sparse(nbkz,nb*nk);
bkprime = zeros(nbkz,nb*nk);

ind0 = (1:length(bkz))';

% assign ll
linind = sub2ind(size(bkprime),ind0,llind);
bkprime(linind) = (1-blottprob).*(1-klottprob);
% assign lh
linind = sub2ind(size(bkprime),ind0,luind);
bkprime(linind) = (1-blottprob).*(klottprob);
% assign hl
linind = sub2ind(size(bkprime),ind0,ulind);
bkprime(linind) = blottprob.*(1-klottprob);
% assign hh
linind = sub2ind(size(bkprime),ind0,uuind);
bkprime(linind) = blottprob.*klottprob;
% have to transpose so multiply by f0 correctly
bkprime = bkprime';

% maybe need to assign first and then make sparse for efficiency
% aprimeind = zeros(length(bkz),ergda);
% % convert to linear indices
% linind = sub2ind(size(aprimeind),ind0,lbind);
% aprimeind(linind) = 1-blottprob;  
% linind = sub2ind(size(aprimeind),ind0,ubind);
% aprimeind(linind) = blottprob;
% aprimeind = sparse(aprimeind);

% generate A matrix
% erase = sparse(kron(eye(nz),ones(nb*nk,nb*nk)));
erase = sparse(kron(eye(nz),ones(nb*nk,nb*nk)));

A = repmat(bkprime,nz,1).*erase;

switch method
case 1
    % initialize uniformly
    % only start agents in feasible gridpoints
    f0 = zeros(nbkz,1);
    f0(~infeas) = 1;
    f0 = f0./sum(f0);
%     
%     f0 = zeros(nbkz,1);
%     bminind = bkz(:,1) == param.bmin;
%     f0(bminind) = ones(nk*nz,1);
%     f0 = f0./sum(f0);
%     f1 = zeros(nbkz,1);
    

    ndif = 1;
    ctr = 0; 
    while ndif > 1e-15
        ctr = ctr+1;
        ftilde = A*f0;
        f1 = kronm({eye(nb*nk),Qzz'},ftilde);
        ndif = norm(abs(f1-f0));
        f0 = f1;
    end
    n = f1;

case 2
    Q2 = kronm({eye(ergda),Qzz'},A);
    Qsparse = sparse(Q2);
    [f,d] = eigs(Qsparse,1,1);
    f = f./sum(f);
    f = max(f,0);
    f = f./(sum(f));
    % ergdif0 = sum(abs(f2 -  Qsparse*f2),1);
end

end

