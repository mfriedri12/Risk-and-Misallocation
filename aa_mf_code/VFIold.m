function [V0] = VFI(param,prices)
% Input: prices and parameters
% Output: optimal policy functions
w = prices.w;
r = prices.r;
th = param.th; lmd = param.lmd; dt = param.dt; bbar = param.bbar; bt = param.bt;
rho = param.rho; sgm = param.sgm; gm = param.gm;

ugrid = linspace(0,1,param.ngrid);
bgrid = param.bmin + (param.bmax - param.bmin).*ugrid.^2; bgrid = bgrid';
kgrid = param.kmax.*ugrid.^2; kgrid = kgrid';
zgrid = linspace(param.zmin,param.zmax,param.ngrid); zgrid = zgrid';
bkz = [repmat(bgrid,param.ngrid*param.ngrid,1),...
      repmat(repelem(kgrid,param.ngrid,1),param.ngrid,1),...
      repelem(zgrid,param.ngrid*param.ngrid,1)];
nbkz = size(bkz,1);      

[X1,X2,X3] = ndgrid(bgrid,kgrid,zgrid);

%% compute value of going up
V0 = griddedInterpolant(X1,X2,X3,zeros(size(X1)),'makima','nearest');

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
dif = 1; tol = 1e-2; gstol = 1e-3;
iter = 1;
[profit,~] = staticchoices(bkz(:,2),bkz(:,3),param,prices);
[ghx,ghw] = GaussHermite(10);
while dif >=tol
    % outer golden search over next period capital for upward adjusters
    kprime(:,1) = (1-dt).*bkz(:,2);
    % some agents cannot adjust upward so take this max for them
    kprime(:,2) = max((profit + (1+r).*bkz(:,1) + (1-dt).*bkz(:,2) + bbar)./(1-th*lmd*(1-dt)),kprime(:,1));
    gk = kprime(:,2) - kprime(:,1);
    x(:,1) = kprime(:,1) + alpha1*gk; 
    x(:,2) = kprime(:,1) + alpha2*gk;
    gskdif = 1;
    
    % check kprime boundaries
    for j = 1:2
        bprime(:,1) = -th*lmd*(1-dt).*kprime(:,j)-bbar + min((1+r).*bkz(:,1),0);
        bprime(:,2) = profit + (1+r).*bkz(:,1) - (kprime(:,j)-(1-dt).*bkz(:,2));
        gb = bprime(:,2) - bprime(:,1);
        y(:,1) = bprime(:,1) + alpha1*gb;
        y(:,2) = bprime(:,1) + alpha2*gb;

       %check bprime boundaries
       for jj = 1:2
       % need to take expectation over z now
       keval = repmat(kprime(:,j),1,size(ghx,1));
       beval = repmat(bprime(:,jj),1,size(ghx,1));
       zrep = repmat(bkz(:,3),1,size(ghx,1));
       ghxeval = repmat(ghx',nbkz,1);
       zeval =  rho.*zrep + sqrt(2)*sgm.*ghxeval;
       V0eval = V0(beval,keval,zeval);
       EV = V0eval*ghw./(sqrt(pi)); % Expected value next period
       c = profit + (1+r).*bkz(:,1) - (kprime(:,j)-(1-dt).*bkz(:,2)) - bprime(:,jj);
       bduu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
       end
        % normal bprime golden search
        gsbdif = 1;
        while gsbdif >= gstol
           for jj = 1:2
           % need to take expectation over z now
           keval = repmat(kprime(:,j),1,size(ghx,1));
           beval = repmat(y(:,jj),1,size(ghx,1));
           zrep = repmat(bkz(:,3),1,size(ghx,1));
           ghxeval = repmat(ghx',nbkz,1);
           zeval =  rho.*zrep + sqrt(2)*sgm.*ghxeval;
           V0eval = V0(beval,keval,zeval);
           EV = V0eval*ghw./(sqrt(pi)); % Expected value next period
           c = profit + (1+r).*bkz(:,1) - (kprime(:,j)-(1-dt).*bkz(:,2)) - y(:,jj);
           uu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
           end

           bprime(:,1) = bprime(:,1).*(uu(:,1)>uu(:,2)) + y(:,1).*(uu(:,1)<=uu(:,2));
           bprime(:,2) = bprime(:,2).*(uu(:,1)<=uu(:,2)) + y(:,2).*(uu(:,1)>uu(:,2));
           gb = bprime(:,2) - bprime(:,1);
           y(:,1) = bprime(:,1) + alpha1*gb; 
           y(:,2) = bprime(:,1) + alpha2*gb;
           gsbdif = max(gb);
        end
        testuu = [uu(:,1),bduu(:,1),bduu(:,2)];
        bdu(:,j) = max([uu(:,1),bduu(:,1),bduu(:,2)],[],2);   
    end
    
    % normal kprime golden search
    while gskdif >= gstol
        for j = 1:2
            I am not sure if lmd needs to be here because adjusting up
            bprime(:,1) = -th*lmd*(1-dt).*x(:,j)-bbar + min((1+r).*bkz(:,1),0);
            bprime(:,2) = profit + (1+r).*bkz(:,1) - (x(:,j)-(1-dt).*bkz(:,2));
            gb = bprime(:,2) - bprime(:,1);
            y(:,1) = bprime(:,1) + alpha1*gb;
            y(:,2) = bprime(:,1) + alpha2*gb;
            
           %check boundaries
           for jj = 1:2
           % need to take expectation over z now
           keval = repmat(x(:,j),1,size(ghx,1));
           beval = repmat(bprime(:,jj),1,size(ghx,1));
           zrep = repmat(bkz(:,3),1,size(ghx,1));
           ghxeval = repmat(ghx',nbkz,1);
           zeval =  rho.*zrep + sqrt(2)*sgm.*ghxeval;
           V0eval = V0(beval,keval,zeval);
           EV = V0eval*ghw./(sqrt(pi)); % Expected value next period
           c = profit + (1+r).*bkz(:,1) - (x(:,j)-(1-dt).*bkz(:,2)) - bprime(:,jj);
           bduu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
           end
            % normal golden search
            gsbdif = 1;
            while gsbdif >= gstol
               for jj = 1:2
               % need to take expectation over z now
               keval = repmat(x(:,j),1,size(ghx,1));
               beval = repmat(y(:,jj),1,size(ghx,1));
               zrep = repmat(bkz(:,3),1,size(ghx,1));
               ghxeval = repmat(ghx',nbkz,1);
               zeval =  rho.*zrep + sqrt(2)*sgm.*ghxeval;
               V0eval = V0(beval,keval,zeval);
               EV = V0eval*ghw./(sqrt(pi)); % Expected value next period
               c = profit + (1+r).*bkz(:,1) - (x(:,j)-(1-dt).*bkz(:,2)) - y(:,jj);
               uu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
               end
               
               bprime(:,1) = bprime(:,1).*(uu(:,1)>uu(:,2)) + y(:,1).*(uu(:,1)<=uu(:,2));
               bprime(:,2) = bprime(:,2).*(uu(:,1)<=uu(:,2)) + y(:,2).*(uu(:,1)>uu(:,2));
               gb = bprime(:,2) - bprime(:,1);
               y(:,1) = bprime(:,1) + alpha1*gb; 
               y(:,2) = bprime(:,1) + alpha2*gb;
               gsbdif = max(gb);
            end
            testuu = [uu(:,1),bduu(:,1),bduu(:,2)];
            u(:,j) = max([uu(:,1),bduu(:,1),bduu(:,2)],[],2);   
        end         
        kprime(:,1) = kprime(:,1).*(u(:,1)>u(:,2)) + x(:,1).*(u(:,1)<=u(:,2));
        kprime(:,2) = kprime(:,2).*(u(:,1)<=u(:,2)) + x(:,2).*(u(:,1)>u(:,2));
        gk = kprime(:,2) - kprime(:,1);
        x(:,1) = kprime(:,1) + alpha1*gk; 
        x(:,2) = kprime(:,1) + alpha2*gk;
        gskdif = max(kprime(:,2) - kprime(:,1));
    end
    testu = [u(:,1),bdu(:,1),bdu(:,2)];
    ufinal = max([u(:,1),bdu(:,1),bdu(:,2)],[],2);
    
    V1 needs to be envelope of Vup and Vdown
    V1eval = reshape(ufinal,param.ngrid,param.ngrid,param.ngrid);
    V1 = griddedInterpolant(X1,X2,X3,V1eval,'makima','nearest');
    
    dif0 = (V1(bkz)-V0(bkz))./((V1(bkz)+V0(bkz))./2);
    [dif,ind] = max(abs(dif0))
    iter = iter + 1;
    V0 = V1;
end

Vup = V0;


end

