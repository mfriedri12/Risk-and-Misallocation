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
zgrid = logspace(log10(param.zmin),log10(param.zmax),param.ngrid); zgrid = zgrid';
bkz = [repmat(bgrid,param.ngrid*param.ngrid,1),...
      repmat(repelem(kgrid,param.ngrid,1),param.ngrid,1),...
      repelem(zgrid,param.ngrid*param.ngrid,1)];
nbkz = size(bkz,1);
  
      
    

[X1,X2,X3] = ndgrid(bgrid,kgrid,zgrid);
V0 = griddedInterpolant(X1,X2,X3,zeros(size(X1)),'spline','nearest');
Vu = griddedInterpolant(X1,X2,X3,zeros(size(X1)),'spline','nearest');
Vd = griddedInterpolant(X1,X2,X3,zeros(size(X1)),'spline','nearest');

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
    kprimea = bkz(:,2);
    % some agents cannot adjust upward so take this max for them
    kprimeb = max((profit + (1+r).*bkz(:,1) + (1-dt).*bkz(:,2) + bbar)./(1-th*lmd*(1-dt)),bkz(:,2));
    gk = kprimeb - kprimea;
    x(:,1) = kprimea + alpha1*gk; 
    x(:,2) = kprimea + alpha2*gk;
    gskdif = 1;
    while gskdif >= gstol
        for j = 1:2
            bprimea = -th*lmd*(1-dt).*x(:,j)-bbar;
            bprimeb = profit + (1+r).*bkz(:,1) - (x(:,j)-(1-dt).*bkz(:,2));
            gb = bprimeb - bprimea;
            y(:,1) = bprimea + alpha1*gb;
            y(:,2) = bprimea + alpha2*gb;
            gsbdif = 1;
            while gsbdif >= gstol
               for jj = 1:2
               % need to take expectation over z now
               beval = repmat(x(:,j),1,size(ghx,1));
               keval = repmat(y(:,jj),1,size(ghx,1));
               zrep = repmat(bkz(:,3),1,size(ghx,1));
               ghxeval = repmat(ghx',nbkz,1);
               zeval =  rho.*zrep + sqrt(2)*sgm.*ghxeval;
               Veval = V0(beval,keval,zeval);  
               EV = Veval*ghw./(sqrt(pi)); % Expected value next period
               c = profit + (1+r).*bkz(:,1) - (x(:,j)-(1-dt).*bkz(:,2)) - y(:,jj);
               uu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
               end
               
               bprimea = bprimea.*(uu(:,1)>uu(:,2)) + y(:,1).*(uu(:,1)<=uu(:,2));
               bprimeb = bprimeb.*(uu(:,1)<=uu(:,2)) + y(:,2).*(uu(:,1)>uu(:,2));
               gb = bprimeb - bprimea;
               y(:,1) = bprimea + alpha1*gb; 
               y(:,2) = bprimea + alpha2*gb;
               gsbdif = max(kprimeb - kprimea);
            end
            u(:,j) = uu(:,1);   
        end         
        kprimea = kprimea.*(u(:,1)>u(:,2)) + x(:,1).*(u(:,1)<=u(:,2));
        kprimeb = kprimeb.*(u(:,1)<=u(:,2)) + x(:,2).*(u(:,1)>u(:,2));
        gk = kprimeb - kprimea;
        x(:,1) = kprimea + alpha1*gk; 
        x(:,2) = kprimea + alpha2*gk;
        gskdif = max(kprimeb - kprimea);
    end
    
    Veval = reshape(u(:,1),param.ngrid,param.ngrid,param.ngrid);
    V1 = griddedInterpolant(X1,X2,X3,Veval,'spline','nearest');
    
    dif0 = (V1(bkz)-V0(bkz))./((V1(bkz)+V0(bkz))./2);
    dif = max(abs(dif0));
    iter = iter + 1;
    V0 = V1;
end

end

