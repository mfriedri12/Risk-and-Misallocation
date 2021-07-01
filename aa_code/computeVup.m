function [Vupeval] = computeVup(V0,param,prices,bkz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
w = prices.w;
r = prices.r;
th = param.th; lmd = param.lmd; dt = param.dt; bbar = param.bbar; bt = param.bt;
rho = param.rho; sgm = param.sgm; gm = param.gm;

nbkz = size(bkz,1);      

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
dif = 1; tol = 1e-2; gstol = 1e-3;
iter = 1;
[profit,~] = staticchoices(bkz(:,2),bkz(:,3),param,prices);
[ghx,ghw] = GaussHermite(10);

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
       if min(c) < -1e-5
           error('negative c')
       else
           c = max(c,0);
       end
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
           if min(c) < -1e-5
               error('negative c')
           else
               c = max(c,0);
           end
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
           if min(c) < -1e-5
               error('negative c')
           else
               c = max(c,0);
           end
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
               if min(c) < -1e-5
                   error('negative c')
               else
                   c = max(c,0);
               end
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
    Vupeval = max([u(:,1),bdu(:,1),bdu(:,2)],[],2);
end

