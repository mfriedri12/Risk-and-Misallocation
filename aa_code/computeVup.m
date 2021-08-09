function [Vupeval,kpol,bpol] = computeVup(V0,param,prices,bkz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
w = prices.w;
r = prices.r;
th = param.th; lmd = param.lmd; dt = param.dt; bbar = param.bbar; bt = param.bt;
rho = param.rho; sgm = param.sgm; gm = param.gm; infeaspen = param.infeaspen;

nbkz = size(bkz,1);      

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
dif = 1; tol = 1e-2; gstol = 1e-3;
iter = 1;
[profit,~] = staticchoices(bkz(:,2),bkz(:,3),param,prices);
[ghx,ghw] = GaussHermite(10);

infeas = bkz(:,1) < -th.*lmd.*(1-dt).*bkz(:,2)-bbar;

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
        % error
%         I think this top version is wrong because you pay back the (1+r)bkz(:,1) before being able to borrow again
%         bprime(:,1) = -th*lmd*(1-dt).*kprime(:,j)-bbar + min((1+r).*bkz(:,1),0);
        bprime(:,1) = -th*lmd*(1-dt).*kprime(:,j)-bbar;
        bprime(:,2) = profit + (1+r).*bkz(:,1) - (kprime(:,j)-(1-dt).*bkz(:,2));
        % doing this to avoid the situation where bprime(:,2) < bprime(:,1)
        % I think this happens because the grid has points that are not
        % feasible. i.e. no one would ever be lent money to end up in some
        % of these bkz combinations. maybe the interpretation of the
        % following is that if you are in one of those positions, the
        % amount you cannot pay back is forgiven,  but you still have to
        % pay the maximum you can pay back. not sure if this will mess with
        % optimal policy as agents engage in moral hazard
        bprime(:,2) = max(bprime,[],2);
        gb = bprime(:,2) - bprime(:,1);
        if min(gb) < 0
            error('bprime2 less than bprime1')
        end
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
       c(infeas) = 0;
       if min(c) < -1e-5
           error('negative c')
       else
           c = max(c,0);
       end
       bduu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
       bduu(infeas,jj) = -infeaspen;
       end
       bbduu = bprime;
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
           c(infeas) = 0;
           if min(c) < -1e-5
               error('negative c')
           else
               c = max(c,0);
           end
           uu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
           uu(infeas,jj) = -infeaspen;
           end

           bprime(:,1) = bprime(:,1).*(uu(:,1)>uu(:,2)) + y(:,1).*(uu(:,1)<=uu(:,2));
           bprime(:,2) = bprime(:,2).*(uu(:,1)<=uu(:,2)) + y(:,2).*(uu(:,1)>uu(:,2));
           gb = bprime(:,2) - bprime(:,1);
           bpoltemp= [y(:,1),bbduu]; % putting this here to back out policy function. see y changes on next line
           y(:,1) = bprime(:,1) + alpha1*gb; 
           y(:,2) = bprime(:,1) + alpha2*gb;
           gsbdif = max(gb);
        end
        [bdu(:,j),maxind] = max([uu(:,1),bduu(:,1),bduu(:,2)],[],2,'linear'); 
        bbdu(:,j) = bpoltemp(maxind);
        kbdu(:,j) = kprime(:,j);
    end
    
    % normal kprime golden search
    while gskdif >= gstol
        for j = 1:2
%             bprime(:,1) = -th*lmd*(1-dt).*x(:,j)-bbar + min((1+r).*bkz(:,1),0);
            bprime(:,1) = -th*lmd*(1-dt).*x(:,j)-bbar;
            bprime(:,2) = profit + (1+r).*bkz(:,1) - (x(:,j)-(1-dt).*bkz(:,2));
            bprime(:,2) = max(bprime,[],2);

            gb = bprime(:,2) - bprime(:,1);
            if min(gb) < 0
                error('bprime2 less than bprime1')
            end
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
           c(infeas) = 0;
           if min(c) < -1e-5
               error('negative c')
           else
               c = max(c,0);
           end
           bduu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
           bduu(infeas,jj) = -infeaspen;
           end
           bbduu = bprime;
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
               c(infeas) = 0;
               if min(c) < -1e-5
                   error('negative c')
               else
                   c = max(c,0);
               end
               uu(:,jj) = 1/(1-gm).*c.^(1-gm) + bt.*EV;
               uu(infeas,jj) = -infeaspen;
               end
               
               bprime(:,1) = bprime(:,1).*(uu(:,1)>uu(:,2)) + y(:,1).*(uu(:,1)<=uu(:,2));
               bprime(:,2) = bprime(:,2).*(uu(:,1)<=uu(:,2)) + y(:,2).*(uu(:,1)>uu(:,2));
               gb = bprime(:,2) - bprime(:,1);
               bpoltemp= [y(:,1),bbduu]; % putting this here to back out  policy function. see y changes on next line
               y(:,1) = bprime(:,1) + alpha1*gb; 
               y(:,2) = bprime(:,1) + alpha2*gb;
               gsbdif = max(gb);
            end
            [u(:,j),maxind] = max([uu(:,1),bduu(:,1),bduu(:,2)],[],2,'linear');
            bu(:,j) = bpoltemp(maxind);
            ku(:,j) = x(:,j);
        end         
        kprime(:,1) = kprime(:,1).*(u(:,1)>u(:,2)) + x(:,1).*(u(:,1)<=u(:,2));
        kprime(:,2) = kprime(:,2).*(u(:,1)<=u(:,2)) + x(:,2).*(u(:,1)>u(:,2));
        gk = kprime(:,2) - kprime(:,1);
        x(:,1) = kprime(:,1) + alpha1*gk; 
        x(:,2) = kprime(:,1) + alpha2*gk;
        gskdif = max(kprime(:,2) - kprime(:,1));
    end
    [Vupeval,maxind] = max([u(:,1),bdu(:,1),bdu(:,2)],[],2,'linear');
    % assign proper capital and bond choice
    bpoltemp= [bu(:,1),bbdu];
    kpoltemp = [ku(:,1),kbdu];
    bpol = bpoltemp(maxind);
    kpol = kpoltemp(maxind);  
    
    % make sure policy functions are correct. i.e. give the same value as Vupeval
    keval = repmat(kpol,1,size(ghx,1));
    beval = repmat(bpol,1,size(ghx,1));
    zrep = repmat(bkz(:,3),1,size(ghx,1));
    ghxeval = repmat(ghx',nbkz,1);
    zeval =  rho.*zrep + sqrt(2)*sgm.*ghxeval;
    V0eval = V0(beval,keval,zeval);
    EV = V0eval*ghw./(sqrt(pi));
    c = profit + (1+r).*bkz(:,1) - (kpol-(1-dt).*bkz(:,2)) - bpol;
    c(infeas) = 0;
    if min(c) < -1e-5
       error('negative c')
    else
       c = max(c,0);
    end
    testV = 1/(1-gm).*c.^(1-gm) + bt.*EV;
    testV(infeas) = -infeaspen;
    check = abs(Vupeval-testV);
    if max(check)> 1e-10
        error('policy functions not giving same value as Vupeval')
    end
end

