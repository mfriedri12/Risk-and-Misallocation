%% update distribution 
%-------------------------------------------------------------------------%
% given controls, update distribuiton
% more useful view that I might rewrite as --
% reshape(solution.distribution.values, solution.states.dimensions.states)
%-------------------------------------------------------------------------%

function [dist] = updatedistribution(param,KPOL,BPOL)

method = 'default';

%% get grid
% rearranged to be zbk, that's how Eugene Tan's algorithm was written and
% I'm too lazy to rearrange, tho I have it on 'default' now which doesn't
% use his memory saving trick because our grids are so small 

ugrid = linspace(0,1,param.ngrid);
grid.b = param.bmin + (param.bmax - param.bmin).*ugrid.^2; grid.b  = grid.b';
grid.k = param.kmax.*ugrid.^2; grid.k = grid.k';
%zgrid = linspace(param.zmin,param.zmax,param.ngrid); zgrid = zgrid';
grid.z = tauchen(0,param.rho, param.sgm, param.ngrid);  % only slightly different than yours and let's me do the distribution
zbk = [ repmat(grid.z,param.ngrid*param.ngrid,1) ,...
    repmat(repelem(grid.b,param.ngrid,1),param.ngrid,1),...
      repelem(grid.k,param.ngrid*param.ngrid,1)
      ];  %rearranged zbk
nzbk = size(zbk,1);      

%[X1,X2,X3] = ndgrid(grid.b,grid.k,grid.z);


%% initialize distribution 
dist = (1/nzbk).*ones(nzbk,1);


%% get exogenous transition matricies 
[~,transition_exogenous] = tauchen(0,param.rho, param.sgm, param.ngrid); % second argument .. might be wrong 



%% get endogenous transition matricies

transition_endogenous = zeros(nzbk);
spacer = eye(prod(param.ngrid));
s = 1; 
policies.k = @(states) KPOL(states(2),states(3),states(1)); %switching order so that it's z
policies.b = @(states) BPOL(states(2),states(3),states(1));
policynames(1) = 'k';
policynames(2) = 'b'; 


for i=1:length(dist)
    
     for j=1:2
        
        % find transition probability from state point i to grid space j  
        transition_i_to_j = zeros(param.ngrid,1) ;
        lookup = policies.(policynames(j))(zbk(i,:));
        h = find(lookup<=grid.(policynames(j)),1); % first index greater 
        if h==1  % lookup value below grid 
            transition_i_to_j(1) = 1; 
        elseif isempty(h) % lookup value above grid 
            transition_i_to_j(end) = 1; 
        else 
            transition_i_to_j(h-1) = 1-((lookup-grid.(policynames(j))(h-1))/(grid.(policynames(j))(h) - grid.(policynames(j))(h-1)));
            transition_i_to_j(h)   = 1-((grid.(policynames(j))(h)-lookup)  /(grid.(policynames(j))(h) - grid.(policynames(j))(h-1)));
        end
       
        % add these probabilities to the previous ones  
        if j==1 
            transition_i_to_i = transition_i_to_j; 
        else
            transition_i_to_i = reshape(transition_i_to_i*transition_i_to_j',[],1);
        end
          
    end   
    
    % add exogenous state zeros or duplicate for full matrix 
    
    switch method
       case 'tan' 
          transition_endogenous(i,:) = reshape(transpose(transition_i_to_i*spacer(s,:)),[],1); 
        case 'default'
          transition_endogenous(i,:) = reshape(transpose(transition_i_to_i*transition_exogenous(s,:)),[],1); 
    end
    if s < param.ngrid 
      s = s + 1;
    else 
      s = 1;
    end
    
end



%% solve for distribution 
% start with some guess, and then iterate forward to solve for distribution 

iter = 0; 
max_iter = 1000; 
dif = 1e5;
tol = 1e-4;
while tol < dif && iter < max_iter 
    
    switch method
        case 'tan'
            intermediate = reshape(transpose(transition_endogenous)*dist, param.ngrid, nzbk/param.ngrid);
            guess = reshape(transpose(transition_exogenous )*intermediate, nzbk,1);
        case 'default'
            guess = transpose(transition_endogenous)*dist;
    end
    
dif = abs(max(guess - dist));
dist = guess;  
 
iter = iter + 1;
fprintf('iteration: %i, precision: %i \n', iter, dif);  
end


%% repermute dist back to bkz 

% switch back to bkz 
indicies = 0:(nzbk/param.ngrid -1); 
permutation = [];
for i = 1:(param.ngrid)
    permutation = [permutation ((param.ngrid).*indicies)+i];
end
dist = dist(permutation'); %test = zbk(permutation',[2 3 1]); %test

% check solution 
% for i=1:nzbk; testb(i) = policies.b(zbk(i,:)); end
% testk = 0;
% for i=1:nzbk; testk(i) = policies.k(zbk(i,:)); end
% [zbk testb' testk' dist dist>0] % looks good 
% sum(reshape(dist,3,nzbk/3),2) % looks good 


 
end