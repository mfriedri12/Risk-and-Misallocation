% GOLDEN Computes local maximum of univariate function on interval via Golden Search
% USAGE
%   [x,fval] = golden(f,a,b,P1,P2,...);
% INPUTS
%   f         : name of function of form fval=f(x)
%   optvar 
%   a,b       : left, right endpoints of interval
%   P1,P2,... : optional additional arguments for f
% OUTPUTS
%   x       : local maximum of f
%   fval    : function value estimate
%
% USER OPTIONS (SET WITH OPSET)
%   tol     : convergence tolerance

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
% vectorized by MRF 

function [x1,f1] = golden_mrf(f, tol, a,b,varargin)

%tol = optget('golden','tol',sqrt(eps));

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
d  = b-a;
x1 = a+alpha1*d;
x2 = a+alpha2*d;
f1 = feval(f,x1,varargin{:});
f2 = feval(f,x2,varargin{:});

d = alpha1*alpha2*d;
while max(d)>tol
  d = d*alpha2;
    
  selector = f2<f1 ; 
  
  newx2 = selector.*x1 + (1-selector).*(x2+d); 
  newx1 = selector.*(x1-d) + (1-selector).*x2; 
  
  x1 = newx1; 
  x2 = newx2;
  
  newf2 = selector.*f1                      + (1-selector).*feval(f,x2,varargin{:}); 
  newf1 = selector.*feval(f,x1,varargin{:}) + (1-selector).*f2;  
  
  f1 = newf1; 
  f2 = newf2; 
  
  %if f2<f1 % x2 is new upper bound
  %  x2 = x1; x1 = x1-d; 
  %  f2 = f1; f1 = feval(f,x1,varargin{:});
  %else     % x1 is new lower bound
  %  x1 = x2; x2 = x2+d; 
  %  f1 = f2; f2 = feval(f,x2,varargin{:});
  %end
  
end

% Return the larger of the two
selector = f2>f1;
x1 = selector.*x2 + (1-selector).*x1;
f1 = selector.*f2 + (1-selector).*f1;


%if f2>f1, x1 = x2; f1 = f2; end
