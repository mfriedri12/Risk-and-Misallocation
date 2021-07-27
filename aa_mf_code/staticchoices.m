function [profit,n] = staticchoices(k,z,param,prices)
% input: talent z(can be scalar or vector), capital stock k (can be scalar
% or vector,param, prices
% output: profit from running business, labor demand
% are these profit functions correct? look at overleaf doc
al = param.al; eta = param.eta; w = prices.w;

n = (eta.*(1-al).*exp(z).^(1-eta).*k.^(al*eta)./w).^(1/(1-(1-al)*eta));
profit =  exp(z).^(1-eta).*(k.^al.*n.^(1-al)).^(eta) - w.*n;
end

