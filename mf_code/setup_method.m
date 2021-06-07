%% setup method 
%-------------------------------------------------------------------------%
% define the grid, method & method paramters (algorithms), convergence thresholds, max iterations  
%-------------------------------------------------------------------------%

function [method] = setup_method


%% grids 

method.grid.endogenous.nodes = 6; %  might want to adjust for new models
method.grid.exogenous.nodes  = 20; %  might want to adjust for new models  

method.grid.endogenous.spacing = 'even'; % griddedInterpolant doesn't like log 
method.grid.exogenous.spacing = 'custom';  % if 'custom' info on how to space is in "model"


%% algorithms
% this defines the methods used in each subscript 

method.update_controls.algorithm = 'BLPT';  %'BLPT' 'EE'  'VFI' 'EGM' 'ECM' 'PFI'
method.update_controls.print = 1; 
method.update_controls.reset = false;

method.interpolate.algorithm = 'griddedInterpolant';  %'griddedInterpolant', 'basis' 'none'; method for *interpolate* 
method.interpolate.basis = 'spline'; % either method for griddedInterpolant or basis for basis; sub-method for *interpolate* (used 'makima' before)
method.interpolate.extrapolation_method = 'makima'; % extrapolation method for griddedInterpolant for *interpolate* (used 'spline' before) 
method.interpolate.degree = 3; % degree for splines or for complete polynomials, only uses if *interpolate* is "basis"
method.interpolate.breakpoints = {}; % MODEL SPECIFIC these are breakpoints used for splines if spliens is  the basis

method.integrate.algorithm = 'gaussian quadrature';  %'tauchen' 'rowenhurst'; method for *integrate* note that gaussian quadrature requires a basis 
method.integrate.gaussian = 'qnwnorm'; % type of gaussian quadrature for *integrate*
method.integrate.nodes = 5; % notes for *integrate* if method is gaussian quadrature 

method.maximize.algorithm = 'golden search';  % hill climbing algorithm for  'newton';

method.update_distribution.algorithm =  'tan'; %  'eigenvector' 'default'
method.update_distribution.exogenous_transition = 'tauchen'; % 'rowenhurst' 'gaussian quadrature'
method.update_distribution.nodes = 5; % nodes for tauchen 
method.update_distribution.reset = false; % reset distribution to ones each time (instead of using previous) 
method.update_distribution.print = 40; 

method.update_prices.interval = 1.25; % determines price search interval -- should be greater than 1 -- needs to be big enough or algorithm will break 
method.update_prices.check_interval = false; 
method.update_prices.algorithm = 'brent-dekker'; % 'newton' 'brent'
method.update_prices.brent_threshold = 1e-10;


%% thresholds  
% these specify thresholds for convergence for the three update algorithms 

method.update_controls.threshold = 1e-5;
method.update_distribution.threshold = 1e-5;
method.update_prices.threshold = 1e-5;


%% iterations 
% these specify maximum iterations for the three update algorithms 

method.update_controls.iterations = 1000; 
method.update_distribution.iterations = 1000; 
method.update_prices.iterations = 100; 


end