%% setup method 
%-------------------------------------------------------------------------%
% define the grid, method & method paramters (algorithms), convergence thresholds, max iterations  
%-------------------------------------------------------------------------%

function [method] = setup_method


%% grids 
method.grid.endogenous.nodes = 10; %  might want to adjust for new models
method.grid.exogenous.nodes  = 8; %  might want to adjust for new models  

%method.grid.exogenous.spacing = 'gaussian quadrature';  % if 'custom' info on how to space is in "model"
method.grid.endogenous.spacing = 'even'; % griddedInterpolant doesn't like log 


%% algorithms
% this defines the methods used in each subscript 
method.update_controls.algorithm = 'VFI';  %'BLPT' 'EE'  'VFI' 'EGM' 'ECM' 'PFI'
method.update_controls.print = 5; 
method.update_controls.reset = false;
method.update_controls.enforce_grid_max = true;
method.update_controls.vtrick = false;  % value function transformation -- not really working right now

method.interpolate.toolkit.value_function = 'griddedInterpolant';  % 'griddedInterpolant', 'basis' 'none'; method for *interpolate* 
method.interpolate.basis.value_function = 'spline'; % either method for griddedInterpolant or basis for basis; sub-method for *interpolate* (used 'makima' before)
method.interpolate.extrapolate.value_function = 'linear'; % extrapolation method for griddedInterpolant for *interpolate* (used 'spline' before) 
method.interpolate.toolkit.control_function = 'griddedInterpolant';  % 'griddedInterpolant', 'basis' 'none'; method for *interpolate* 
method.interpolate.basis.control_function = 'spline'; % either method for griddedInterpolant or basis for basis; sub-method for *interpolate* (used 'makima' before)
method.interpolate.extrapolate.control_function = 'linear'; % extrapolation method for griddedInterpolant for *interpolate* (used 'spline' before) 
method.interpolate.toolkit.portfolio_problem = 'griddedInterpolant';  % 'griddedInterpolant', 'basis' 'none'; method for *interpolate* 
method.interpolate.basis.portfolio_problem = 'spline'; % either method for griddedInterpolant or basis for basis; sub-method for *interpolate* (used 'makima' before)
method.interpolate.extrapolate.portfolio_problem = 'linear'; % extrapolation method for griddedInterpolant for *interpolate* (used 'spline' before) 
method.interpolate.degree = 3; % degree for splines or for complete polynomials, only uses if *interpolate* is "basis" %method.interpolate.order = 3; % order for chebeshev polynomials? not sure what order means for splines? 
method.interpolate.breakpoints = {}; % MODEL SPECIFIC these are breakpoints used for splines if splines is  the basis

method.integrate.algorithm = 'gaussian quadrature';  %'tauchen' 'rowenhurst'; method for *integrate* note that gaussian quadrature requires a basis 
method.integrate.gaussian = 'qnwnorm'; % type of gaussian quadrature for *integrate*
method.integrate.nodes = 12; % nodes for *integrate* if method is gaussian quadrature 

method.maximize.algorithm = 'golden search';  % hill climbing algorithm for  'newton';
method.maximize.threshold = 1e-5;  % hill climbing algorithm for  'newton';
method.maximize.checkendpoint = true;

method.update_distribution.algorithm = 'tan'; % 'tan' 'eigenvector' 'default'
method.update_distribution.exogenous_transition = 'gaussian quadrature'; % 'tauchen' 'rowenhurst' 'gaussian quadrature'
method.update_distribution.nodes = 10; % nodes for tauchen 
method.update_distribution.reset = true; % reset distribution to ones each time (instead of using previous) 
method.update_distribution.print = 50; 

method.update_prices.interval = 1.25; % determines price search interval -- should be greater than 1 -- needs to be big enough or algorithm will break 
method.update_prices.check_interval = true; % false doesn't work right now, need defaults for previous prices 
method.update_prices.algorithm = 'brent-dekker'; % 'newton' 'brent'
method.update_prices.brent_threshold = 1e-10;


%% thresholds  
% these specify thresholds for convergence for the three update algorithms 
method.update_controls.threshold = 1e-5;
method.update_distribution.threshold = 1e-5;
method.update_prices.threshold = 1e-5;


%% iterations 
% these specify maximum iterations for the three update algorithms 
method.update_controls.iterations = 500; %if not converged by here probably not going to 
method.update_distribution.iterations = 1000; 
method.update_prices.iterations = 100; 


end