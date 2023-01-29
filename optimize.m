function [x, stats] = optimize(handles, x0, gtol, alg)
% [x,stat] = optimize(handles, gtol, alg)
% Minimize a convex smooth function specified by handles. 
% Input arguments:
%    handles: a struct with the following fields (plus others not used)
%          f: a function that takes x and returns f(x)
%          g: a function that takes x and returns grad f(x)
%          fg: a function that takes x and returns [f(x), grad f(x)]
%          L, ell: moduli of smoothness and strong convexity resp.  These
%             should be set to NaN if the optimizer is supposed to estimate
%             them.
%    x0: initial guess
%    gtol: termination when norm(grad f(x_k)) <= gtol.  Infinity norm used
%       for CG-Descent; 2-norm for others
%    alg: which algorithm to use.  This is a struct with
%           entries as follows.
%      'method': '+' means C+AG, 'a' means AG, 'c' means cgdescent
% Return arguments:
%    x: optimizer
%    stats: a struct with the following fields
%       'success': true or false (was tolerance achieved)
%       'algname': string with name of algorithm including options
%       'numit': number of iterations
%       'other': further information about the algorithm

n = length(x0);
algparams = alg.params;
algparams.gtol = gtol;
algparams_cell = namedargs2cell(algparams);
if strcmpi(alg.method, '+')
    [x, stats] = cplusag(handles, x0, algparams_cell{:});
elseif strcmpi(alg.method, 'a')
    [x, stats] = ag(handles, x0, algparams_cell{:});
elseif strcmpi(alg.method, 'c')  
    [x, status, stats1] = cg_descent(x0, 3*gtol/sqrt(n), handles.f, ...
        handles.g, handles.fg);
    stats.numit = stats1.iter;
    stats.success = status == 0;
    stats.algname = 'CGDescent';
    stats.other = '';
else
    error('unrecognized method')
end

        
    
    
    
    
    
    
    
        


end

