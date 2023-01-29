function [x, stats] = ag(handles, x0, algparams)
arguments
    handles          (1,1) struct
    x0               (:,1) double
    algparams.maxit  (1,1) double = Inf
    algparams.maxfg  (1,1) double = Inf
    algparams.gtol   (1,1) double = 1e-6
end 
% [x,stat] = ag(handles, x0, algparams)
% Nesterov's accelerated gradient to minimize f(x)
% Input arguments:
%  handles: a struct with the following fields (plus others not used here):
%     fg: a function handle that takes x, an n-vector, and:
%         returns [f(x), grad f(x)]
%     L: a positive number, the smoothness modulus of f
%        If L is NaN, this means true L is unknown, and ag should
%        estimate it.
%     ell: a nonnegative number, the strong convexity parameter of f
%        If ell is NaN, this means it is unknown.
%     fgcount(): a function with no arguments that returns the number of
%        function-gradient evaluations so far.  This is invoked only if
%        algparams.maxfg < Inf.
%  x0: an n-vector, the initial guess for the solution
%
%  The following are optional keyword arguments
%
%  maxit: a positive integer, the iteration limit.  Default: Inf
%    Inf means no limit.
%  maxfg: a positive integer, the max number of fg evaluations
%    default: Inf.  Note: if this option is used, the handle
%    must provide handle.fgcount() to count fg evaluations.
%    Inf means no limit.
%  gtol: a positive number, the termination criterion, i.e., stop
%        when norm(grad f(x_k)) <= gtol.  Default: 1e-6.
% Output arguments:
%   x: an n-vector, the approximate optimizer
%   stats: a struct with the following fields
%      numit: number of iterations
%      algname: a string with the name of the algorithm used
%      success: was convergence achieved (boolean)
%      other: unused
% Note that this function does not count function/gradient evaluations.
% It is the responsibility of handles.fg to count these.
%

[f,g] = handles.fg(x0);
L = handles.L;
ell = handles.ell;
if ~isinf(algparams.maxfg)
    fgcount = @()(handles.fgcount());
else
    fgcount = @()0;
end
fginit = fgcount();
y = x0; %y0
v = y; % v0
x = y;  % x0
mnv = f;
maxit = algparams.maxit;
Lnan = isnan(L) || isnan(ell) || isinf(L) || isinf(ell); 
stats = struct('numit', 0, ...
    'other', '', ...
    'success', false, ...
    'algname', sprintf('AG(HaveL=%s)', string(~Lnan)));
if Lnan
    L = estimateL2(x,f,g,handles.fg,1);
    if isnan(L) || isinf(L)
        error('Unable to estimate L from initial x')
    end
    if ~isnan(ell) && ell ~= 0
        error('Input error: cannot specify ell when L=NaN')
    end
    ell = 0;
end 
gamma = L;
it = 0;
while it < algparams.maxit && fgcount() - fginit < algparams.maxfg
    it = it + 1;
    stats.numit = it;    
    AA = L;
    BB = gamma - ell;
    CC = -gamma;
    discr = BB^2 - 4 * AA * CC;
    assert(discr >= 0)
    theta = (2 * CC) / (-BB - sqrt(discr));
    oldgamma = gamma;  %gamma_k
    gamma = (1 - theta) * oldgamma + theta * ell; %gamma_{k+1}    
    y = (theta * oldgamma * v + gamma * x) / (oldgamma + theta * ell); %y_k
    [f,g] = handles.fg(y);
    if norm(g) <= algparams.gtol
        stats.success = true;
        return
    end
    yk = y;
    gk = g;
    if Lnan
        L = estimateL1(y, f, g, handles.fg, L);
    end
    x = y - g / L;  % x_{k+1}
    k1 = ell * norm(y - v) ^ 2 / 2 + g' * (v - y);
    mnv = (1 - theta) * mnv + theta * f - ...
        theta^2 * norm(g)^2 / (2 * gamma) + ...
        theta * (1-theta) * oldgamma * k1 / gamma;    
    v = ((1 - theta)*oldgamma*v + ell*theta*yk - theta*gk) / gamma;
    if any(isnan(v)) || any(isinf(v))
        error('NaN/Inf encountered in AG')
    end
end
stats.numit = maxit + 1;
end
