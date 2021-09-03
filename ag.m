function [x, stat] = ag(handles, x0, algparams)
% [x,stat] = ag(handles, x0, algparams)
% Nesterov's accelerated gradient to minimize f(x)
% See cplusag.m for description of the arguments
f = handles.f(x0);
L = handles.L;
ell = handles.ell;
y = x0; %y0
v = y; % v0
x = y;  % x0
mnv = f;
gamma = L;  %gamma0
maxit = algparams.maxit;
stat = struct('numit', 0);
for it = 1 : maxit
    stat.numit = it;    
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
        return
    end
    yk = y;
    gk = g;
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


end
