function [x, stats] = cplusag(handles, x0, algparams)
% [x,stats] = cplusag(handles, x0, algparams)
% Minimize a smooth, convex function f(x) using the C+AG method described
% in Karimi & Vavasis, 'Nonlinear conjugate gradient for smooth convex
% functions', 2021
% Input arguments:
%  handles: a struct with the following fields (plus others not used here):
%     fg: a function handle that takes x, an n-vector, and:
%         returns [f(x), grad f(x)]
%     L: a positive number, the smoothness modulus of f
%     ell: a nonnegative number, the strong convexity parameter of f
%  x0: an n-vector, the initial guess for the solution
%  algparams: a struct with the following fields
%     maxit: a positive integer, the iteration limit
%     gtol: a positive number, the termination criterion, i.e., stop
%        when norm(grad f(x_k)) <= gtol
%     use_restart_z: a logical variable; true means that xbar in the 
%        progress measure should use the more complicated formula
%        involving z described in the paper
% Output arguments:
%   x: an n-vector, the approximate optimizer
%   stats: a struct with the following field
%      numit: number of iterations
%      numag: number of iterations that are AG
% Note that this function does not count function/gradient evaluations.
% It is up to handles.fg to count these.

fg = handles.fg;
L = handles.L;
ell = handles.ell;
x = x0;
n = length(x0);
[f,g] = fg(x);
p = -g;
gamma = L;
v = x0;
phistar = f;
ag_it = 0;
norminitg = norm(g);
consec_cg = 0;
only_ag = false;
z = [];
zAz = nan(1);
xbar = x;
fbar = f;
gbar = g;

for it = 0 : algparams.maxit - 1
    % invariant is:
    % f(x_{k}) <= phi^*_{k}
    % where phi_k(x) = phi^*_k + gamma_k/2 * norm(x - v_k)^2   
    
    % Restart CG after 10*n+1 consecutive CG iterations
    if consec_cg >= 10*n + 1
        p = -g;
        consec_cg = 0;
    end
    
    
    [theta, newgamma] = ...
        compute_theta_gamma(L, gamma, ell);  %theta_{k}, gamma_{k+1}
    
    % compute x_{k+1}, y_k, v_{k+1}, phistar_{k+1} in this loop
    % maintain invariant that phistar_{k+1} >= f(x_{k+1}).
  
 
    % stepcount = 1 means: CG
    % stepcount = 2 means: restarted CG (i.e., steepest descent direction)
    % stepcount = 3 means: AG
    for stepcount = 1 : 3
        % Don't attempt CG or restarted CG if the 'only_ag' flag is set
        if stepcount < 3 && only_ag
            continue;
        end
        if stepcount == 2
            p = -g;
        end
        xsol = [];
        if stepcount <= 2
            % Take a CG step.  This section of code computes alphacg,
            % the CG line-search step length.
            xtilde = x + p/L;
            [ftilde, gtilde] = fg(xtilde);
            if norm(gtilde) < algparams.gtol
                xsol = xtilde;
                break
            end
            
            Ap = L * (gtilde - g);
            pAp = p' * Ap;
            alphacg = -g' * p / pAp;
            xcg = x + alphacg * p;  % x_{k+1}
            [fcg, gcg] = fg(xcg);
            if norm(gcg) < algparams.gtol
                xsol = xcg;
                break
            end
             % Update zupd if length(z) > 0
    
             prevxbar = xbar;
             prevfbar = fbar;
             prevgbar = gbar;
             
             % Compute xbar, f(xbar) (denoted fbar), and grad f(xbar)
             % (denoted gbar).  Two separate rules for xbar depending
             % on whether z has been selected.
             
             if length(z) > 0
                 zAp = z' * Ap;
                 delta = zAp / pAp;
                 z = z - delta * p;
                 zAz = zAz - 2 * delta * zAp + delta^2 * pAp;
                 alpha3 = -(gcg' * z) / zAz;
                 xbar = xcg + alpha3 * z;
                 [fbar, gbar] = fg(xbar);
                 if norm(gbar) < algparams.gtol
                     xsol = xbar;
                     break
                 end
             else
                 xbar = xcg;
                 fbar = fcg;
                 gbar = gcg;
             end           
            
             % Get the new function phi_{k+1}
             
            [newv, newphistar] = ...
                compute_v_phistar(theta, gamma, newgamma, ell, ...
                v, phistar, prevxbar, prevfbar, prevgbar); %v_{k+1}, phistar_{k+1} 
            
            % If the progress measure is satisfied by either fcg or
            % (in the case that z is set) fbar, then we accept the
            % step because adequate progress is made.
            
            if fbar <= newphistar || fcg <= newphistar    
                v = newv;
                phistar = newphistar;
                if stepcount == 1
                    consec_cg = consec_cg + 1;
                else
                    consec_cg = 0;
                end
                
                % This block of code computes the Hager-Zhang formula
                % for beta
                haty = gcg - g;
                assert(norm(haty) > 0)
                hatyp = haty' * p;
                betahz1 = ((haty - p * (2 * norm(haty)^2 / hatyp))' * gcg) / hatyp;
                betahz2 = -1 / (norm(p) * min(.01*norminitg, norm(gcg)));
                betahz = max(betahz1, betahz2);
                if isnan(betahz) || isinf(betahz)
                    betahz = 0;
                end
                if betahz == 0
                    consec_cg = 0;
                end
                newp = -gcg + betahz * p;
                break
            end                
        else   
            % We reach this point if stepcount == 3, i.e., an AG step
            
            ag_it = ag_it + 1;
            consec_cg = 0;
            
            % If this is the first AG step in a sequence, then we
            % set the only_ag flag and also compute the norm of g
            % for use in the test to terminate the block of AG steps.
            if ~only_ag
                only_ag = true;
                norm_ag_initg = norm(g);
            end           
            
            % Use y_k := nesterov formula
            yy = (theta * gamma * v + newgamma * x) / (gamma + theta * ell); %y_{k}
            [fyy, gyy] = fg(yy);
            
            if norm(gyy) < algparams.gtol
                xsol = yy;
                break
            end    
            
            xcg = yy - gyy/L;
            
            [v, phistar] = ...
                compute_v_phistar(theta, gamma, newgamma, ell, ...
                v, phistar, yy, fyy, gyy); %v_{k+1}, phistar_{k+1} 
            
            % Terminate AG when the gradient norm is reduced by a factor
            % of 1/4 compared to its value at the start of the block of AG
            % steps.
            if norm(gyy) < .25 * norm_ag_initg
                [fcg,gcg] = fg(xcg);
                xbar = xcg;
                fbar = fcg;
                gbar = gcg;
                if algparams.use_restart_z
                    z = v - xcg;
                    [fv,gv] = fg(v);
                    zAz = z' * (gv - gcg);
                end
                newp = -gcg;
                only_ag = false;
            end
        end
    end  
    if length(xsol) > 0
        break
    end    
    gamma = newgamma;
    x = xcg;
    f = fcg;
    g = gcg;
    p = newp;
end
if length(xsol) == 0
    x = xcg;
else
    x = xsol;
end
stats = struct('numit', it, 'numag', ag_it);
end
   

function [theta, newgamma] = ...
    compute_theta_gamma(L, gamma, ell) %theta_{k-1}, gamma_{k}
AA = L;
BB = gamma - ell;
CC = -gamma;
discr = BB^2 - 4 * AA * CC;
theta = (2 * CC) / (-BB - sqrt(discr));
newgamma = (1 - theta) * gamma + theta * ell;  % gamma_{k}
end


function [v, phistar] = ...
    compute_v_phistar(theta, gamma, newgamma, ell, v, phistar, yy, fyy, gyy) %v_k, phistar_k
oldv = v; % v_{k-1}
v = (1/newgamma) * ((1 - theta) * gamma * v + theta * ell * yy ...
    - theta * gyy);  % v_k
phistar = (1-theta) * phistar + theta * fyy - theta^2/(2*newgamma) * norm(gyy)^2 + ...
    (theta*(1-theta)*gamma/newgamma) * (ell/2 * norm(yy - oldv)^2 + ...
    gyy' * (oldv - yy)); % phistar_k
end
        
