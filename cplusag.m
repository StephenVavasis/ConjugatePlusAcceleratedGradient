function [x, stats] = cplusag(handles, x0, algparams)
arguments
    handles                  (1,1)   struct
    x0                       (:,1)   double
    % The following are optional keyword arguments.  For Matlab releases
    % prior to 2022b, usage is: 
    %      cplusag(handles, x0, 'maxit', 500, 'switchfac', 0.7) 
    % For Matlab releases starting with 2022b, usage is:
    %      cplusag(handles, x0, maxit = 500, switchfac = 0.7)     
    algparams.maxit          (1,1)   double   = Inf
    algparams.maxfg          (1,1)   double   = Inf
    algparams.gtol           (1,1)   double   = 1e-6
    algparams.use_elab_progr (1,1)   logical  = false
    algparams.switch_fac     (1,1)   double   = 0.8
    algparams.check_ag_term  (1,1)   double   = 8
    algparams.bp_c1          (1,1)   double   = Inf
    algparams.bp_c2          (1,1)   double   = 0
    algparams.bp_c3          (1,1)   double   = Inf
end
% [x,stats] = cplusag(handles, x0, algparams)
% Minimize a smooth, convex function f(x) using the C+AG method described
% in Karimi & Vavasis, 'Nonlinear conjugate gradient for smooth convex
% functions', 2023
% Version 2.0
% Input arguments:
%  handles: a struct with the following fields (plus others not used here):
%     fg: a function handle that takes x, an n-vector, and:
%         returns [f(x), grad f(x)]
%     L: a positive number, the smoothness modulus of f
%        If L is NaN, this means true L is unknown, and cplusag should
%        estimate it.
%     ell: a nonnegative number, the strong convexity parameter of f
%        If ell is NaN, this means it is unknown.
%     fgcount(): a 0-argument function that counts the number of
%     function-gradient evaluations.  Used only if algparams.maxfg < Inf
%  x0: an n-vector, the initial guess for the solution
%  maxit (keyword): a positive integer, the iteration limit.  Default: Inf
%     (no limit).
%  maxfg (keyword): a positive integer, a limit on the number of
%      function-gradient evaluations; Inf for no limit; the functionn is 
%      is expected to provide this via handles.fgcount().
%  gtol (keyword): a positive number, the termination criterion, i.e., stop
%      when norm(grad f(x_k)) <= gtol.  Default: 1e-6.
%  use_elab_progr (keyword): a logical variable; 'true' means that xbar in the 
%      progress measure should use the more elaborate formula
%      involving w_k,z_k described in the paper. Default: false
%  switch_fac (keyword): a number in (0,1) that is the coefficient in the 
%      sufficient-decrease test for AG steps before switching back to
%      CG.  Default: 0.8.
%  check_ag_term (keyword): a positive integer that specifies the interval between
%      AG iterations when to check the sufficient-decrease termination
%      test.  Default: 8
%  For the following constants, refer to Dai & Yuan, Convergence Properties
%  of Beale-Powell Restart Algorithm, 1998.
%  bp_c1: A number in [0,1). Beale-Powell constant c1 in BP restart 
%     test #1 which is: restart if 
%         abs(prevg' * g) >= c1 * norm(grad)^2
%     Default: Inf (i.e., test is inactive)
%  bp_c2: A number in [0,1).  Beale-Powell constant c2 in BP restart
%     test #2 which is: restart if
%         p'*g >= -c2 * norm(g)^2
%     Default: 0 (i.e., restart only if p fails to be a descent direction)
%  bp_c3: A number in (1,infty).  Beale-Powell constant C3 in BP restart
%         p'*g <= -c3 * norm(g)^2
%     Default: Inf (i.e., condition is never active).
% Output arguments:
%  x: an n-vector, the approximate optimizer
%  stats: a struct with the following fields
%      numit: number of iterations
%      algname: a string with the name of the algorithm used
%      success: was convergence achieved (boolean)
%      other: string reporting the number of iterations that are AG

fg = handles.fg;
L = handles.L;
ell = handles.ell;
Lnan = isnan(L) || isnan(ell) || isinf(L) || isinf(ell); 
if ~isinf(algparams.maxfg)
    fgcount = @()(handles.fgcount());
else
    fgcount = @()0;
end
fginit = fgcount();
x = x0;
n = length(x0);
[f,g] = fg(x);
norminitg = norm(g);
p = -g; % First CG search direction is negative gradient
Ap = [];
v = x0;
phistar = f;
ag_it = 0;
consec_cg = 0;
consec_ag = 0;
xsol = [];
only_ag = false;
if Lnan
    L = estimateL2(x,f,g,fg,1);
    if isnan(L) || isinf(L)
        error('Unable to estimate L from initial x')
    end
    if ~isnan(ell) && ell ~= 0
        error('Input error: cannot specify nonzero ell when L=NaN')
    end
    ell = 0;
end 
gamma = L;

% wz contains the additional scalars and vectors of the
% elaborate progress measure that are needed for recurrent updates
wz = struct('baralpha', 0, ...
    'barbeta', 0, ...
    'w', [], ...
    'Aw', [], ...
    'z', [], ...
    'Az', []);

it = 0;
while it < algparams.maxit && fgcount() - fginit < algparams.maxfg
    it = it + 1;

    % compute x_{k+1}, y_k, v_{k+1}, phistar_{k+1} in this loop
    % maintain invariant that phistar_{k+1} >= f(x_{k+1}).
    % where phistar_k(x) = phi^*_k + gamma_k/2 * norm(x - v_k)^2   
    [theta, gammanew] = compute_theta_gamma(L, gamma, ell);  
    xnew = [];
    fnew = NaN;
    gnew = [];
    xbarnew = [];
    fbarnew = NaN;
    gbarnew = [];
    wznew = struct();
    pnew = [];
    
 
    % whichsteptype = 1 means: CG
    % whichsteptype = 2 means: restarted CG (i.e., steepest descent direction)
    % whichsteptype = 3 means: AG
    
    for whichsteptype = 1 : 3
        % Within this loop: 'break' with xsol = [] means successful step;
        % go to next iteration.  'continue' means unsuccessful step; try
        % the next higher value of whichsteptype.  Finally,
        % 'break' with length(xsol)>0 means convergence attained.
        
        if only_ag && whichsteptype < 3
            continue  % continue to whichsteptype=3
        end
        
        % Restart CG every 6*n+1 iterations.
        if whichsteptype == 1 && consec_cg >= 6*n + 1
            continue  % continue to whichsteptype=2
        end
        
        if whichsteptype == 2
            p = -g;
            consec_cg = 0;
        end
        
        
        if whichsteptype <= 2
            % This block of code attempts a CG iteration.
            % If success, then success is true and xnew is set.  If failure,
            % then whichsteptype is incremented by 1.
            consec_ag = 0;
            consec_cg = consec_cg + 1;
            
            if consec_cg == 1
                
                % On the first CG iteration after restarting, re-estimate L
                % if necessary.
                if Lnan && it > 0
                    Lnew = estimateL1(x,f,g,fg,L);
                    if ~isnan(Lnew)
                        L = Lnew;
                    end
                end
                
                % If using the more elaborate progress measure, then must
                % initialize variables associated with progress on the
                % first CG iteration.  Note: use the more elaborate test
                % only if v is not equal to x.
                use_elab = norm(v - x) > 0 && algparams.use_elab_progr;
                if use_elab
                    wz.w = v - x;
                    [~,gv] = fg(v);
                    wz.Aw = gv - g;
                    wz.baralpha = (g' * wz.w) / (wz.w' * wz.Aw);
                    xbar = x - wz.baralpha * wz.w;
                    [fbar, gbar] = fg(xbar);
                else
                    xbar = x;
                    fbar = f;
                    gbar = g;
                end
            end
            
            % Take a CG step.  These statements compute alphacg,
            % the CG line-search step length.
            
            xtilde = x + p/L;
            [~, gtilde] = fg(xtilde);
            if norm(gtilde) < algparams.gtol
                xsol = xtilde;
                break
            end
            Aoldp = Ap;
            Ap = L * (gtilde - g);
            pAp = p' * Ap;
            gp = g' * p;
            normg = norm(g);
            if pAp <= 0 || gp >= -algparams.bp_c2 * normg^2 || ...
                    gp <= -algparams.bp_c3 * normg^2
                % Failures for CG step if g'*p >= -c2*norm(g)^2 (i.e., p 
                % not a sufficient descent direction) or excessive descent                 
                % or if p'*Ap <= 0 (i.e., local quadratic
                % estimate for f has failed).  In this case, try another
                % kind of step (restarted CG or AG).
                continue
            end
            % May also restart if other Beale-Powell condition not
            % satisfied.
            if consec_cg > 1 && ~isinf(algparams.bp_c1) && ...
               abs(g'* prevg) >= algparams.bp_c1 * normg^2
                   continue
            end                
            prevg = g;
            alphacg = -gp / pAp;
            xnew = x + alphacg * p;  % x_{k+1}
            [fnew, gnew] = fg(xnew);
            if norm(gnew) < algparams.gtol
                xsol = xnew;
                break
            end
            
            [vnew, phistarnew] = ...
                compute_v_phistar(theta, gamma, gammanew, ell, ...
                v, phistar, xbar, fbar, gbar);
            
            wznew = wz;
            if use_elab
                sigma = (1 - theta) * gamma / gammanew;
                tau = theta / gammanew;
                if consec_cg == 1
                    Awz = wz.baralpha * wz.Aw;
                    xtmp = xnew + Awz / L;
                    [~,gtmp] = fg(xtmp);
                    AAwz = L * (gtmp - gnew);
                    wnew2 = (sigma - (1 - sigma) * wz.baralpha) * wz.w + tau * Awz;
                    Awnew2 = (sigma - (1 - sigma) * wz.baralpha) * wz.Aw + tau * AAwz;
                    znew2 = -wz.baralpha * wz.w + Awz / L;
                    Aznew2 = -wz.baralpha * wz.Aw + AAwz / L;
                else
                    Awz = wz.baralpha * wz.Aw + wz.barbeta * wz.Az;
                    xtmp = xnew + Awz / L;
                    [~,gtmp] = fg(xtmp);
                    AAwz = L * (gtmp - gnew);
                    wnew1 = (sigma - (1 - sigma) * wz.baralpha) * wz.w - ...
                        (1 - sigma) * wz.barbeta * wz.z + ...
                        tau * Awz;
                    Awnew1 = (sigma - (1 - sigma) * wz.baralpha) * wz.Aw - ...
                        (1 - sigma) * wz.barbeta * wz.Az + ...
                        tau * AAwz;
                    znew1 = -wz.baralpha * wz.w - wz.barbeta * wz.z + Awz / L;
                    Aznew1 = -wz.baralpha * wz.Aw - wz.barbeta * wz.Az + AAwz / L;
                    delta1 = (wnew1' * Aoldp) / (oldp' * Aoldp);
                    wnew2 = wnew1 - delta1 * oldp;
                    Awnew2 = Awnew1 - delta1 * Aoldp;
                    eps1 = (znew1' * Aoldp) / (oldp' * Aoldp);
                    znew2 = znew1 - eps1 * oldp;
                    Aznew2 = Aznew1 - eps1 * Aoldp;
                end
                delta2 = (wnew2' * Ap) / (p' * Ap);
                wznew.w = wnew2 - delta2 * p;
                wznew.Aw = Awnew2 - delta2 * Ap;
                eps2 = (znew2' * Ap) / (p' * Ap);
                znew3 = znew2 - eps2 * p;
                Aznew3 = Aznew2 - eps2 * Ap;
                mu = (znew3' * wznew.Aw) / (wznew.w' * wznew.Aw);
                wznew.z = znew3 - mu * wznew.w;
                wznew.Az = Aznew3 - mu * wznew.Aw;
                wznew.baralpha = (gnew' * wznew.w) / (wznew.w' * wznew.Aw);
                wznew.barbeta = (gnew' * wznew.z) / (wznew.z' * wznew.Az);
                xbarnew = xnew - wznew.baralpha * wznew.w - ...
                    wznew.barbeta * wznew.z;
                [fbarnew,gbarnew] = fg(xbarnew);
                if norm(gbarnew) < algparams.gtol
                    xsol = xbarnew;
                    break
                end
            else
                xbarnew = xnew;
                fbarnew = fnew;
                gbarnew = gnew;
            end
            if fbarnew <= phistarnew || fnew <= phistarnew
                % Successful CG step; compute new p and go to next iteration.
                [~,pnew] = hager_zhang_beta(p, g, gnew, norminitg);
                break
            end 
        else  % whichsteptype == 3 (AG step)
            
            if ~only_ag  %If this is the first AG step in a sequence, then
            % we set the only_ag flag.
                consec_ag = 0;
                only_ag = true;
            end   
            
            consec_ag = consec_ag + 1;
            ag_it = ag_it + 1;
            consec_cg = 0;
            
            % Use y_k := nesterov formula
            yy = (theta * gamma * v + gammanew * x) / (gamma + theta * ell); 
            [fyy, gyy] = fg(yy);
            if norm(gyy) < algparams.gtol
                xsol = yy;
                break
            end    
            if Lnan
                [L,fa,ga] = estimateL1(yy, fyy, gyy, fg, L);
            end
            % Test whether to terminate ag; this test happens only
            % when the ag iteration count is divisible by
            % algparams.check_ag_term.
            checkterm = rem(consec_ag, algparams.check_ag_term) == 0;
            xnew = yy - gyy/L;
            if ~Lnan && checkterm
                % If checking termination, need an extra function
                % evaluation.  This evaluation was already done a few lines
                % ago in the cases that Lnan holds.
                [fa,ga] = fg(xnew);
            end                       
            [vnew, phistarnew] = ...
                compute_v_phistar(theta, gamma, gammanew, ell, ...
                v, phistar, yy, fyy, gyy); %v_{k+1}, phistar_{k+1}                         
            if checkterm 
                q = -gyy'/(2*L) * (ga + gyy);
                if fa - fyy <= algparams.switch_fac * q
                    only_ag = false;
                    [fnew,gnew] = fg(xnew);
                    pnew = -gnew;
                    xbarnew = xnew;
                    fbarnew = fnew;
                    gbarnew = gnew;
                end    
            end
         end
    end  
    if ~isempty(xsol) %Test for convergence
        break
    end    
    gamma = gammanew;
    x = xnew;
    f = fnew;
    g = gnew;
    xbar = xbarnew;
    fbar = fbarnew;
    gbar = gbarnew;
    oldp = p;
    p = pnew;
    v = vnew;
    wz = wznew;
    phistar = phistarnew;
end

if isempty(xsol)
    x = xnew;
    it = algparams.maxit;
    success = false;
else
    x = xsol;
    success = true;
end
it = it + 1;
% Indicate in algorithm name whether BP parameters are default or not.
if isinf(algparams.bp_c1) && algparams.bp_c2 == 0 && isinf(algparams.bp_c3)
    bp = 'D';
else
    bp = '+';
end
stats = struct('numit', it, ...
    'success', success,...
    'algname', sprintf('C+AG(HaveL=%s,u=%s,s=%f,c=%d,BP=%s)',...
    string(~Lnan), ...
    string(algparams.use_elab_progr), ...
    algparams.switch_fac, algparams.check_ag_term, bp), ...
    'other', sprintf('num_ag=%d', ag_it));
end


function [theta, gammanew] = ...
    compute_theta_gamma(L, gamma, ell) %theta_{k-1}, gamma_{k}
% Compute theta and the next value of gamma according to 
% Nesterov's formulas
AA = L;
BB = gamma - ell;
CC = -gamma;
discr = BB^2 - 4 * AA * CC;
theta = (2 * CC) / (-BB - sqrt(discr));
gammanew = (1 - theta) * gamma + theta * ell;  % gamma_{k}
end

function [v, phistar] = ...
    compute_v_phistar(theta, gamma, newgamma, ell, v, phistar, yy, fyy, gyy) %v_k, phistar_k
% Compute the update to v and phistar according to Nesterov's formulas
oldv = v; % v_{k-1}
v = (1/newgamma) * ((1 - theta) * gamma * v + theta * ell * yy ...
    - theta * gyy);  % v_k
phistar = (1-theta) * phistar + theta * fyy - theta^2/(2*newgamma) * norm(gyy)^2 + ...
    (theta*(1-theta)*gamma/newgamma) * (ell/2 * norm(yy - oldv)^2 + ...
    gyy' * (oldv - yy)); % phistar_k
end

function [betahz,newp] = hager_zhang_beta(p, g, newg, norminitg)
% Compute the new search direction p according to Hager-Zhang's formulas.
haty = newg - g;
assert(norm(haty) > 0)
hatyp = haty' * p;
betahz1 = ((haty - p * (2 * norm(haty)^2 / hatyp))' * newg) / hatyp;
betahz2 = -1 / (norm(p) * min(.01*norminitg, norm(newg)));
betahz = max(betahz1, betahz2);
if isnan(betahz) || isinf(betahz)
    betahz = 0;
end
newp = -newg + betahz * p;
end
