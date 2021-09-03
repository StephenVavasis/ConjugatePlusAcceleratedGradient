function do_one_problem_all_algs(whichtest, n, params, gtol, ...
    maxits, outfilename)
% do_one_problem_all_algs(whichtest, n, paramrange, maxits, outfilename_
% This m-file runs four algorithms (or a subset of the four) on a single
% problem.
%
% Input arguments:
%   whichtest: either 'b' for ABPDN, 'L' for logistic, or 'u' for Huber.
%    (see paper)
%   n: problem size
%   params: All three problems need one or more parameter.  See comments below.
%      This is a vector of parameter values.
%   gtol: convergence tolerance
%   maxits:  A vector of four nonnegative integers.  This is the maximum
%      number of iterations for each algorithm.  The four algorithms in 
%      order are:
%        (1) C+AG, with z-test enabled (see paper; z-test is a more complicated
%            formula for xbar in the case that CG restarts after AG
%            iterations)
%        (2) C+AG with z-test disabled
%        (3) AG
%        (4) CG-Descent (Hager & Zhang)
%      To omit one of the four algorithms, set the corresponding entry of
%      maxits to 0.
%   outfilename: a text file that holds the result of all runs



if whichtest == 'b'
    % lambda is the coefficient of the regularization term.
    % delta is the smoothing parameter for the approximation to the abs
    % value
    lambda = params(1);
    delta = params(2);
    m = round(sqrt(n));
    k = m;
    while 1
        p0 = primes(k);
        if length(p0) >= m
            break
        end
        k = k * 2;
    end
    dctrows = p0(1:m);
    b = sin((1:m).^2)';
    handles = makeorthol1test(dctrows, n, b, lambda, delta);
elseif whichtest == 'L'
    % lambda is the coefficient of the regularization term
    % meandist is the mean distance of the points to the supporting
    % hyperplane
    % sigma is the standard deviation
    lambda = params(1);
    meandist = params(2);
    sigma = params(3);
    m = 0;
    d = round(n/2);
    A = zeros(n,d);
    x = meandist * ones(1,d) / sqrt(d);
    rng(19);
    for j = 1 : n
        A(j,:) = x + sigma * randn(1,d);
    end
    handles = makehandles_logistic(A, lambda);
else  
    assert(whichtest == 'u')
    % tau is the cutoff parameter of the Huber function
    tau = params(1);  
    is = zeros(2*n,1);
    js = zeros(2*n,1);
    es = zeros(2*n,1);
    is(1:2:2*n) = 1:n;
    is(2:2:2*n) = 2:n+1;
    js(1:2:2*n) = 1:n;
    js(2:2:2*n) = 1:n;
    es(1:2:2*n) = 1;
    es(2:2:2*n) = -1;
    A = sparse(is,js,es,n+1,n);
    m = n+1;
    b = (1:m)';
    handles = make_huber_regression(A,b,tau);
end

x0 = handles.goodstart();
outfile = fopen(outfilename, 'a');    
fprintf(outfile, '----------\n');
fprintf(outfile, 'begin at %s\n', datetime);
fprintf(outfile, 'code = %s  whichtest = %s m = %d n = %d gtol = %e\n', ...
    mfilename, whichtest, m, n, gtol);
fprintf(outfile,  'L = %e  ell = %e\n', handles.L, handles.ell);
fprintf(outfile, '  params = ');
for i = 1 : length(params)
    fprintf(outfile, '%e ', params(i));
end
fprintf(outfile, '\n');
for alg = 1 : 4
    maxit = maxits(alg);
    if maxit == 0
        continue
    end
    use_z = false;
    if alg == 1
        whichalg = 'CPLUSAG/NO_Z';
    elseif alg == 2
        whichalg = 'CPLUSAG/Z';
        use_z = true;
    elseif alg == 3
        whichalg = 'AG';
    else
        whichalg = 'CGDESCENT';
    end
    algparams = struct('gtol', gtol, 'maxit', maxit, ...
        'use_restart_z', use_z);
    
    fprintf('Alg %d=%s\n', alg, whichalg)
    fprintf(outfile, '++++++Alg = %s', whichalg);
    % reset the count of function/gradient evaluations to 0
    handles.resetcounts()
    if alg == 4
        % Convergence tolerance for CG_DESCENT uses the infinity-norm
        % instead of the 2-norm.  Experiments indicate that 
        % when the infinity-norm tolerance is set to 3*gtol/sqrt(n)
        % this yields a 2-norm of approximately gtol.        
        [x,status, stats] = cg_descent(x0, 3*gtol/sqrt(n), handles.f, ...
            handles.g, handles.fg);
        if status ~= 0
            fprintf(outfile, '   CG_DESCENT failed to converge\n');
        else
            fprintf(outfile, '   CG_DESCENT reports: numit = %d  numf = %d  numg = %d\n', ...
                stats.iter, stats.nfunc, stats.ngrad);
            numit = stats.iter;
        end
    elseif alg == 3
        [x, stats] = ag(handles, x0, algparams);
        numit = stats.numit;
    else
        [x, stats] = cplusag(handles, x0, algparams);
        numit = stats.numit;
    end
    fprintf(outfile,  '  numit = %d maxit = %d err/gtol = %f \n', ...
        numit, maxit, norm(handles.g(x))/gtol);
    if strcmpi(whichalg(1:min(7,length(whichalg))), 'CPLUSAG')
        fprintf(outfile, '  number of AG steps = %d\n', ...
            stats.numag);
    end
    counts = handles.getcounts();
    fprintf(outfile, '   fcount = %d fgcount = %d gcount = %d\n', ...
        counts.fcount, counts.fgcount, counts.gcount);
    fprintf(outfile, '   time is now %s\n', datetime);
    
end
fclose(outfile);
end
