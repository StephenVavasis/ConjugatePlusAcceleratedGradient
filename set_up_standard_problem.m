function [handles, whichtest_str, n, postprocess] = set_up_standard_problem(whichprob)
% [handles, whichprob_str] = set_up_standard_problem(whichprob)
% Input argument:
%  whichprob is a struct whose fields are:
%   'family': either 'b' for ABPDN, 'L' for logistic, 'u' for Huber,
%   'q' for quadratic with spread-out eigenvalues, 'Q' for quadratic with
%       clustered eigenvalues. (see paper)
%   'n' problem size
%     and additional parameters depending on the family (see code below)
% Output argument:
%   handles: Handles for function and gradient evaluation.  
%      See optimize.m for more details. 
%   whichprob_str: a string containing a printable description of the
%     problem.
whichtest = whichprob.family;
n = whichprob.n;
postprocess = @(y)("none");
if whichtest == 'b'
    % lambda is the coefficient of the regularization term.
    % delta is the smoothing parameter for the approximation to the abs
    % value
    lambda = whichprob.lambda;
    delta = whichprob.delta;
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
    whichtest_str = sprintf('b(lambda=%e,delta=%e)', lambda, delta);
    postprocess = @(y)(sprintf("numappxzero=%d", sum(abs(y) <= sqrt(delta))));
elseif whichtest == 'L'
    % lambda is the coefficient of the regularization term
    % meandist is the mean distance of the points to the supporting
    % hyperplane
    % sigma is the standard deviation
    lambda = whichprob.lambda;
    meandist = whichprob.meandist;
    sigma = whichprob.sigma;
    d = round(n/2);
    A = zeros(n,d);
    x = meandist * ones(1,d) / sqrt(d);
    rng(19);
    for j = 1 : n
        A(j,:) = x + sigma * randn(1,d);
    end
    handles = makehandles_logistic(A, lambda);
    whichtest_str = sprintf('L(lambda=%e,meandist=%e,sigma=%e)', lambda, meandist, sigma);
elseif whichtest == 'u'
    % tau is the cutoff parameter of the Huber function
    tau = whichprob.tau;  
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
    b = 500*[ones(n,1);-n*1.1];
    handles = make_huber_regression(A,b,tau);
    whichtest_str = sprintf('u(tau=%e)', tau);
elseif whichtest == 'q'
    is = 1:n;
    js = 1:n;
    es = (1:n).^2;
    A = sparse(is,js,es,n,n);
    b = sin(1:n)';
    handles = make_quadratic(A,b,false);
    handles.L = n^2;
    handles.ell = 1;
    whichtest_str = 'q';
elseif whichtest == 'Q'
    is = 1:n;
    js = 1:n;
    es = zeros(1,n);
    % eigenvalues are 1..sqrt(n) rounded down (so many repeats) plus noise
    for i = 1 : n
        es(i) = floor(sqrt(i)) + 1e-5 * sin(i);
    end
    A = sparse(is,js,es,n,n);
    b = sin(1:n)';
    handles = make_quadratic(A,b,false);
    handles.L = n^2;
    handles.ell = 1;
    whichtest_str = 'Q';
else
    error('unrecognized problem family')
end

    


end
