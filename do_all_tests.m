function do_all_tests(whichsize, testparams)
% This function tests all optimizers on eight standard problems for
% several sizes and parameter ranges.  The results are printed in the
% command window and saved to a file.
% Arguments 
%  whichsize should be set to either 'quick' (smaller sized
%    problems; should run in a few minutes on a modern laptop) or 'full'
%    (larger sized problems that take several hours to complete).
%  testparams is logical (either true or false) and indicates whether
%    fuller tests exercising all parameters of C+AG should be tried.
%    Default value is true if whichsize is 'quick' and false otherwise.

assert(strcmp(whichsize,'quick') || strcmp(whichsize,'full'))
quick = strcmp(whichsize,'quick');
if nargin < 2
    if quick
        testparams = true;
    else
        testparams = false;
    end
end


outfilename = [whichsize,'tests1a.out'];
outfile = fopen(outfilename, 'w');
fclose(outfile);
% The following values are used for all test runs
if quick
    maxfg = 100000;
else
    maxfg = 1000000;
end
algs = {
    struct('method', '+', 'useL', false, ...
      'params', struct( 'maxfg', maxfg, 'use_elab_progr', false))
    struct('method', '+', 'useL', false, ...
      'params', struct('maxfg', maxfg, 'use_elab_progr', true))
    struct('method', '+', 'useL', true, ...
      'params', struct('maxfg', maxfg, 'use_elab_progr', false))
    struct('method', '+', 'useL', true, ...
      'params', struct('maxfg', maxfg, 'use_elab_progr', true))
    struct('method', '+', 'useL', false, ...
      'params', struct( 'maxfg', maxfg, 'use_elab_progr', false, ...
      'bp_c1', 0.2, 'bp_c2', 0.8, 'bp_c3', 1.2))
    struct('method', 'a', 'useL', true, ...
      'params', struct('maxfg', maxfg))
    struct('method', 'a', 'useL', false, ...
      'params', struct('maxfg', maxfg))
    struct('method', 'c', 'useL', false, ...
      'params', struct())
    };
if quick
    n1 = 64;
    n2 = 200;
    n3 = 1000;
else
    n1 = 65536;
    n2 = 6000;
    n3 = 10000;    
end
if ~testparams
    algs = algs([1,6,7,8]);
end
for probnum = 1 : 9
    if probnum == 1
        gtol = 1e-8;
        whichprob = struct('family', 'q', 'n', 1000);
    elseif probnum == 2    
        gtol = 1e-8;
        whichprob = struct('family', 'b', 'n', n1, 'lambda', 1e-3, ...
            'delta', 1e-4);
    elseif probnum == 3
        gtol = 1e-8;
        whichprob = struct('family', 'b', 'n', n1, 'lambda', 1e-3, ...
            'delta', 5e-6);
    elseif probnum == 4   
        gtol = 1e-8;
        whichprob = struct('family', 'b', 'n', n1*4, 'lambda', 1e-3, ...
            'delta', 1e-4);
    elseif probnum == 5
        gtol = 1e-8;     
        whichprob = struct('family', 'b', 'n', n1*4, 'lambda', 1e-3, ...
            'delta', 5e-6);
    elseif probnum == 6
        gtol = 1e-8;
        whichprob = struct('family', 'L', 'n', n2, 'lambda', 1e-4, ...
            'meandist', 1, 'sigma', .4);
    elseif probnum == 7
        gtol = 1e-8;
        whichprob = struct('family', 'L', 'n', n2, 'lambda', 5e-6, ...
            'meandist', 1, 'sigma', .4);
    elseif probnum == 8
        gtol = 1e-6;
        whichprob = struct('family', 'u', 'n', n3, 'tau', 250);
    elseif probnum == 9
        gtol = 1e-6;
        whichprob = struct('family', 'u', 'n', n3, 'tau', 1000);
    end
    [handles, whichtest_str, ~, postprocess] = ...
        set_up_standard_problem(whichprob);
    do_one_problem_all_algs(handles, whichtest_str, gtol, algs, ...
        postprocess, outfilename);
end
end

