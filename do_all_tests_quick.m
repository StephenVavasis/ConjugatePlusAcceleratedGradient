function do_all_tests_quick()
% This function tests all four optimizers on all three problems for
% several sizes and parameter ranges.  The sizes are set so that
% the total run time is on the order of a few minutes.
outfilename = 'quicktests.out';
outfile = fopen(outfilename, 'w');
fclose(outfile);
% The following values are used for all test runs
gtol = 1e-8;
maxits = [50000, 50000, 100000, 50000];
% ABPDN tests
n1 = 64;
n2 = 256;
lambda = 1e-3;
delta1 = 1e-4;
delta2 = 5e-6;
params1 = [lambda,delta1];
params2 = [lambda,delta2];
disp('ABPDN 1')
do_one_problem_all_algs('b', n1, params1, gtol, maxits, outfilename)
disp('ABPDN 2')
do_one_problem_all_algs('b', n1, params2, gtol, maxits, outfilename)
disp('ABPDN 3')
do_one_problem_all_algs('b', n2, params1, gtol, maxits, outfilename)
disp('ABPDN 4')
do_one_problem_all_algs('b', n2, params2, gtol, maxits, outfilename)
% Logistic tests
lambda1 = 1e-4;
lambda2 = 5e-6;
meandist = 1;
sigma = .4;
params1 = [lambda1, meandist, sigma];
params2 = [lambda2, meandist, sigma];
n = 200;
disp('Logistic 1')
do_one_problem_all_algs('L', n, params1, gtol, maxits, outfilename)
disp('Logistic 2')
do_one_problem_all_algs('L', n, params2, gtol, maxits, outfilename)
% Huber tests
tau1 = 250;
tau2 = 1000;
params1 = tau1;
params2 = tau2;
n = 1000;
disp('Huber 1')
do_one_problem_all_algs('u', n, params1, gtol, maxits, outfilename)
disp('Huber 2')
do_one_problem_all_algs('u', n, params2, gtol, maxits, outfilename)
end

