function do_one_problem_all_algs(handles, whichtest_str, gtol, algs, ...
    postprocess, outfilename)
% do_one_problem_all_algs(whichtest, n, paramrange, whichalgs, outfilenams
% This m-file runs solver algorithms on a single problem.
%
% Input arguments:
% handles: objective function (see optimize.m)
% gtol: convergence tolerance (gradient 2-norm)
% algs: a cell array with algorithms to use.  See the documentation
%   of 'optimize.m' for interpretation.  An additional field of each entry 
%   of algs, 'useL', indicates whether to use the original L/ell returned by
%   'handles' or set them to Nan.
% postprocess: a function handle invoked on the solution vector to produce
%    a string with additional information.  Use @(x)('') if unneeded.
% outfilename: a text file that holds the result of all runs

n = handles.n;
if isfield(handles, 'goodstart')
    x0 = handles.goodstart();
else
    x0 = zeros(n,1);
end

outfile = fopen(outfilename, 'a');
fprintf(outfile, '----------\n');   
fprintf('----------\n');
fprintf(outfile, 'begin at %s\n', datetime);
fprintf(outfile, 'code = %s  whichtest = %s n = %d gtol = %e\n', ...
    mfilename, whichtest_str, n, gtol);
fprintf('code = %s  whichtest = %s n = %d gtol = %e\n', ...
    mfilename, whichtest_str, n, gtol);
fprintf(outfile, '\n');
origL = handles.L;
origell = handles.ell;
for whichalg = 1 : length(algs)
    handles.resetcounts()
    alg = algs{whichalg};
    if alg.useL
        handles.L = origL;
        handles.ell = origell;
    else
        handles.L = nan(1);
        handles.ell = nan(1);
    end
    [x, stat2] = optimize(handles, x0, gtol, alg);
    counts = handles.getcounts();
    out = sprintf(['%-30s, L=%e, ell=%e, success=%-5s numit=%7d %-12s\n', ...
        '   totalfgcount=%8d (%8d) norm(grad f(x^*))=%8.2e sol_info=%s\n'], ...
        stat2.algname, handles.L, handles.ell, ...
        string(stat2.success), stat2.numit, stat2.other, ...
        counts.fcount + counts.gcount + counts.bothfgcount, ...
        handles.fgcount(), ...
        norm(handles.g(x)), ...
        postprocess(x));
    if whichalg == 1
        firstx = x;
    else
        out = [out, sprintf('   norm(x-firstx) = %e\n', norm(x - firstx))];
    end
    fprintf('%s', out);
    fprintf(outfile, '%s', out);
    fprintf(outfile, '   time is now %s\n', datetime);  
end
fclose(outfile);
end
