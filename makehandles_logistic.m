function handles = makehandles_logistic(A, lambda)

L = norm(A)^2 + lambda;
ell = lambda;
n = size(A, 2);

fcount = 0;
gcount = 0;
bothfgcount = 0;
    function f = fval(x, A, lambda)
        z = A * x;
        f = sum(log(1 + exp(-z))) + norm(x)^2 * lambda / 2;
        fcount = fcount + 1;
    end

    function [f,g] = fg(x, A, lambda)
        z = A * x;
        ee = exp(-z);
        f = sum(log(1 + ee)) + norm(x)^2 * lambda / 2;
        g = -A' * (ee ./ (1 + ee)) + lambda * x; 
        bothfgcount = bothfgcount + 1;
    end

    function g = gval(x, A, lambda)
        [~, g] = fg(x, A, lambda);
        gcount = gcount + 1;
    end


    function reset_counts()
        fcount = 0;
        gcount = 0;
        bothfgcount = 0;
    end

    function s = getcounts()
        s = struct('fcount', fcount, 'gcount', gcount, ...
            'bothfgcount', bothfgcount);
    end


    function s = fgcount()
        s = fcount + gcount + bothfgcount;
    end

localf = @(x)fval(x, A, lambda);
localg = @(x)gval(x, A, lambda);
localfg = @(x)fg(x, A, lambda);
goodstart = @()(zeros(n,1));
handles =  struct('fg', localfg, ...
    'f', localf, ...   
    'g', localg, ...
    'L', L, ...
    'ell', ell, ...
    'n', n, ...
    'goodstart', goodstart, ...
    'getcounts', @getcounts, ...
    'fgcount', @fgcount, ...
    'resetcounts', @reset_counts);
end
