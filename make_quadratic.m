function handles = make_quadratic(A,b,computeL)

fcount = 0;
gcount = 0;
bothfgcount = 0;
    function f = func(x)
        f = x' * A * x / 2 - b' * x;
        fcount = fcount + 1;
    end

    function g = grad(x)
        g = A * x - b;
        gcount = gcount + 1;
    end

    function [f,g] = funcgrad(x)
        Ax = A * x;
        f = x' * Ax / 2 - b' * x;
        g = Ax - b;
        bothfgcount = bothfgcount + 1;
    end



    function resetcounts()
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
        

A = (A + A') / 2;
if computeL
    ee = eig(A);
    ell = min(ee);
    L = max(ee);
    assert(L > 0)
    assert(ell >= -L * 1e-15)
    ell = max(ell, 0);
else
    L = nan(1);
    ell = nan(1);
end


n = size(A,1);
handles = struct('f',@func, ...
    'g', @grad, ...
    'fg', @funcgrad, ...
    'L', L, ...
    'ell', ell, ...
    'n', n, ...
    'goodstart', zeros(n,1), ...
   'getcounts', @getcounts, ...
   'fgcount', @fgcount, ...
    'resetcounts', @resetcounts);
end

