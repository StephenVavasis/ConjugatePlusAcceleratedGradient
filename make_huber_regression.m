function handles = make_huber_regression(A,b,tau)

fcount = 0;
gcount = 0;
fgcount = 0;
[m,n] = size(A);
ell = 0;
At = A';
L = 2*norm(At*A,1);

    function f = localf(x)
        res = A*x - b;
        mask1 = res < -tau;
        mask2 = res > tau;
        mask3 = ~(mask1 | mask2);
        f = sum(res(mask3).^2) - 2 * tau * sum(res(mask1)) + ...
            2 * tau * sum(res(mask2)) - (sum(mask1) + sum(mask2)) * tau^2;
        fcount = fcount + 1;
    end
    function [f,g] = localfg(x)
        res = A*x - b;
        mask1 = res < -tau;
        mask2 = res > tau;
        mask3 = ~(mask1 | mask2);
        f = sum(res(mask3).^2) - 2 * tau * sum(res(mask1)) + ...
            2 * tau * sum(res(mask2)) - (sum(mask1) + sum(mask2)) * tau^2;
        fderiv = zeros(m,1);
        fderiv(mask1) = -2 * tau;
        fderiv(mask2) = 2 * tau;
        fderiv(mask3) = 2 * res(mask3);
        g = At * fderiv;
        fgcount = fgcount + 1;
    end
    function g = localg(x)
        [~,g] = localfg(x);
        gcount = gcount + 1;
    end


    function reset_counts()
        fcount = 0;
        gcount = 0;
        fgcount = 0;
    end
    function s = getcounts()
        s = struct('fcount', fcount, 'gcount', gcount, ...
            'fgcount', fgcount);
    end
goodstart = @()(zeros(n,1));
handles =  struct('fg', @localfg, ...
    'f', @localf, ...   
    'g', @localg, ...
    'L', L, ...
    'ell', ell, ...
    'n', n, ...
    'goodstart', goodstart, ...
    'getcounts', @getcounts, ...
    'resetcounts', @reset_counts);
end


   