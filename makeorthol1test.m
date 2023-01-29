function handles = makeorthol1test(dctrows, n, b, lambda, delta)
fcount = 0;
bothfgcount = 0;
gcount = 0;
    function resetcounts()
        fcount = 0;
        bothfgcount = 0;
        gcount = 0;
    end


    function s = getcounts()
        s = struct('fcount', fcount, 'gcount', gcount, ...
            'bothfgcount', bothfgcount);
    end
    function s = fgcount()
        s = fcount + gcount + bothfgcount;
    end

        

    function y = Ax(x)
        y0 = fdct(x);
        y = y0(dctrows);
    end
    


    function y = At(z)
        v = zeros(n,1);
        v(dctrows) = z;
        y = fdct(v);
    end

    function fval = localfuneval(y)
        fval = norm(Ax(y) - b)^2 / 2 + lambda * sum(sqrt(y.^2 + delta));
        fcount = fcount + 1;
    end

    function [fval,grad]=localfg(y)
        Ayb = Ax(y) - b;
        sqrty2 = sqrt(y.^2 + delta);
        fval = norm(Ayb)^2 / 2 + lambda * sum(sqrty2);
        grad = At(Ayb) + lambda * (y ./ sqrty2);
        bothfgcount = bothfgcount + 1;
    end

    function grad = localg(y)
        [~,grad] = localfg(y);
        gcount = gcount + 1;
    end

 
        


L = 1 + lambda/sqrt(delta);
%mu = lambda * delta  / norm(b,'inf');
%y1max = norm(b)^2 / (2 * lambda) + n * sqrt(delta);
%mu = lambda * delta / (y1max + delta)^3;
bfull = zeros(n,1);
bfull(dctrows) = b;
ystart = fdct(bfull);
btest = fdct(ystart);
assert(norm(btest(dctrows) - b) <= 1e-15 * norm(b))
y1obj = lambda * sum(sqrt(ystart.^2 + delta));
mxentry = sqrt(y1obj^2 / lambda^2 - delta);
mu = lambda * delta / (mxentry^2 + delta)^1.5;

    function y = goodstart()
        y = ystart;
    end

handles = struct('fg', @localfg, ...
    'f', @localfuneval, ... 
    'g', @localg, ...
    'L', L, ...
    'ell', mu, ...
    'n', n, ...
    'goodstart', @goodstart, ...
    'resetcounts', @resetcounts, ...
    'fgcount', @fgcount, ...
    'getcounts', @getcounts);

end
