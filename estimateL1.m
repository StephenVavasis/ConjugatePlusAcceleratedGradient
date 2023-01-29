function [L,fa,ga] = estimateL1(x, f, g, fg, Linit)
itcount = 0;
ng = norm(g)^2;
L = Linit;
while true
    itcount = itcount + 1;
    if itcount > 60
        L = nan(1);
        break
    end
    [fa, ga] = fg(x - g / L);
    if fa > f - ng / (2 * L) && abs(fa - f) > 1e-11 * abs(f) 
        L = L * sqrt(2);
    else
        break
    end
end
end