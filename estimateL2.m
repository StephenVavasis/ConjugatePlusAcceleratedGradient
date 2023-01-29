function L = estimateL2(x, f, g, fg, Linit)
itcount = 0;
ng = norm(g)^2;
L = Linit;
fa = fg(x - g / L);
while true
    itcount = itcount + 1;
    if itcount > 100
        error('Function may be unbounded below');
    end
    [fa, ~] = fg(x - g / L);
    if fa < f - ng / (2 * L)
        L = L / sqrt(2);
    else
        L = estimateL1(x, f, g, fg, L);
        return
    end
end
end