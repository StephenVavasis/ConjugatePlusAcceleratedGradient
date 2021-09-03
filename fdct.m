function y = fdct(x)
N = length(x);
N2 = 2*N;
x2 = [x;-flipud(x)] .* exp(1i*pi*(0:2*N-1)'/(2*N)) * (exp(1i*pi/(4*N))/sqrt(N/2));
q = ifft(x2) * N;
y = real(q(1:N) .* exp(1i*pi * (0:N-1)'/(2*N)));
end

