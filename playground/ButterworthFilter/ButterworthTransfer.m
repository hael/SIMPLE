function Hn = ButterworthTransfer(s, n, fc)
Bn = ones(size(s));
s  = 1i*s/fc;
for k = 0:(n/2-1)
    Bn = Bn.*(s.^2 - 2*s*cos(2*pi*(2*k+n+1)/(4*n)) + 1);
end
Hn = abs(1./Bn);
end