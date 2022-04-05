function [Hn, dHn] = ButterworthTransferShort(s, n, fc)
Bn  = zeros(size(s));
dBn = zeros(size(s));
s   = 1i*s/fc;
an = [1	5.1258	13.1371	21.8462	25.6884	21.8462	13.1371	5.1258	1];
for k = 0:n
    Bn  = Bn  +   an(k+1)*s.^k;
    dBn = dBn + k*an(k+1)*s.^k;
end
dBn = -dBn/fc;
Kn  = 1./Bn;
dKn = -dBn/Bn.^2;
Hn  = abs(Kn);
dHn = real( Kn.*conj(dKn) )./Hn;
end