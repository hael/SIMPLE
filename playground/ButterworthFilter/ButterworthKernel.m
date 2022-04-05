function [ker, d_ker] = ButterworthKernel(width, n, fc)
w_over_2 = width/2;
ker      = zeros(width);
d_ker    = zeros(width);
for k = 1:width
    for l = 1:width
        dist     = sqrt((k-w_over_2).^2 + (l - w_over_2).^2);
        [ker(k,l), d_ker(k,l)] = ButterworthTransferShort(dist, n, fc);
    end
end
end