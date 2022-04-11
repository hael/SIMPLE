function [ker, d_ker] = ButterworthKernel(width, n, fc)
w_over_2 = width/2;
ker      = zeros(width,width,width);
d_ker    = zeros(width,width,width);
for k = 1:width
    for l = 1:width
        for j = 1:width
            dist     = sqrt((k-w_over_2).^2 + (l - w_over_2).^2 + (j - w_over_2).^2);
            [ker(k,l,j), d_ker(k,l,j)] = ButterworthTransferShort(dist, n, fc);
        end
    end
end
end