%%
X = 202;
bin_image = zeros(X,X);
max_r = 70;
min_r = 35;
for x = 1:X
    for y = 1:X
        dist = sqrt((x - X/2)^2 + (y - X/2)^2);
        if (dist < max_r && dist > min_r)
            bin_image(x,y) = 1;
        end
    end
end

% adding noises
noisy_image = imnoise(bin_image,'gaussian', 0, 0.01);

%%
figure;
subplot(141); imagesc(bin_image);   colormap gray; axis image;
subplot(142); imagesc(noisy_image); colormap gray; axis image;

%%
x  = 1:300;
fc = 2;
Hn      = ButterworthTransfer(     x, 8, fc);
HnShort = ButterworthTransferShort(x, 8, fc);
Hn_rev = Hn(end:-1:1);
figure;  plot(abs(Hn));
hold on; plot(abs(HnShort));

%%
ButterKer = ButterworthKernel(202, 8, fc);
figure; plot(abs(ButterKer(101,:)));

ImgConv = abs(IFT(FT(ButterKer).*FT(bin_image)));
%ImgConv = ImgConv/max(ImgConv(:));
figure;
subplot(121); imagesc(noisy_image); colormap gray; axis image;
subplot(122); imagesc(ImgConv);     colormap gray; axis image;


%%
% options = optimoptions('fmincon',...
%     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',false);
% x0 = [7];
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = [];
% ub = [];
% [x,fval] = fmincon(@(x) L2NormCost(x, bin_image, ImgConv),x0,A,b,Aeq,beq,lb,ub,[],options);

l = -inf;
u =  inf;
opts    = struct( 'x0', 7);
opts.printEvery     = 1;
opts.m  = 5;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1e7;

[x,f,info] = lbfgsb( @(x) L2NormCost(x, bin_image, ImgConv) , l, u, opts );
%[x,f,info] = lbfgsb( @lbfgsb_test_func , l, u, opts );

    function [y, dy] = L2NormCost(x, bin_image, ImgConv)
        [ButterKer, d_ButterKer] = ButterworthKernel(202, 8, x);
        img       = abs(IFT(FT(ButterKer)  .*FT(bin_image)));
        df        = abs(IFT(FT(d_ButterKer).*FT(bin_image)));
        temp      = abs(img - ImgConv).^2;
        y         = sum(temp(:));
        diff      = (img - ImgConv).*df;
        dy        = 2*sum(diff(:));
        disp('x = ')
        disp(x)
    end