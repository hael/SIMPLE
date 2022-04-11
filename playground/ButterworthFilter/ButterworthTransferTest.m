%%
X = 202;
bin_image = zeros(X,X,X);
max_r = 70;
min_r = 35;
for x = 1:X
    for y = 1:X
        for z = 1:X
            dist = sqrt((x - X/2)^2 + (y - X/2)^2 + (z - X/2)^2);
            if (dist < max_r && dist > min_r)
                bin_image(x,y,z) = 1;
            end
        end
    end
end

% adding noises
noisy_image1 = imnoise(bin_image,'gaussian', 0, 0.01);
noisy_image2 = imnoise(bin_image,'gaussian', 0, 0.01);

%%
figure;
subplot(141); imagesc(bin_image(:,:,101));   colormap gray; axis image;
subplot(142); imagesc(noisy_image1(:,:,101)); colormap gray; axis image;

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
figure; plot(abs(ButterKer(101,:,101)));

ImgConv = abs(IFT(FT(ButterKer).*FT(bin_image)));
%ImgConv = ImgConv/max(ImgConv(:));
figure;
subplot(121); imagesc(noisy_image1(:,:,101)); colormap gray; axis image;
subplot(122); imagesc(ImgConv(:,:,101));     colormap gray; axis image;


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

l = 1;
u = 50;
opts    = struct( 'x0', 1);
opts.printEvery     = 1;
opts.m  = 5;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1e7;

[x,f,info] = lbfgsb( @(x) L2NormCost(x, noisy_image1, noisy_image2, ImgConv) , l, u, opts );
%[x,f,info] = lbfgsb( @lbfgsb_test_func , l, u, opts );


% plot the restoration
ButterKer = ButterworthKernel(202, 8, x);
ImgConv = abs(IFT(FT(ButterKer).*FT(noisy_image1)));
%ImgConv = ImgConv/max(ImgConv(:));
figure;
subplot(121); imagesc(noisy_image1(:,:,101)); colormap gray; axis image;
subplot(122); imagesc(ImgConv(:,:,101));     colormap gray; axis image;

    function [y, dy] = L2NormCost(x, image1, image2, ImgConv)
        [ButterKer, d_ButterKer] = ButterworthKernel(202, 8, x);
        img1      = abs(IFT(FT(ButterKer)  .*FT(image1)));
        df1       = abs(IFT(FT(d_ButterKer).*FT(image1)));
        img2      = abs(IFT(FT(ButterKer)  .*FT(image2)));
        df2       = abs(IFT(FT(d_ButterKer).*FT(image2)));
        %temp      = abs(img1 - image2).^2 + abs(img2 - image1).^2;
        temp      = abs(img1 - ImgConv).^2;
        y         = sum(temp(:));
        %diff      = (img1 - image2).*df1 + (img2 - image1).*df2;
        diff      = (img1 - ImgConv).*df1;
        dy        = 2*sum(diff(:));
        disp(max(d_ButterKer(:)))
        disp('x = ')
        disp(x)
        disp('cost = ')
        disp(y)
        disp('dcost = ')
        disp(dy)
    end