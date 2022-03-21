%%
X = 102;
bin_image = zeros(X,X);
max_r = 40;
min_r = 25;
for x = 1:X
    for y = 1:X
        dist = sqrt((x - 51)^2 + (y - 51)^2);
        if (dist < max_r && dist > min_r)
            bin_image(x,y) = 1;
        end
    end
end

% adding noises
noisy_image = imnoise(bin_image,'gaussian', 0, 0.01);

%%
figure;
subplot(141); imagesc(bin_image); colormap gray; axis image;
subplot(142); imagesc(noisy_image); colormap gray; axis image;

% doing TVD restoration
res_image = tvd_2D(noisy_image, 0.2, 100);
subplot(143); imagesc(res_image); colormap gray; axis image;

% doing TVD restoration with dense matrix
res_image = tvd_2D_dense(noisy_image, 0.2, 100);
subplot(144); imagesc(res_image); colormap gray; axis image;