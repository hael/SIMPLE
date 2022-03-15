%%
X = 31;
bin_image = zeros(X,X);
max_r = 10;
min_r = 7;
for x = 1:X
    for y = 1:X
        dist = sqrt((x - 16)^2 + (y - 16)^2);
        if (dist < max_r && dist > min_r)
            bin_image(x,y) = 1;
        end
    end
end

%%
figure;
subplot(131); imagesc(bin_image); colormap gray; axis image;

% adding noises
noisy_image = imnoise(bin_image,'gaussian', 0, 0.01);
subplot(132); imagesc(noisy_image); colormap gray; axis image;

% doing TVD restoration
res_image = tvd_2D(noisy_image, 0.2, 50);
subplot(133); imagesc(res_image); colormap gray; axis image;