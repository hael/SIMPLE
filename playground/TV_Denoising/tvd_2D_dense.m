% This implements the algorithm in
% 'On total-variation denoising: A new majorization-minimization
% algorithm and an experimental comparison with wavalet denoising.'
% M. Figueiredo, etc.
function ret = tvd_2D_dense(image, lambda, Niter)
[X, Y]    = size(image);
image_1D  = image(:);
N         = length(image_1D);
id_mat    = eye(N);
der_mat   = id_mat(2:N, :) - id_mat(1:N-1, :);
adjoint   = der_mat * der_mat';

cur   = image_1D;
der_x = der_mat*cur;
der_y = der_mat*image_1D;

for k = 1:Niter    
    banded_mat = zeros(N-1);
    for j = 1:N-1
        banded_mat(j,j) = abs(der_x(j))/lambda;
    end
    banded_mat = banded_mat + adjoint; % banded matrix
    cur        = image_1D - der_mat'*(banded_mat\der_y);            % solving the banded linear system
    der_x      = der_mat*cur;
    cost       = 0.5*sum(abs(cur - image_1D).^2) + lambda*sum(abs(der_x)); % current cost value
end
ret = reshape(cur, [X, Y]);
end