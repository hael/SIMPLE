function F = IFT(A)
if (length(size(squeeze(A))) >= 3)
    F = ifftshift( ifftn( ifftshift( A ) ) );
else
    F = ifftshift( ifft2( ifftshift( A ) ) );
end