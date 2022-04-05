function F = FT(A)
if (length(size(squeeze(A))) >= 3)
    F = fftshift( fftn( fftshift( A ) ) );
else
    F = fftshift( fft2( fftshift( A ) ) );
end