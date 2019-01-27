function im = ifft2c(d, N)
% im = ifft2c(d)
% im = ifft2c(d, N)
%
% ifft2c performs a centered ifft2
%
% d - matrix to ifft
% N - matrix size of output (optional)

if (nargin < 2)
  im = fftshift(ifft2(fftshift(d)));
else
  im = fftshift(ifft2(fftshift(d), N, N));
end  
