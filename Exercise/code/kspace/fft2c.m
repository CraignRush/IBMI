function im = fft2c(d, N)
% im = fft2c(d)
% im = fft2c(d, N)
%
% fft2c performs a centered fft2
%
% d - matrix to fft
% N - matrix size of output (optional)

if (nargin < 2)
  im = fftshift(fft2(fftshift(d)));
else
  im = fftshift(fft2(fftshift(d), N, N));
end