%% TUTORIAL 2 CODE=========================================================
%==========================================================================
%==========================================================================
clear all
close all
%% this part is relevant for ASSIGNMENT 1==================================
%==========================================================================
% Phantom creation --------------------------------------------------------
phntm = phantom(); % Shepp-Logan phantom; you can use it for your own phantom,
% but experiment with the parameters

figure; % visualize the phantom
imagesc(phntm);
axis equal tight;
colormap gray;
colorbar
title('Shepp-Logan Phantom')
%% Generate sinogram ------------------------------------------------------
theta = 0:1:179; % specify projection angles
sinogram = radon(phntm, theta); % create a sinogram for the specified angles;
figure; % visualize the sinogram
imagesc(sinogram);
title('Sinogram')
xlabel('angles')
colormap gray;
colorbar
%% Reconstruct using (un)filtered packprojection --------------------------
recon = iradon(sinogram, theta, 'linear', 'none'); % reconstruct from the
%sinogram using unfiltered backproj.
recon_filtered = iradon(sinogram, theta, 'Ram-Lak'); % reconstruct using ramp filter
figure; % visualize reconstruction results
subplot(1,2,1); imagesc(recon); axis equal tight; colormap gray;
title('Unfiltered BP')
subplot(1,2,2); imagesc(recon_filtered); axis equal tight; colormap gray;
title('Filtered BP')
%% Reconstructing noisy data ----------------------------------------------
sinogram_noisy = awgn(sinogram, 30, 'measured'); % add zero-mean gaussian noiseto projections (30 dB SNR)
figure; % plot noisy and noiseless projection taken at 0 degrees
plot(sinogram(:,1)); hold on
plot(sinogram_noisy(:,1)); hold off
% now reconstruct the noisy data with various filters:
recon = iradon(sinogram_noisy, theta, 'linear', 'none'); % unfiltered backproj
% after fft, ramp filter multiplies every frequency w by its absolute value |w|.
% Hence, low frequencies are dampened and high frequencies are amplified
recon_filtered_ramp = iradon(sinogram_noisy, theta, 'linear', 'Ram-Lak'); %backproj with the ramp filter
% hamming filter can be obtained by element-wise multiplication of the rampfilter by the hamming window
% Low frequencies are dampened; mid. frequencies are amplified; high frequenciesare dampened
recon_filtered_ham = iradon(sinogram_noisy, theta, 'linear', 'Hamming'); %backproj with the hamming filter
figure; % visualize and compare the results
subplot(1,3,1); imagesc(recon); axis equal tight; colormap gray;
title('Unfiltered BP') % blurry
subplot(1,3,2); imagesc(recon_filtered_ramp); axis equal tight; colormap gray;
title('Filtered BP, ramp filter') % noisy, because high frequencies are amplified
subplot(1,3,3); imagesc(recon_filtered_ham); axis equal tight; colormap gray;
title('Filtered BP, cosine filter') % does a better job under noise
%% TO DO: imlement filtered backprojection with ramp, hamming and cosine filters
%==========================================================================
% you will need to complete the following code as a part of Assignment 1 in
% Homework 1. Alternatively, you can write your own code without using this
% template. DO NOT use iradon() in your implementation of the
% backprojection algorithm.
% see also: https://en.wikipedia.org/wiki/Fourier_transform
% how to use fftshift correctly:
%https://de.mathworks.com/matlabcentral/newsreader/view_thread/285244

% see slide 19 for visual representation of the filters

%ramp = ; % create a ramp filter

% NOTE: because we will be using naturally ordered fft, the DC component (or
% '0 freqiency') is in the CENTER of the fft. See slide 19 for more clues.
% NOTE 2: theory says that ramp filter should go as high as the maximum
% frequency. We will be using a scaled version with a maximum value of 1

%hmf = ; % create hamming filter by element-wise multiplication of the ramp
%filter by the hamming window
%cosf = ; % create cosine filter by element-wise multiplication of the ramp
%filter by the cosine window
% hint: you can get the cosine window by using cos(). Make sure it is
% symmetric and the values are in the [0, 1] range; check slide 19 for how
% it should look like
%recon = ; % initialize the reconstructed image with zeros
%for i = 1:length(theta) % cycle through projections (views)
%line = sinogram(:,i);
% see slide 16-17
%line= ; % 1. perform fft,
% 2. filter the fft by multiplying element-wise by the chosen filter
% 3. perform ifft
%image= ; % backproject; you may want to use repmat. Don't forget to rotate
%to the correct angle afterwards
%recon= ; % update your reconstruction; you may want to use + (yes that's a
%plus sign)
%end
%% This part is relevant for ASSIGNMENT 2 =================================
%==========================================================================
L = 1; % assume (arbitrary) resolution of the discretization
A = [L, L, 0, 0;... % set up the model matrix as shown in the slides
0, 0, L, L; ... % the model matrix relates the values of absorption mu =
                %[mu1, mu2, mu3, mu4];
L, 0, L, 0;... % to the measurements p = [P1, P2, P3, P4]' as in : A*mu = p
0, L, 0, L];
mu = rand([4,1]); % assume (random) values of absorption mu
p = A*mu; % simulate corresponding measurements b
% now we try to recover mu back from our measurements:
% A*mu = p => mu = inv(A)*p - we try to solve for mu that is assumed unknown 
%inv(A)*p % doesn't work, rank of A is 3!
mu_est = lsqr(A, p) % min||A*mu-p||.^2 - we try to find a solution using minmization procedure
abs(mu - mu_est)./mu*100 % looking at the relative error in percent, we're pretty far off from the real values
%% RECONSTRUCTION USING FOURIER SLICE THEOREM==============================
%=========================================================================
% now we go back to the CT reconstruction and demonstrate that Fourier
% slice theorem actually works. This is not used in any assignment.
proj_fft = fftshift( fft( ifftshift( sinogram, 1 ), [], 1 ), 1 ); % perform 1D fft of the projections - all at once
figure; % see how an fft of a projection looks like
subplot(1, 2, 1); plot(sinogram(:,1)); title('Projection')
subplot(1, 2, 2); plot(real(proj_fft(:,1))); title('FFT');
N = size(proj_fft, 1); % keep the length of the fft
%% create a radial grid for the FFT values, create the cartesian grid for The
FFT values
if mod(N,2) == 1
omega_r=(-N/2:(N-1)/2).*(2*pi/N); %define x axis of the FFT,
else
omega_r=(-(N-1)/2:(N-1)/2).*(2*pi/N); %define x axis of the FFT
end
omega_theta=theta*pi/180; % angles of the radial grid of FFT in rad
omega_xy=omega_r; % we keep the same resolution in Fx and Fy axes
% create a radial grid for the FFT values that we have - now we know their
% radial coordinates in the 2D Fourier space!
[theta_grid, r_grid] = meshgrid(omega_theta, omega_r);
% now create the cartesian grid for the 2D FFT values - we need the values
% of the 2D fft at those coordinates in order to compute inverse FFT and
% reconstruct the image
[omega_grid_x, omega_grid_y] = meshgrid(omega_xy, omega_xy);
% because the coordinates of the values that we obtained by taking the fft of
% the projections are polar, and the coordinates of the points that we need
% the values of are cartesian, we need to convert one of those sets of
% coordinates. I chose to convert cartesian to polar
[coord_th_fft2, coord_r_fft2] = cart2pol(omega_grid_x, omega_grid_y);
coord_r_fft2 = coord_r_fft2.*sign(coord_th_fft2); % if theta is negative or greater than pi
coord_th_fft2(coord_th_fft2<0)=coord_th_fft2(coord_th_fft2<0)+pi;
% now we have:
% proj_fft - values of the fourier transform on a radial grid
% theta_grid, r_grid - polar coordinates of those values
% coord_th_fft2, coord_r_fft2 - polar coordinates of the points that we
% need to know the values of in order to have a complete 2D FFT of our image
% we find the values we need by interpolation
FFT2 = interp2(theta_grid,r_grid,proj_fft,coord_th_fft2,coord_r_fft2,'bilinear',(0+1i.*0)); 
%interpolate coefficients that we have to the grid points Fx, Fy
FFT2 = flipud(FFT2);
I = fftshift(ifft2(ifftshift(FFT2))); % perform inverse FFT to get the image
pad = ceil( (N - size(phntm, 1))/2 ); % compute padding introduced by radon
figure;
subplot(1,2,1); imagesc( phntm ); colormap gray; axis equal tight
title('Original') % originall image
subplot(1,2,2); imagesc(abs(I(pad:end-pad, pad:end-pad))); colormap gray; axis
equal tight
title('Reconstructed') % cropped reconstructed image