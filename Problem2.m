%% Problem 2

close all;
clear all;
load braint1data.mat;
%% reconstruct image using complex k-space data
ima = ifft2c(rawkspace);

figure;
imagesc(abs(ima));
colormap(gray);axis equal tight;
title('Reconstruction from k-space data');
xlabel('x');
ylabel('y ');

% using only magnitude of k-space
ima = ifft2c(abs(rawkspace));
figure;
imagesc(abs(ima));
colormap(gray);axis equal tight;
title('Reconstruction using only magnitude of k-space');
xlabel('x');
ylabel('y ');

% using only phase of k-space
ima = ifft2c(angle(rawkspace));
figure;
imagesc(abs(ima));
colormap(gray);axis equal tight;
title('Reconstruction using only phase of k-space');
xlabel('x');
ylabel('y ');
%% remove every other data point from k-space along x-axis
reduced_rawkspace = rawkspace;
reduced_rawkspace( :, 1:2:end)=[];

%reconstruct image from reduced k-space data
reduced_ima = ifft2c(reduced_rawkspace);
figure;
imagesc(abs(reduced_ima));
colormap(gray);axis equal tight;
title(['Reconstruction using only every other data point from k-space' ...
'along horizontal direction']);
xlabel('x');
ylabel('y ');

%% disturb pixel (160,160)
disturbed_kspace = rawkspace;
disturbed_kspace(160,160) = 10^4;
disturbed_ima = ifft2c(disturbed_kspace);

figure;
imagesc(abs(disturbed_ima));

colormap(gray);axis equal tight;
title('Reconstruction with disturbed pixel in k-space [160,160]');
xlabel('x');
ylabel('y ');
%% disturb column 160
disturbed_kspace2 = rawkspace;
disturbed_kspace2(:,160) = 10^2;
disturbed_ima2 = ifft2c(disturbed_kspace2);

figure;
imagesc(abs(disturbed_ima2));
colormap(gray);axis equal tight;
title('Reconstruction with distrubed column in k-space 160');
xlabel('x');
ylabel('y ');

%% SNR in signal
ima = ifft2c(rawkspace);

% calculate SNR from small window with no signal
yinn=1; youtn=50;
xinn=1; xoutn=50;
small_window = abs(ima(yinn:youtn,xinn:xoutn));
standard_deviation = std(std(small_window));
overall_SNR = mean(mean(abs(ima)))/standard_deviation

% % SNR measurements
% % original data
yins=110; youts=115; xins=90; xouts=95;
yinn=1; youtn=50; xinn=1; xoutn=50;

SNRorig = mean(mean(abs(ima(yins:youts,xins:xouts))))/...
    std(std(abs(ima(yinn:youtn,xinn:xoutn))))

%SNRtriang = mean(mean(abs(imtriang(yins:youts,xins:xouts))))/...
%    std(std(abs(imtriang(yinn:youtn,xinn:xoutn))))

% SNR has been defined as the ratio of the average signal value to the
% standard deviation
% average signal over 5x5 window
average_mask = fspecial('average', 5);
averaged_image = imfilter(abs(ima), average_mask);
% devide by standard_deviation
SNRorig = averaged_image/standard_deviation;

figure;
imagesc(SNRorig);
colormap(jet);axis equal tight;
colorbar
title('Image SNR over 5x5 window average');
xlabel('x');
ylabel('y ');
textColor = 'black';
textBackground = 'white';

text(60, 15, ...
['Overall SNR: ', num2str(overall_SNR) ], ...
'Color', textColor, ...
'BackgroundColor', textBackground, ...
'HorizontalAlignment', 'Center');
%% zero-filling to 512x512
% image 512
rawkspacefill = zeros(512,512);
rawkspacefill([1:256] + 128,[1:256] + 128) = rawkspace;
imfill = ifft2c(rawkspacefill);

% reconstruct the image
figure;
imagesc(abs(imfill));
colormap(gray);axis equal tight;
title('Reconstruction with zero spacing');
xlabel('x');
ylabel('y ');

% now SNR of image --> double size of window and filter
yinn=1; youtn=100;
xinn=1; xoutn=100;
small_window = abs(imfill(yinn:youtn,xinn:xoutn));
standard_deviation_2 = std(std(small_window));

%overall SNR
overall_SNR_imfill = mean(mean(abs(imfill)))/standard_deviation_2
average_mask = fspecial('average', 10);
averaged_imfill = imfilter(abs(imfill), average_mask);

% devide by standard_deviation
SNRimfill = averaged_imfill/standard_deviation_2;

figure;
imagesc(SNRimfill);
colormap(jet);axis equal tight;
colorbar
title('Zero spacing: Recontructed image SNR');
xlabel('x');
ylabel('y ');
textColor = 'black';
textBackground = 'white';
text(120, 25, ...
['Overall SNR: ', num2str(overall_SNR_imfill) ], ...
'Color', textColor, ...
'BackgroundColor', textBackground, ...
'HorizontalAlignment', 'Center');
%% Hanning window
% hanning window in k-space
rawkspace_hann = rawkspace.* (hann(256)* hann(256)');
ima_hann = ifft2c(rawkspace_hann);

%reconstruct image
figure;
imagesc(abs(ima_hann));
colormap(gray);axis equal tight;
title('Reconstruction with Hanning filter');
xlabel('x');
ylabel('y ');

% calculate SNR from small window with no signal
yinn=1; youtn=50;
xinn=1; xoutn=50;
small_window_hann = abs(ima_hann(yinn:youtn,xinn:xoutn));
standard_deviation_hanning = std(std(small_window_hann));
overall_SNR_hanning = mean(mean(abs(ima_hann)))/standard_deviation_hanning

% SNR has been defined as the ratio of the average signal value to the
% standard deviation
% average signal over 5x5 window
average_mask_hann = fspecial('average', 5);
averaged_image = imfilter(abs(ima_hann), average_mask);

% devide by standard_deviation
SNRhann = averaged_image/standard_deviation_hanning;

figure;
imagesc(SNRhann);
colormap(jet);axis equal tight;
colorbar
title('Hanning Filter: Reconstructed image SNR');
xlabel('x');
ylabel('y ');
textColor = 'black';
textBackground = 'white';
text(60, 15, ...
['Overall SNR: ', num2str(overall_SNR_hanning) ], ...
'Color', textColor, ...
'BackgroundColor', textBackground, ...
'HorizontalAlignment', 'Center');