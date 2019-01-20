%% Assignment 1
close all;
clear all;
%% 1) Create a phantom for simulating XCT measurements.
% The phantom should contain % ellipses and polygons of varying intensity
% (use phantom() and augment the resulting
% image). Show an image of your phantom.

% Size of quadratic Phantom in pixels
img_width = 256;
img_length = img_width;

% Determine Ellipe count (min = 3)
num_ellipses = 6;

% Initialize array
E = zeros(num_ellipses,6);

% Get basic skull structure from phantom
[~,base] = phantom('Modified Shepp-Logan', img_width);
E(1:2,:) = base(1:2,:);

% Get remaining array size
[M,N] = size(E(3:end,:));

% Create random ellipse along boundaries
E(3:end,1) = rand(M,1);
E(3:end,2) = -0.2 + (0.3+0.3)*rand(M,1);
E(3:end,3) = -0.4 + (0.4+0.4)*rand(M,1);
E(3:end,4) = -0.4 + (0.4+0.4)*rand(M,1);
E(3:end,5) = -0.45 + (0.45+0.45)*rand(M,1);
E(3:end,6) =  180*rand(M,1);
% Build phantom only of ellipses
phntm = phantom(E,img_width);

% Increase inensity of pixels in a square of the matrix
phntm(img_width/2-16:img_width/2+16,img_width/2-16:img_width/2+16) =...
    phntm(img_width/2-16:img_width/2+16,img_width/2-16:img_width/2+16) + 0.5;

% visualize the phantom
figure;
imagesc(phntm);
axis equal tight;
colormap gray;
colorbar
xlabel('x');
ylabel('y');
title('Random Shepp-Logan Phantom');
%% 2) Compute the views (projections) for the range of angles from 0° to 179° with spacing
% of 1°, 5° and 10°. Show the projections at 0°, 30°, 45° and 90° (in one axes). Show the
% sinogram with the most angles/projections.

% specify projection angles
theta = {0:1:179;
    0:5:179;
    0:10:179};

% pad the image with zeros so nothing gets lost during rotation
img_diag = sqrt(img_length^2 + img_width^2);
padding = ceil(img_diag - img_width) + 4;
pad_img = zeros(img_length + padding, img_width + padding);
pad_img(ceil(padding/2):(ceil(padding/2) + img_length - 1),...
    ceil(padding/2):(ceil(padding/2) + img_width - 1)) = phntm;

% loop over the number of angles and summarize
th = theta{1};
n = length(th);
img_pr = zeros(size(pad_img,2), n);
for i = 1:n
    tmp_img = imrotate(pad_img,180+th(i), 'bilinear', 'crop');
    img_pr(:,i) = (sum(tmp_img))';
end

% create a sinograms for the specified angles
for i = 1:length(theta)
    sinogram{i}(:,:) = radon(phntm, theta{i});
end

% Plot sinogram data at specific points in the same axes
figure;
plot(sinogram{1}(:,1),'DisplayName',[num2str(theta{1}(:,1))  ' degrees']);
hold on
plot(sinogram{1}(:,31),'DisplayName',[num2str(theta{1}(:,31)) ' degrees']);
plot(sinogram{1}(:,46),'DisplayName',[num2str(theta{1}(:,46)) ' degrees']);
plot(sinogram{1}(:,91),'DisplayName',[num2str(theta{1}(:,91)) ' degrees']);
title('Radon Transform at specific angles');
legend

% visualize the sinogram
figure;
imagesc(sinogram{1});
title(['Sinogram @ ' num2str(length(theta{1})) 'angles (' num2str(theta{1}(:,1)) '°-' num2str(theta{1}(:,end)) '°)']);
xlabel('Angle')
colormap gray;
colorbar

%% 3) Implement the backprojection algorithm (i.e. inverse radon transform) according to
% slide 16 of the Tutorial slides. You are not allowed to use iradon().

% see slide 19 for visual representation of the filters
% NOTE: because we will be using naturally ordered fft, the DC component (or
% '0 freqiency') is in the CENTER of the fft. See slide 19 for more clues.
% NOTE 2: theory says that ramp filter should go as high as the maximum
% frequency. We will be using a scaled version with a maximum value of 1
% hint: you can get the cosine window by using cos(). Make sure it is
% symmetric and the values are in the [0, 1] range; check slide 19 for how
% it should look like
% create cosine filter by element-wise multiplication of the ramp
%filter by the cosine window

% see backprojection.m

%% 4) Reconstruct the phantom data with the specified angular spacings using your
% backprojection algorithm without filtering. Show the obtained reconstructions.
recon = backprojection(sinogram{1},theta{1});
% visualize reconstruction results
figure;
subplot(1,2,1);
imagesc(phntm);
axis equal tight;
colormap gray;
title('Original'); xlabel('x'); ylabel('y');
subplot(1,2,2);
imagesc(recon);
axis equal tight;
colormap gray;
title('Unfiltered BP'); xlabel('x'); ylabel('y');
%% 5) Incorporate filtering in your backprojection. Implement 3 filters (ramp, cosine and
% hamming) and test their influence on the reconstruction using your phantom data (pick
% a single angle spacing). Show the reconstruction results
for i = 1:3
    % Compute filtered backprojections
    recon_cos = backprojection(sinogram{i},theta{i},'Cos');
    recon_ramp = backprojection(sinogram{i},theta{i},'Ramp');
    recon_hamming = backprojection(sinogram{i},theta{i},'Hamming');
    % visualize reconstruction results
    figure;
    subplot(2,2,1);
    imagesc(phntm);
    title(['w/o Filtering,' num2str(length(theta{i})) ' angles']);
    colormap gray; colorbar; xlabel('x'); ylabel('y');
    subplot(2,2,2);
    imagesc(recon_ramp);
    title(['w/ Ramp Filtering, ' num2str(length(theta{i})) ' angles']);
    colormap gray; colorbar; xlabel('x'); ylabel('y');
    subplot(2,2,3);
    imagesc(recon_hamming);
    title(['w/ Hamming Filtering, '  num2str(length(theta{i})) ' angles']);
    colormap gray; colorbar; xlabel('x'); ylabel('y');
    subplot(2,2,4);
    imagesc(recon_cos);
    title(['w/ Cosine Filtering, ' num2str(length(theta{i})) ' angles']);
    colormap gray; colorbar; xlabel('x'); ylabel('y');
end
%% 6) Reconstruct the provided datasets (CT_2018.mat) with your backprojection algorithm
% without filtering and with each of the implemented filters, respectively. Show the
% reconstructed images.

% load sample data
S = load('CT_2018.mat');
% loop over sample data, compute backprojections and plot them each in a
% single figure
for i = 1:3
    % change variable name dynamically
    sino_name = S.(['sino' num2str(i-1)]);
    angle = S.(['theta' num2str(i-1)]);
    % visualize
    fig_name = sprintf('%s with %d angles from %d° to %d°', ['sino' num2str(i-1)],length(angle), angle(1),angle(end));
    figure('Name',fig_name);
    subplot(2,2,1);
    recon = backprojection(sino_name,angle);
    imagesc(recon);
    title('w/o Filtering');
    colormap gray; colorbar; xlabel('x'); ylabel('y');
    subplot(2,2,2);
    recon = backprojection(sino_name,angle,'Ramp');
    imagesc(recon);
    title('w/ Ramp Filtering');
    colormap gray; colorbar; xlabel('x'); ylabel('y');
    subplot(2,2,3);
    recon = backprojection(sino_name,angle,'Hamming');
    imagesc(recon);
    title('w/ Hamming Filtering');
    colormap gray; colorbar; xlabel('x'); ylabel('y');
    subplot(2,2,4);
    recon = backprojection(sino_name,angle,'Cos');
    imagesc(recon);
    title('w/ Cosine Filtering');
    colormap gray; colorbar; xlabel('x'); ylabel('y');
end
disp('7) Shortly interpret your results.');
disp('a. What is the effect of different angular spacings on the reconstruction?');
disp('Finer Spacings <==> More Information, therefore better reconstrucion, more details visible');
disp('b. How do the different filters change the reconstruction results?');
disp('Filters affect the edge sharpness, the redsidual noise and the overall brightness of an image');
disp('c. Which filter performs best? Why? Under which circumstances?');
disp(['The Hamming filter performs best overall, because he does not damp the low frequencies to zero,' ...
    'cosine and ramp filter do.']);
