%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Backprojection Algorithm %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
% sinogram: Input 2D - sinogram matrix
% theta: corresponding angles of sinogram
% filter_shape('none','Ramp','Cos','Hamming'): Selection of applied filter
% method
%
% Output: 
% rec: Reconstructed image

function rec = backprojection(sinogram,theta,filter_shape)

if nargin <3
    filter_shape = 'none';
end


% figure out how big the picture is going to be.
n = size(sinogram,1);
sideSize = n;

%filter setups
x = linspace(-1,1,sideSize);
ramp = abs(x);

switch filter_shape
    case 'Ramp'
        filter = ramp;
    case 'Hamming'
        filter = ramp .* hamming(sideSize)';
    case 'Cos'
        filter = ramp .* cos((x./2).*pi).^2;
    otherwise
        filter = 1;
end

% set up the image
rec = zeros(sideSize,sideSize);

for i = 1:length(theta)
    
    line = sinogram(:,i);
    
    % perform fft
    line_fft = fftshift(fft(ifftshift(line)));
        
    % filter in frquency space
    line_fft_filtered = line_fft .* filter';
    
    % transform back into time domain
    line_filtered = ifftshift(ifft(fftshift(line_fft_filtered)));
    
    % backproject
    image = repmat(line_filtered,1,sideSize);
    image = imrotate(image,theta(i),'crop');
    
    % sum up final picture
    rec = rec + image;
end

% rotate image as the original
rec = real(imrotate(rec,90))./length(theta);
end