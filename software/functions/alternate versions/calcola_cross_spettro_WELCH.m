% Alternate version that has been converted to MATLAB code from IDL
% Originally used for cross spectra cases but this can be used to
% calculate the observed SAR spectrum
function [cross_spec] = calcola_cross_spettro_WELCH(img1, img2, n, dx, dim_filt, width, kx, ky)

    n21 = n / 2 + 1;
    s = 2 * pi / (n * dx);
    tmp = -pi / dx + (0:n) * s;
    % kx = tmp(2:n+1);
    % ky = kx;
    
    
    % % Calculate the backscatter coefficient
    % sigma0 = mean(img1(:));
    % dev_sigma0 = std(img1(:));
    % fprintf('sigma0 = %.2f [dB]\n', 10 * log10(sigma0));
    % fprintf('dev sigma0 = %.2f [dB]\n', 10 * (dev_sigma0 / sigma0) / log(10));
    % 
    
    scale = 512 / n;
   
    % Define the Hamming window
    w = hamming(n) * hamming(n)';
    
    imfin = zeros(n, n);
    for i = 1:4
        for j = 1:4
            
            image1 = w .* double(img1((i-1)*128+1:i*128, (j-1)*128+1:j*128));
            image2 = w .* double(img2((i-1)*128+1:i*128, (j-1)*128+1:j*128));
            ave1 = mean(image1(:));
            ave2 = mean(image2(:));
            image1 = image1 / ave1 - 1; % -1 to remove the new average which = 1 
            image2 = image2 / ave2 - 1;
            % image = (fft2(image1) ) .* conj(fft2(image2));
            image = abs(fft2(image1)).^2; % Standard procedure (find in literature)
            image = fftshift(image);
            imfin = imfin + image;
        end
    end
    % imfin = ((dx / (2 * pi * n))^2) * imfin / (2 * scale - 1)^2;
    imfin = (dx / (2 * pi * n))^2 * imfin; 
    imfin = imfin / (scale)^2;
    % imfin = imfin/16;
    
    % % Estimate the raw cross-covariance function
    % cross_var_RAW = fftshift(ifft2(fftshift(imfin)));
    
    

    % Butterworth band-pass filter
    % [kr, ka] = meshgrid(kx, ky);
    kr = kx;
    ka = ky;
    k = sqrt(kr.^2 + ka.^2);
    k(k == 0) = 1e-7;
    cutoff = 500; % This is dependent on the sea state and used to remove the low frequencies. Maybe see if 2 * swell will work.
    k0 = 2 * pi / cutoff;
    butterworth = 1 ./ (1 + (k0 ./ k).^20);
    cross_spec = butterworth .* imfin;

    % Smoothing with a Gaussian filter
    in = dim_filt;
    x = (1:in) - (in) / 2;
    [nn, mm] = meshgrid(x, x);
    k_filt = exp(-(nn.^2 + mm.^2) / width^2);
    k_filt = k_filt/sum(k_filt(:));
    cross_spec = conv2(cross_spec, k_filt, 'same');
    
    % % Estimate the smoothed cross-covariance function
    % cross_var_SMOOTHED = fftshift(ifft2(fftshift(cross_spec)));
end