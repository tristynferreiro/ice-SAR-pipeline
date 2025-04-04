function [sar_spectrum] = observedSARSpectrum(sar_transect, sar_transect_size, sar_sub_transect_size, sar_azimuth_resolution, k, size_of_filter_window, width_of_gaussian_lobe, cutoff_wavelength,first_guess_sar_spectrum, kx, ky)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    
    dimension = sar_transect_size / sar_sub_transect_size;
    
    % Define the Hamming window
    window = hamming(sar_sub_transect_size) .* hamming(sar_sub_transect_size)';
    
    sar_sub_transect = zeros(sar_sub_transect_size);
    for i = 1:dimension
        for j = 1:dimension
            
            filtered_image = window .* double(sar_transect((i-1)*128+1:i*128, (j-1)*128+1:j*128));
            
            average = mean(filtered_image(:));
      
            filtered_image = filtered_image / average - 1; % -1 to remove the new average which = 1 
          
            filtered_image = abs(fftshift(fft2(filtered_image))).^2; % Standard procedure (find in literature)
            sar_sub_transect = sar_sub_transect + filtered_image;
        end
    end

    sar_sub_transect = (sar_azimuth_resolution / (2 * pi * sar_sub_transect_size))^2 * sar_sub_transect; 
    % sar_sub_transect = sar_sub_transect ./ (dimension)^2;
    
    % figure; plotLibrary().generalSpectrumPlots(0,sar_sub_transect, kx, ky, "After Mean Filter");

    if (cutoff_wavelength ~= -1)
        % Butterworth high-pass filter: however, due to the symmetry of k in
        % both positive and negative directions, the filter appears as a
        % bandstop
        k(k == 0) = 1e-7; % to avoid division by zero
        k0 = (2 * pi) / cutoff_wavelength; % define the cutoff wavenumber

        % [temp1, temp2] = find(first_guess_sar_spectrum==max(first_guess_sar_spectrum(:)));
        % k0 = k(temp1(1),temp2(1)) / 2;

        butterworth = 1 ./ (1 + (k0 ./ k).^(2*10)); % Defined in Frequency Domain :)
        % figure; surf(butterworth); title("Butterworth Response");
        % anello  = 1 - 1./(1+(k0./modk).^5); % MIGHT NEED TO ADJUST POWER (5 in this case) FOR DIFFERENT CASES
        
        sar_spectrum = butterworth .* sar_sub_transect;

        % figure; surf(k,k,butterworth); title("Butterworth filter for cutoff = 500 --> k0 = 2pi/500)"); xlabel("Wavenumber,k"); ylabel("Wavenumber,k");
        % figure; surf(butterworth); title("Butterworth filter response"); colorbar;
        % figure; surf(butterworth); title("Butterworth filter for cutoff = 500 --> k0 = 2pi/500)"); colorbar; hold on; surf(k); colorbar; hold off;
        % figure; surf(k); colorbar; title("Wave numbers"); c = colorbar();c.Label.String = 'Wave Number';
    else
        sar_spectrum = sar_sub_transect;
    end

    % figure; plotLibrary().generalSpectrumPlots(0,sar_spectrum, kx, ky, "After Butterworth Filter");

    
    % Smoothing with a Gaussian filter
    x = (1:size_of_filter_window) - (size_of_filter_window) / 2; % create filter window coords centered around 0
    [xx, yy] = meshgrid(x, x); % create filter grid
    k_filt = exp(-(xx.^2 + yy.^2) / (width_of_gaussian_lobe)^2); % width_of_gaussian_lobe: Controls the spread of the Gaussian. A larger width means the Gaussian will be wider and smoother.
    k_filt = k_filt/sum(k_filt(:));
    % k_filt_large = zeros(128,128); % Initialize with zeros
    % k_filt_large(1:size_of_filter_window, 1:size_of_filter_window) = k_filt; % Place small filter in top-left
    % 
    % % Shift filter to center (for correct frequency response)
    % k_filt_large = fftshift(k_filt_large);
    % 
    % sar_spectrum = sar_spectrum .* k_filt_large;
    sar_spectrum = conv2(sar_spectrum, k_filt, 'same');
    
    % figure; surf(k_filt); title("Gaussian Response");
    % figure; plotLibrary().generalSpectrumPlots(0,sar_spectrum, kx, ky, "After Gaussian Filter");

    
end