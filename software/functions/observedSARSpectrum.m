function [sar_spectrum_lvl3] = observedSARSpectrum(plotsON,sar_transect, sar_transect_size, sar_sub_transect_size, sar_azimuth_resolution, k, size_of_filter_window, width_of_gaussian_lobe, cutoff_wavelength,kx,ky)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    
    dimension = sar_transect_size / sar_sub_transect_size;
    
    % Windowed 2D power spectrum estimation process:
    % Not exactly a filter, more a windowing function. It is a pre-FFT window that reduces spectral leakage and side lobes in the frequency domain.
    
    % Define the Hamming window
    window = hamming(sar_sub_transect_size) .* hamming(sar_sub_transect_size)';

    sar_spectrum_lvl1 = zeros(sar_sub_transect_size);
    for i = 1:dimension
        for j = 1:dimension
            
            filtered_image = window .* double(sar_transect((i-1)*128+1:i*128, (j-1)*128+1:j*128));
            
            average = mean(filtered_image(:));
      
            filtered_image = filtered_image / average - 1; % -1 to remove the new average which = 1 
          
            filtered_image = abs(fftshift(fft2(filtered_image))).^2; % Get the power spectral density [Monaldo, 1987]
            sar_spectrum_lvl1 = sar_spectrum_lvl1 + filtered_image;
        end
    end
    
    sar_spectrum_lvl1 = sar_spectrum_lvl1 ./ (dimension)^2; % Divide by the number of images summed
    sar_spectrum_lvl1 = (sar_azimuth_resolution / (2 * pi * sar_sub_transect_size))^2 * sar_spectrum_lvl1; % Convert from pixels to spatial frequency units
    
    if plotsON
        figure; plotLibrary().generalSpectrumPlots(0, sar_spectrum_lvl1, kx, ky, "After Mean Filter","?");
    end

    if (cutoff_wavelength ~= -1) % -1 condition means do not apply
        % Butterworth band pass filter: due to the symmetry of k in
        % both positive and negative directions, the filter appears as a
        % bandpass. Filters out all the wavenumbers below the cutoff
        % wavenumber
        k0 = (2 * pi) / cutoff_wavelength; % define the cutoff wavenumber

        butterworth = 1 ./ (1 + (k0 ./ k).^(2*10)); % Defined in Frequency Domain :)
        butterworth(isinf(butterworth)) = 0; butterworth(isnan(butterworth)) = 0;
        figure; contour(kx,ky,butterworth); title("Butterworth Response");

        sar_spectrum_lvl2 = butterworth .* sar_spectrum_lvl1;

        % figure; surf(k,k,butterworth); title("Butterworth filter for cutoff = 500 --> k0 = 2pi/500)"); xlabel("Wavenumber,k"); ylabel("Wavenumber,k");
        % figure; surf(butterworth); title("Butterworth filter response"); colorbar;
        % figure; surf(butterworth); title("Butterworth filter for cutoff = 500 --> k0 = 2pi/500)"); colorbar; hold on; surf(k); colorbar; hold off;
        % figure; surf(k); colorbar; title("Wave numbers"); c = colorbar();c.Label.String = 'Wave Number';
    else
        sar_spectrum_lvl2 = sar_spectrum_lvl1;
    end

    if plotsON
        figure; plotLibrary().generalSpectrumPlots(0,sar_spectrum_lvl2, kx, ky, "After Butterworth Filter","?");
    end
    
    % Smoothing with a Gaussian filter
    x = (0:size_of_filter_window) - (size_of_filter_window / 2); % create filter window coords centered around 0
    [xx, yy] = meshgrid(x, x); % create filter grid
    k_filt = exp(-(xx.^2 + yy.^2) ./ (width_of_gaussian_lobe)^2); % width_of_gaussian_lobe: Controls the spread of the Gaussian. A larger width means the Gaussian will be wider and smoother.
    k_filt = k_filt/sum(k_filt(:));

    sar_spectrum_lvl3 = conv2(sar_spectrum_lvl2, k_filt, 'same');
    
    
    if plotsON
       figure; contour(k_filt); title("Gaussian Response");
       figure; plotLibrary().generalSpectrumPlots(0,sar_spectrum_lvl3, kx, ky, "After Gaussian Filter","?");
    end
    
end