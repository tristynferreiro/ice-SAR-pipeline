function waves = waveFunctions
% Wave Function Library 
%   This contains all of the functions relating to waves

waves.waveNumberSpectrum = @waveNumberSpectrum; % E(f, direction) -> E(kx,ky)
waves.waveSpectrum = @waveSpectrum; % E(kx,ky) -> E(f, direction)

waves.interpolateWaveNumberSpectrum = @interpolateWaveNumberSpectrum;
waves.calculateSpectrumCharacteristics = @calculateSpectrumCharacteristics;
end

function [E_kx_ky, Jacobian,kx_matrix,ky_matrix] = waveNumberSpectrum(wave_spectrum,omega,k, direction_bins_degrees)
%WaveNumberSpectrum convert a wave spectrum E(frequency, direction) to k domain, E(kx,ky)
%   According to Holthuijsen: "Remote-sensing and numerical wave models can 
%   estimate the full two-dimensional spectrum, usually the wave-number 
%   spectrum E (k x , k y )" [pg.52 Holthuijsen] 
    
% Wave Number Spectrum E(kx,ky)
    % First convert to omega from frequency
    era5_d2wd = wave_spectrum ./ (2 * pi); % conver E(f, theta) to E(w, theta)

    % Wave phase speed (c_w) for deep water
    % era5_cw_wave_speed = gravity ./ (era5_omega); % c_0, [Eq.5.4.24 Holthuijsen]
    era5_cw_wave_speed = omega ./k; % c, [Eq.3.5.36 Holthuijsen] same as above
    
    % Group wave speed (c_g)
    era5_n_dispersion_scaling_factor = 0.5; % 0.5 is used for deep water, [Eq.5.4.32 Holthuijsen]
    era5_cg_group_wave_speed = era5_n_dispersion_scaling_factor * era5_cw_wave_speed; % [Eq.5.4.32 Holthuijsen]
    
    % Calculate the wave number spectrum E(kx,ky)
    Jacobian = ((era5_cw_wave_speed .* era5_cg_group_wave_speed) ./ omega ) ; % [Eq.3.5.36 Holthuijsen]
    E_kx_ky = Jacobian .* era5_d2wd; % [Eq.3.5.36 Holthuijsen]
    
% kx and ky variables
    % This is H&H convention:
    % era5_kx_matrix_hhdefintion = era5_k .* cosd(era5_direction_bins); % [rad/m] [Eq.3.5.19b Holthuijsen]
    % era5_ky_matrix_hhdefintion = era5_k .* sind(era5_direction_bins); % [rad/m] [Eq.3.5.19b Holthuijsen]
    % k = sqrt(era5_kx_matrix_hhdefintion.^2 + era5_ky_matrix_hhdefintion.^2);
    % [era5_E_kx_ky_hhdefinition, Jacobian] = waveNumberSpectrum(era5_d2fd,era5_omega,k);
    % From the above equations, and because we want to orient the grid in line
    % with the SAR image where the vertical axis (ky) is azimuth and the
    % horizontal axis (kx) is range. THIS IS DIFFERENT TO THE CASE OF H&H BUT
    % CHANGING THINGS HERE WILL MEAN NO ISSUES GOING FORWARD. IT IS IMPORTANT
    % TO NOTE THIS CONVENTION IN THE BEGINNING.
    kx_matrix = k .* sind(direction_bins_degrees);
    ky_matrix = k .* cosd(direction_bins_degrees);

% Plot the two versions of the wave number spectrum
    % CHECK THE PLOTTING DIRECTION COMPARED TO POLAR & OCEANOGRAPHY STANDARDS
    % figure('Position', [100, 100, 800, 300]);
    % subplot(1,2,1); contour(era5_kx_matrix_hhdefintion, era5_ky_matrix_hhdefintion, era5_E_kx_ky_hhdefinition,40);
    % xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    % xlabel("kx"); ylabel("ky"); title("OLD Wave Number spectrum E(K) where kx = k * cos(theta) and y = k * sin(theta)");
    % grid on;
    % subplot(1,2,2); contour(era5_kx_matrix, era5_ky_matrix, era5_E_kx_ky,40);
    % xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    % xlabel("kx"); ylabel("ky"); title("Wave Number spectrum E(K) where kx = k * sin(theta) and y = k * cos(theta)");
    % grid on; colorbar;
end

function [wave_spectrum] = waveSpectrum(wave_number_spectrum, omega, k, direction_bins_degrees,sar_transect_size)
%WaveSpectrum convert E(kx,ky) -> E(frequency, theta)

    % Calculate wave phase speed 
    cw = omega./k;
    % Calculate group wave phase speed 
    cg = 0.5 .* cw;
    
    % Calculate the E(w,theta) spectrum
    d2wd = wave_number_spectrum .* omega ./ (cw .* cg);
    % Calculate the E(f, theta) spectrum
    wave_spectrum = d2wd ./ (2 * pi);
    
    final_direction_bins = linspace(min(direction_bins_degrees),max(direction_bins_degrees),sar_transect_size);
    final_frequency_bins = linspace(min(omega(:)./(2*pi)),max(omega(:)./(2*pi)),sar_transect_size);
    [F, D] = meshgrid(final_frequency_bins, final_direction_bins');
    
    figure;
    contour(F,D, wave_spectrum);
    xlabel("Frequency"); ylabel('Direction [degrees]');title("Final 2-D Wave Spectrum, E(f,theta)", "calculated from E(k)"); c = colorbar();c.Label.String = '[m^2 s / degree]"';
end

function [Hs,Tm,wave_direction,total_variance_or_energy] = calculateSpectrumCharacteristics(direction_bins_degrees,frequency_bins,wave_spectrum,sar_transect_size)
%CalculateSpectrumCharacteristics Calculate the integral values of the wave
%spectrum
%   wave_spectrum = E(f, theta)
    
    direction_bins = direction_bins_degrees * pi/180; % Convert to rads
    
% Significant wave height

    % Calculate the 0th order moments (m_n) [Eq.8]
    m0 = trapz(frequency_bins ,trapz(direction_bins,1*wave_spectrum,1),2); %[Eq.8], dimension 2,3 = direc,freq

    % Calculate the significant wave height (Hs) [Eq.4]
    Hs = 4*sqrt(m0);
    
% Mean Wave period
    % Calculate the -1st order moments (m_n)
    m_neg1 = trapz(frequency_bins,trapz(direction_bins,(1./frequency_bins).* wave_spectrum,1),2); %[Eq.8]
    
    % Calculate the mean wave period (Tm) [Eq.5]
    Tm = m_neg1/m0; % [Eq.5]
    
% Mean Wave direction
    % Calculate the angular moments [Eq.9&10]
    a1 = trapz(frequency_bins,trapz(direction_bins,cos(direction_bins_degrees).*wave_spectrum,1),2); %[Eq.9]
    b1 = trapz(frequency_bins,trapz(direction_bins,sin(direction_bins_degrees).*wave_spectrum,1),2); %[Eq.10]
   
    % Calculate the mean Direction of the wave
    wave_direction = atand(b1./a1); %[Eq.6]
    % wave_direction_radians = wave_direction * pi/180;
    
   
% Spectral Denisty
    % Total variance/energy of the waves = spectral density integrated over all
    % frequencies and wave numbers at specified time of interest
    total_variance_or_energy = trapz(trapz(wave_spectrum,2),1);

end

function [interp_E_kx_ky,interp_kx_matrix, interp_ky_matrix] = interpolateWaveNumberSpectrum(sar_transect_size,sar_range_resolution, sar_azimuth_resolution, wave_number_spectrum, kx_matrix, ky_matrix)
%interpolateWaveNumberSpectrum create a 128x128 wave number spectrum
%   The original spectrum is smaller than the SAR image; interpolation is
%   used to increase the size.

    % % Resize the kx and ky matrices.
    % extrap_kx = linspace((min(min(abs(era5_kx_matrix)))),(max(max(abs(era5_kx_matrix)))),sar_sub_transect_size/2);
    % extrap_kx = [-1*fliplr(extrap_kx) extrap_kx];
    interp_kx = (2*pi/(sar_range_resolution*sar_transect_size)) * ([-1*sar_transect_size/2:1:-1 1:1:sar_transect_size/2]);
    
    interp_ky = (2*pi/(sar_azimuth_resolution*sar_transect_size)) * ([-1*sar_transect_size/2:1:-1 1:1:sar_transect_size/2]);
    % extrap_ky = linspace((min(min(abs(era5_ky_matrix)))),(max(max(abs(era5_ky_matrix)))),sar_sub_transect_size/2);
    % extrap_ky = [-1*fliplr(extrap_ky) extrap_ky];
    
    [interp_kx_matrix, interp_ky_matrix] = meshgrid(interp_kx,interp_ky);
    F = scatteredInterpolant(kx_matrix(:) , ky_matrix(:) , wave_number_spectrum(:), 'natural');
    interp_E_kx_ky = F(interp_kx_matrix,interp_ky_matrix);

    % The code below is a sanity check to understand the effects of the interpolation on the spectrum. The peak wave period of the original (old) and interpolated (new) are calculated and compared.
    % % Check the quality of the extrapolation (internal sanity check)
    
    % [test_index1_new, test_index2_new] = find(interp_E_kx_ky==(max(max(interp_E_kx_ky))));
    % test_k_value_new = sqrt(interp_kx_matrix(test_index1_new, test_index2_new).^2 + interp_ky_matrix(test_index1_new, test_index2_new).^2);
    % test_omega_value_new = sqrt(gravity * test_k_value_new);
    % test_period_new = 2*pi/test_omega_value_new
    % [test_index1_old, test_index2_old] = find(era5_wave_number_spectrum==(max(max(era5_wave_number_spectrum))));
    % test_k_value_old = sqrt(era5_kx_matrix(test_index1_old, test_index2_old).^2 + era5_ky_matrix(test_index1_old, test_index2_old).^2);
    % test_omega_value_old = sqrt(gravity .* test_k_value_old);
    % test_period_old = 2*pi/test_omega_value_old
    % % TODO: COMPARE TO THE PEAK WAVE PERIOD
end


