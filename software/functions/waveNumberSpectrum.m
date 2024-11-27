function [era5_E_kx_ky, Jacobian,era5_kx_matrix,era5_ky_matrix] = waveNumberSpectrum(era5_d2fd,era5_omega,era5_k, era5_direction_bins)
%WaveNumberSpectrum convert a wave spectrum E(frequency, direction) to k domain, E(kx,ky)
%   According to Holthuijsen: "Remote-sensing and numerical wave models can 
%   estimate the full two-dimensional spectrum, usually the wave-number 
%   spectrum E (k x , k y )" [pg.52 Holthuijsen] and Ryan 4022 pg 20 
%   onwards:
    
    %% Wave Number Spectrum E(kx,ky)
    % First convert to omega from frequency
    era5_d2wd = era5_d2fd ./ (2 * pi); % conver E(f, theta) to E(w, theta)

    % Wave phase speed (c_w) for deep water
    % era5_cw_wave_speed = gravity ./ (era5_omega); % c_0, [Eq.5.4.24 Holthuijsen]
    era5_cw_wave_speed = era5_omega ./era5_k; % c, [Eq.3.5.36 Holthuijsen] same as above
    
    % Group wave speed (c_g)
    era5_n_dispersion_scaling_factor = 0.5; % 0.5 is used for deep water, [Eq.5.4.32 Holthuijsen]
    era5_cg_group_wave_speed = era5_n_dispersion_scaling_factor * era5_cw_wave_speed; % [Eq.5.4.32 Holthuijsen]
    
    % Calculate the wave number spectrum E(kx,ky)
    Jacobian = ((era5_cw_wave_speed .* era5_cg_group_wave_speed) ./ era5_omega ) ; % [Eq.3.5.36 Holthuijsen]
    era5_E_kx_ky = Jacobian .* era5_d2wd; % [Eq.3.5.36 Holthuijsen]
    
    %% kx and ky variables
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
    era5_kx_matrix = era5_k .* sind(era5_direction_bins);
    era5_ky_matrix = era5_k .* cosd(era5_direction_bins);

    %% Plot the two versions of the wave number spectrum
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
