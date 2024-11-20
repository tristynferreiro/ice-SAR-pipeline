function [era5_E_kx_ky, Jacobian] = waveNumberSpectrum(era5_d2wd,era5_omega,era5_k)
%WaveNumberSpectrum convert a wave spectrum E(frequency, direction) to k domain, E(kx,ky)
%   According to Holthuijsen: "Remote-sensing and numerical wave models can 
%   estimate the full two-dimensional spectrum, usually the wave-number 
%   spectrum E (k x , k y )" [pg.52 Holthuijsen] and Ryan 4022 pg 20 
%   onwards:

    % Wave phase speed (c_w) for deep water
    % era5_cw_wave_speed = gravity ./ (era5_omega); % c_0, [Eq.5.4.24 Holthuijsen]
    era5_cw_wave_speed = era5_omega ./era5_k; % c, [Eq.3.5.36 Holthuijsen] same as above
    
    % Group wave speed (c_g)
    era5_n_dispersion_scaling_factor = 0.5; % 0.5 is used for deep water, [Eq.5.4.32 Holthuijsen]
    era5_cg_group_wave_speed = era5_n_dispersion_scaling_factor * era5_cw_wave_speed; % [Eq.5.4.32 Holthuijsen]
    
    % Calculate the wave number spectrum E(kx,ky)
    Jacobian = ((era5_cw_wave_speed .* era5_cg_group_wave_speed) ./ era5_omega ) ; % [Eq.3.5.36 Holthuijsen]
    era5_E_kx_ky = Jacobian .* era5_d2wd; % [Eq.3.5.36 Holthuijsen]
end
