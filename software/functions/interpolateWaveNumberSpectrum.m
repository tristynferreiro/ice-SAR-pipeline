function [interp_E_kx_ky,interp_kx_matrix, interp_ky_matrix] = interpolateWaveNumberSpectrum(sar_transect_size,sar_range_resolution, sar_azimuth_resolution, era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix)
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
    F = scatteredInterpolant(era5_kx_matrix(:) , era5_ky_matrix(:) , era5_wave_number_spectrum(:), 'natural');
    interp_E_kx_ky = F(interp_kx_matrix,interp_ky_matrix);

    %%
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