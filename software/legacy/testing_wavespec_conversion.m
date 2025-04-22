TEST_ORIGINAL = 0;
TEST_AZIMUTH_ANGLE = 0;
TEST_INTERPOLATION = 1;
TEST_INTERPOLATION_AND_AZIMUTH = 0;
TEST_GIACOMO_INTERPOLATION = 0;

% Import Library
waveFunctions = waveLibrary;
handh = hAndH1991Library;
transectFunctions = transectLibrary;
plotFunctions = plotLibrary;

%% Import Sentinel data
filepath = "/Users/tris/Documents/MSc/Data/Cape Point/subset_0_of_S1A_IW_GRDH_1SDV_20241014T173427_T173452_tnr_Cal_msk.nc";

% Create a sentinel 1A satellite object to store all important information.
sentinel = Sentinel1A(filepath);

% Import the data.
sar_latGrid = sentinel.getLatitudeGrid;
sar_lonGrid = sentinel.getLongitudeGrid;
sar_data = sentinel.getCalibratedSARImage("VV");

sar_start_datetime = sentinel.AcquisitionStartDatetime; 
sar_stop_datetime = sentinel.AcquisitionStopDatetime; 

useOriginalS1AHeading = false;
if useOriginalS1AHeading
    sar_azimuth_to_north_angle = sentinel.PlatformHeading;
else 
    sar_azimuth_to_north_angle = sentinel.azimuthToNorthAngleConversion;
end

sar_dy_azimuth = sentinel.AzimuthResolution;
sar_dx_range = sentinel.RangeResolution;

%% Define OO transect
wave_buoy_lat = -34.204;
wave_buoy_lon = 18.28666944;

sar_transect_size = 512; % This is the standard size used.
sar_sub_transect_size = 128;

latitude_of_interest = wave_buoy_lat;
longitude_of_interest = wave_buoy_lon;

[sar_transect_data,sar_transect_lonGrid,sar_transect_latGrid,sar_transect_center_longitude, sar_transect_lon_indices,sar_transect_center_latitude,sar_transect_lat_indices] = transectFunctions.transectFromSARImage(latitude_of_interest,longitude_of_interest,sar_latGrid,sar_lonGrid,sar_data, sar_transect_size);

figure; plotFunctions.sarTransectOnSARImageWithBuoyLocation(sar_data(2000:4000,3000:6000), sar_latGrid(2000:4000, 3000:6000), sar_lonGrid(2000:4000, 3000:6000), sar_transect_data, sar_transect_latGrid, sar_transect_lonGrid,latitude_of_interest,longitude_of_interest, "SAR transect plotted on original SAR image with buoy location");

%% Import the ERA5 model data

% Get Parameters of the complete
filepath = "/Users/tris/Documents/MSc/Data/Cape Point/CapePoint_ERA5-2Dws_20221014.nc";

era5 = ERA5(filepath);

% Get attributes
era5_longitudes = era5.getLongitude;
era5_latitudes = era5.getLatitude;
era5_time = era5.getTime;
% era5_d2fd_all_5dims = era5.getAllWaveSpectrumD2FD;

era5_freq_bins = era5.FrequencyBins;

era5_direction_bins_original_degrees = era5.DirectionBins;

gravity = 9.81;
era5_omega = 2 * pi * era5_freq_bins; % Angular velocity, omega = 2pi * frequency; also calculated as era5_omega = sqrt(gravity .* era5_k)
era5_k = (era5_omega).^2 ./ gravity; % from deep water approximation of omega = sqrt(gk), [Eq.5.4.18 Holthuijsen]

%% Slice and Adjust the ERA5 data
% Slice the spectrum
[era5_d2fd, era5_lat, era5_lon, era5_datetime] = era5.getSlicedWaveSpectrumD2Fd(sar_transect_center_latitude,sar_transect_center_longitude,sar_start_datetime);

era5_direction_bins_degrees = era5_direction_bins_original_degrees;

figure('Position', [100, 100, 800, 300]);
subplot(1,2,1); plotFunctions.waveSpectrum(1,era5_d2fd, era5_freq_bins, era5_direction_bins_original_degrees, 0.3, "2-D Wave Spectrum, E(f,theta)");
subplot(1,2,2);
    % TO DO: FIX THIS METHOD
    plotFunctions.waveSpectrumAdjustedToSAR(era5_d2fd, era5_freq_bins, era5_direction_bins_degrees, 0.3, "Adjusted to SAR scene geometry 2-D Wave Spectrum, E(f,theta)"); % ADD AZ AND RANGE AXIS LABELS

%% Calculate wave number spectrum
% Calculate wave number spectrum E(f,direction) --> E(k)
[era5_wave_number_spectrum, Jacobian, era5_kx_matrix,era5_ky_matrix] = waveFunctions.waveNumberSpectrum(era5_d2fd,era5_omega,era5_k, era5_direction_bins_degrees);
[era5_wave_number_spectrum_t, Jacobian_t, era5_kx_matrix_t,era5_ky_matrix_t] = waveFunctions.waveNumberSpectrum(era5_d2fd,era5_omega,era5_k, era5_direction_bins_original_degrees);

subplot(1,3,1);
    plotFunctions.waveSpectrum(0,era5_d2fd, era5_freq_bins, era5_direction_bins_original_degrees, 0.3, "2-D Wave Spectrum, E(f,theta)");
subplot(1,3,2);
    plotFunctions.waveNumberSpectrum(era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix, "2-D Wave number Spectrum, E(k_x,k_y)");

%% TEST: No angle adjustment or Interpolation
if TEST_ORIGINAL
% Convert back to Wave Spectrum [ E(k) --> E(f,theta)]

% Specify which of the inversion outputs to use.
final_wave_number_spectrum = era5_wave_number_spectrum;

% Reverse the wave number spectrum conversion
% NB: if Step 0 is used then k and omega need to be recalculated.
final_omega = era5_omega;
final_k = era5_k;
final_direction_bins = era5_direction_bins_original_degrees;
final_frequency_bins = era5_freq_bins;


[final_d2fd,inverse_Jacobian] = waveFunctions.waveSpectrum(final_wave_number_spectrum, final_omega, final_k, final_direction_bins,final_frequency_bins);

% Compare the characteristics
[final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(final_direction_bins,final_frequency_bins,final_d2fd,sar_sub_transect_size);
[original_Hs_derived,original_Tm_derived,original_direction_derived,original_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(era5_direction_bins_original_degrees,era5_freq_bins,era5_d2fd,sar_sub_transect_size);

disp("----------------------------------------------")
disp("Testing basic method (no azimuth angle adjustment or interpolation")
disp(" ")
disp("ERA5 original Hs: " + original_Hs_derived)
disp("(After) Hs: " + final_Hs_derived)
disp("ERA5 original Tm: " + original_Tm_derived)
disp("(After) Tm: " + final_Tm_derived)
disp("ERA5 original direction: " + original_direction_derived)
disp("(After) direction: " + final_direction_derived)

end

%% Testing Wave spectrum rotation (azimuth adjustment)
if TEST_AZIMUTH_ANGLE
% We need to rotate the wave spectrum so that it is aligned with the angle
% at which the SAR data has been taken.
era5_direction_bins_degrees = era5.getDirectionsInSARGeometry(sar_azimuth_to_north_angle);

figure('Position', [100, 100, 800, 300]);
subplot(1,2,1); plotFunctions.waveSpectrum(1,era5_d2fd, era5_freq_bins, era5_direction_bins_original_degrees, 0.3, "2-D Wave Spectrum, E(f,theta)");
subplot(1,2,2);
    % TO DO: FIX THIS METHOD
    plotFunctions.waveSpectrumAdjustedToSAR(era5_d2fd, era5_freq_bins, era5_direction_bins_degrees, 0.3, "Adjusted to SAR scene geometry 2-D Wave Spectrum, E(f,theta)"); % ADD AZ AND RANGE AXIS LABELS

% Calculate wave number spectrum E(f,direction) --> E(k)
[era5_wave_number_spectrum, Jacobian, era5_kx_matrix,era5_ky_matrix] = waveFunctions.waveNumberSpectrum(era5_d2fd,era5_omega,era5_k, era5_direction_bins_degrees);
[era5_wave_number_spectrum_t, Jacobian_t, era5_kx_matrix_t,era5_ky_matrix_t] = waveFunctions.waveNumberSpectrum(era5_d2fd,era5_omega,era5_k, era5_direction_bins_original_degrees);

subplot(1,3,1);
    plotFunctions.waveSpectrum(0,era5_d2fd, era5_freq_bins, era5_direction_bins_original_degrees, 0.3, "2-D Wave Spectrum, E(f,theta)");
subplot(1,3,2);
    plotFunctions.waveNumberSpectrum(era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix, "2-D Wave number Spectrum, E(k_x,k_y)");

% Convert back to Wave Spectrum [ E(k) --> E(f,theta)]
final_wave_number_spectrum = era5_wave_number_spectrum;

final_omega = era5_omega;
final_k = era5_k;
final_direction_bins = era5_direction_bins_degrees;
final_frequency_bins = era5_freq_bins;
% final_direction_bins = linspace(min(era5_direction_bins_original_degrees),max(era5_direction_bins_original_degrees),sar_sub_transect_size);
% final_frequency_bins = linspace(min(final_omega(:)./(2*pi)),max(final_omega(:)./(2*pi)),sar_sub_transect_size);
% final_frequency_bins = linspace(min(era5_freq_bins),max(era5_freq_bins),sar_sub_transect_size);
    
    % final_direction_bins = linspace(min(direction_bins_degrees),max(direction_bins_degrees),sar_transect_size);
    % final_frequency_bins = linspace(min(omega(:)./(2*pi)),max(omega(:)./(2*pi)),sar_transect_size);
    % [F, D] = meshgrid(final_frequency_bins, final_direction_bins');
    

[final_d2fd,inverse_Jacobian] = waveFunctions.waveSpectrum(final_wave_number_spectrum, final_omega, final_k, final_direction_bins,final_frequency_bins);

% Compare characteristics
[final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(final_direction_bins,final_frequency_bins,final_d2fd,sar_sub_transect_size);
[original_Hs_derived,original_Tm_derived,original_direction_derived,original_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(era5_direction_bins_original_degrees,era5_freq_bins,era5_d2fd,sar_sub_transect_size);

disp("----------------------------------------------")
disp("Testing with azimuth angle adjustment")
disp(" ")
disp("ERA5 original Hs: " + original_Hs_derived)
disp("(After changing angle) Hs: " + final_Hs_derived)
disp("ERA5 original Tm: " + original_Tm_derived)
disp("(After changing angle) Tm: " + final_Tm_derived)
disp("ERA5 original direction: " + original_direction_derived)
disp("(After changing angle) direction: " + final_direction_derived)

corrected_direction_bins = era5.removeDirectionsInSARGeometry(era5_direction_bins_degrees,sar_azimuth_to_north_angle);
[final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(corrected_direction_bins,final_frequency_bins,final_d2fd,sar_sub_transect_size);
disp("*(After removing azimuth angle ) direction: " + final_direction_derived)

end

%% Testing Interpolation

if TEST_INTERPOLATION
% Apply SAR image angle
era5_direction_bins_degrees =  era5_direction_bins_original_degrees;
% % % Check that the number spectrum was calculated correctly: we should be
% % % able to calculate the original energy from the wave number spectrum.
% % isequal(era5_E_kx_ky(5,5)/scaling_factor(5), era5_d2fd(5,5))

[interp_E_kx_ky,interp_kx_matrix,interp_ky_matrix, interp_k_matrix] = waveFunctions.interpolateWaveNumberSpectrum(sar_sub_transect_size,sar_dx_range, sar_dy_azimuth, era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix);

figure('Position', [200, 200, 1200, 400]);
subplot(1,2,1); plotFunctions.waveNumberSpectrum(era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix, "Original Wave Number Spectrum");
subplot(1,2,2); plotFunctions.waveNumberSpectrum(interp_E_kx_ky, interp_kx_matrix, interp_ky_matrix, "Interpolated Wave Number Spectrum")

plotFunctions.interpolationPlots(era5_k,era5_direction_bins_degrees,era5_kx_matrix,era5_ky_matrix,era5_wave_number_spectrum,interp_k_matrix,interp_kx_matrix,interp_ky_matrix,interp_E_kx_ky);

interp_omega = sqrt(gravity .* interp_k_matrix);

% Convert back to Wave Spectrum [ E(k) --> E(f,theta)]
final_wave_number_spectrum = interp_E_kx_ky;
% final_omega = interp_omega(1,sar_sub_transect_size/2:end);
% final_k = interp_k_matrix(1,sar_sub_transect_size/2:end);
final_omega = interp_omega(1,:);
final_k = interp_k_matrix(1,:);

final_direction_bins = linspace(min(era5_direction_bins_degrees),max(era5_direction_bins_degrees),sar_sub_transect_size);
final_frequency_bins = linspace(min(final_omega(:)./(2*pi)),max(final_omega(:)./(2*pi)),sar_sub_transect_size);
[F, D] = meshgrid(final_frequency_bins, final_direction_bins');
    
[final_d2fd,inverse_Jacobian] = waveFunctions.waveSpectrum(final_wave_number_spectrum, final_omega, final_k, final_direction_bins,final_frequency_bins);

% Compare characteristics
[final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(D(:,1),F(1,:),final_d2fd,sar_sub_transect_size);
[original_Hs_derived,original_Tm_derived,original_direction_derived,original_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(era5_direction_bins_original_degrees,era5_freq_bins,era5_d2fd,sar_sub_transect_size);


disp("----------------------------------------------")
disp("Testing with interpolation")
disp(" ")
disp("ERA5 original Hs: " + original_Hs_derived)
disp("(After Interpolation) Hs: " + final_Hs_derived)
disp("ERA5 original Tm: " + original_Tm_derived)
disp("(After Interpolation) Tm: " + final_Tm_derived)
disp("ERA5 original direction: " + original_direction_derived)
disp("(After Interpolation) direction: " + final_direction_derived)

end

%% Testing Interpolation + Azimuth Angle

if TEST_INTERPOLATION_AND_AZIMUTH
% Apply SAR image angle
era5_direction_bins_degrees = era5.getDirectionsInSARGeometry(sar_azimuth_to_north_angle);
% era5_direction_bins_degrees = era5_direction_bins_original_degrees;
% % % Check that the number spectrum was calculated correctly: we should be
% % % able to calculate the original energy from the wave number spectrum.
% % isequal(era5_E_kx_ky(5,5)/scaling_factor(5), era5_d2fd(5,5))

[interp_E_kx_ky,interp_kx_matrix,interp_ky_matrix, interp_k_matrix] = waveFunctions.interpolateWaveNumberSpectrum(sar_sub_transect_size,sar_dx_range, sar_dy_azimuth, era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix);

figure('Position', [200, 200, 1200, 400]);
subplot(1,2,1); plotFunctions.waveNumberSpectrum(era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix, "Original Wave Number Spectrum");
subplot(1,2,2); plotFunctions.waveNumberSpectrum(interp_E_kx_ky, interp_kx_matrix, interp_ky_matrix, "Interpolated Wave Number Spectrum")

plotFunctions.interpolationPlots(era5_k,era5_direction_bins_degrees,era5_kx_matrix,era5_ky_matrix,era5_wave_number_spectrum,interp_k_matrix,interp_kx_matrix,interp_ky_matrix,interp_E_kx_ky);

interp_omega = sqrt(gravity .* interp_k_matrix);

% Convert back to Wave Spectrum [ E(k) --> E(f,theta)]
final_wave_number_spectrum = interp_E_kx_ky;
% final_omega = interp_omega(1,sar_sub_transect_size/2:end);
% final_k = interp_k_matrix(1,sar_sub_transect_size/2:end);
final_omega = interp_omega(1,:);
final_k = interp_k_matrix(1,:);

final_direction_bins = linspace(min(era5_direction_bins_degrees),max(era5_direction_bins_degrees),sar_sub_transect_size);
final_frequency_bins = linspace(min(final_omega(:)./(2*pi)),max(final_omega(:)./(2*pi)),sar_sub_transect_size);
[F, D] = meshgrid(final_frequency_bins, final_direction_bins');
    
[final_d2fd,inverse_Jacobian] = waveFunctions.waveSpectrum(final_wave_number_spectrum, final_omega, final_k, final_direction_bins,final_frequency_bins);

% Compare characteristics
[final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(D(:,1),F(1,:),final_d2fd,sar_sub_transect_size);
[original_Hs_derived,original_Tm_derived,original_direction_derived,original_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(era5_direction_bins_original_degrees,era5_freq_bins,era5_d2fd,sar_sub_transect_size);


disp("----------------------------------------------")
disp("Testing with azimuth angle adjustment & interpolation")
disp(" ")
disp("ERA5 original Hs: " + original_Hs_derived)
disp("(After Interpolation + azimuth angle) Hs: " + final_Hs_derived)
disp("ERA5 original Tm: " + original_Tm_derived)
disp("(After Interpolation + azimuth angle) Tm: " + final_Tm_derived)
disp("ERA5 original direction: " + original_direction_derived)
disp("(After Interpolation + azimuth angle) direction: " + final_direction_derived)

corrected_direction_bins = era5.removeDirectionsInSARGeometry(era5_direction_bins_degrees,sar_azimuth_to_north_angle);
final_direction_bins = linspace(min(corrected_direction_bins),max(corrected_direction_bins),sar_sub_transect_size);
final_frequency_bins = linspace(min(final_omega(:)./(2*pi)),max(final_omega(:)./(2*pi)),sar_sub_transect_size);

[final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(D(:,1),F(1,:),final_d2fd,sar_sub_transect_size);
disp("*(After removing azimuth angle) Hs: " + final_Hs_derived)
disp("*(After removing azimuth angle) Tm: " + final_Tm_derived)
disp("*(After removing azimuth angle ) direction: " + final_direction_derived)

end


%% Testing Giacomo's Interpolation

if TEST_GIACOMO_INTERPOLATION
% Apply SAR image angle
era5_direction_bins_degrees =  era5_direction_bins_original_degrees;
% era5_direction_bins_degrees = era5.getDirectionsInSARGeometry(sar_azimuth_to_north_angle);
% era5_direction_bins_degrees = era5_direction_bins_original_degrees;
% % % Check that the number spectrum was calculated correctly: we should be
% % % able to calculate the original energy from the wave number spectrum.
% % isequal(era5_E_kx_ky(5,5)/scaling_factor(5), era5_d2fd(5,5))

[A_cart1,A_cart2,A_cart3,A_cart4] = GiacomosInterpolation(sar_dx_range,sar_sub_transect_size, era5_freq_bins,era5_direction_bins_degrees , gravity,era5_d2fd,'ASC');
figure; contour(A_cart1);
figure; contour(A_cart2);
figure; contour(A_cart3);
figure; contour(A_cart4);


% [interp_E_kx_ky,interp_kx_matrix,interp_ky_matrix, interp_k_matrix] = waveFunctions.interpolateWaveNumberSpectrum(sar_sub_transect_size,sar_dx_range, sar_dy_azimuth, era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix);
% 
% figure('Position', [200, 200, 1200, 400]);
% subplot(1,2,1); plotFunctions.waveNumberSpectrum(era5_wave_number_spectrum, era5_kx_matrix, era5_ky_matrix, "Original Wave Number Spectrum");
% subplot(1,2,2); plotFunctions.waveNumberSpectrum(interp_E_kx_ky, interp_kx_matrix, interp_ky_matrix, "Interpolated Wave Number Spectrum")
% 
% plotFunctions.interpolationPlots(era5_k,era5_direction_bins_degrees,era5_kx_matrix,era5_ky_matrix,era5_wave_number_spectrum,interp_k_matrix,interp_kx_matrix,interp_ky_matrix,interp_E_kx_ky);
% 
% interp_omega = sqrt(gravity .* interp_k_matrix);
% 
% % Convert back to Wave Spectrum [ E(k) --> E(f,theta)]
% final_wave_number_spectrum = interp_E_kx_ky;
% % final_omega = interp_omega(1,sar_sub_transect_size/2:end);
% % final_k = interp_k_matrix(1,sar_sub_transect_size/2:end);
% final_omega = interp_omega(1,:);
% final_k = interp_k_matrix(1,:);
% 
% final_direction_bins = linspace(min(era5_direction_bins_degrees),max(era5_direction_bins_degrees),sar_sub_transect_size);
% final_frequency_bins = linspace(min(final_omega(:)./(2*pi)),max(final_omega(:)./(2*pi)),sar_sub_transect_size);
% [F, D] = meshgrid(final_frequency_bins, final_direction_bins');
% 
% [final_d2fd,inverse_Jacobian] = waveFunctions.waveSpectrum(final_wave_number_spectrum, final_omega, final_k, final_direction_bins,final_frequency_bins);
% 
% % Compare characteristics
% [final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(D(:,1),F(1,:),final_d2fd,sar_sub_transect_size);
% [original_Hs_derived,original_Tm_derived,original_direction_derived,original_total_variance_or_energy_derived] = waveFunctions.calculateSpectrumCharacteristics(era5_direction_bins_original_degrees,era5_freq_bins,era5_d2fd,sar_sub_transect_size);
% 
% 
% disp("----------------------------------------------")
% disp("Testing with azimuth angle adjustment & interpolation")
% disp(" ")
% disp("ERA5 original Hs: " + original_Hs_derived)
% disp("(After) Hs: " + final_Hs_derived)
% disp("ERA5 original Tm: " + original_Tm_derived)
% disp("(After) Tm: " + final_Tm_derived)
% disp("ERA5 original direction: " + original_direction_derived)
% disp("(After) direction: " + final_direction_derived)

end

disp("----------------------------------------------")

