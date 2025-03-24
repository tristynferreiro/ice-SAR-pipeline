%% ==================== INITIALISATION ====================
clear; close all; clc;
waveLibrary = waveLibrary;
handh = hAndH1991Library;
transectFunctions = transectLibrary;
plotFunctions = plotLibrary;

gravity = 9.81;
g = 9.81;

sarPlotsON = 0;
plotERA5 = 0;
plotsON = 0;
plotSARSpecGen = 1;
inversionPlots = 1;

nonlinearity_order = 3;
applyStepZero = 0;
inversion_iterations = 50;

%% ==================== SENTINEL DATA ====================
% Import from file
filepath = "/Users/tris/Documents/MSc/Data/Cape Point/subset_0_of_S1A_IW_GRDH_1SDV_20241014T173427_T173452_tnr_Cal_msk.nc";

% Import metadata information for ease of figuring out where everything is
% stored (this is just a dev step)
% ncInfo_SAR = ncinfo(filepath); 

% Create a sentinel 1A satellite object to store all important information.
sentinel = Sentinel1A(filepath);

% Import the data.
sar_latGrid = sentinel.getLatitudeGrid;
sar_lonGrid = sentinel.getLongitudeGrid;
sar_data = sentinel.getCalibratedSARImage("VV");

sar_start_datetime = sentinel.AcquisitionStartDatetime; 
sar_stop_datetime = sentinel.AcquisitionStopDatetime; 

% TO DO: Investigate the affect of using the original vs using the SNAP
% values. Does it make that big of a difference?
useOriginalS1AHeading = false;
if useOriginalS1AHeading
    sar_azimuth_to_north_angle = sentinel.PlatformHeading;
else 
    sar_azimuth_to_north_angle = sentinel.azimuthToNorthAngleConversion;
end

sar_dy_azimuth = sentinel.AzimuthResolution;
sar_dx_range = sentinel.RangeResolution;

% Plot the data
if sarPlotsON
    figure; plotFunctions.sarImage(1, sar_data, "Sentinel 1A SAR Image - 14 October 2024");
end

% Define TRANSECT
wave_buoy_lat = -34.204;
wave_buoy_lon = 18.28666944;

sar_transect_size = 512; % This is the standard size used.
sar_sub_transect_size = 128;

latitude_of_interest = wave_buoy_lat;
longitude_of_interest = wave_buoy_lon;

[sar_transect_data,sar_transect_lonGrid,sar_transect_latGrid,sar_transect_center_longitude, sar_transect_lon_indices,sar_transect_center_latitude,sar_transect_lat_indices] = transectFunctions.transectFromSARImage(latitude_of_interest,longitude_of_interest,sar_latGrid,sar_lonGrid,sar_data, sar_transect_size);

if sarPlotsON
    figure; plotFunctions.sarTransectOnSARImageWithBuoyLocation(sar_data(2000:4000,3000:6000), sar_latGrid(2000:4000, 3000:6000), sar_lonGrid(2000:4000, 3000:6000), sar_transect_data, sar_transect_latGrid, sar_transect_lonGrid,latitude_of_interest,longitude_of_interest, "SAR transect plotted on original SAR image with buoy location");
end

% CHOOSE WHICH TRANSECT TO USE (for testing)
chosenTransect = 1;
if chosenTransect == 1
    latitude_of_interest = wave_buoy_lat;
    longitude_of_interest = wave_buoy_lon;

    clear sar_transect_data; clear sar_transect_lonGrid; clear sar_transect_latGrid; clear sar_transect_center_longitude; clear sar_transect_lon_indices; clear sar_transect_center_latitude; clear sar_transect_lat_indices; 

    [sar_transect_data,sar_transect_lonGrid,sar_transect_latGrid,sar_transect_center_longitude, sar_transect_lon_indices,sar_transect_center_latitude,sar_transect_lat_indices] = transectFunctions.transectFromSARImage(latitude_of_interest,longitude_of_interest,sar_latGrid,sar_lonGrid,sar_data, sar_transect_size);
    
    % Get the incidence angle at the center of the SAR transect.
    sar_transect_incidence_angle_degrees = sentinel.getIncidenceAngle(sar_transect_lat_indices, sar_transect_lon_indices);
    sentinel.TransectIncidenceAngles = table("open ocean", sar_transect_incidence_angle_degrees, 'VariableNames', {'Field', 'Value'});
    
    % Get the slant range at the center of the SAR transect
    % sar_transect_slant_range = sentinel.getSlantRange(sar_transect_lat_end_index(1),sar_transect_lat_start_index(1), sar_transect_lon_start_index(1),sar_transect_lon_end_index(1));
    sar_transect_slant_range = sentinel.getSlantRange(sar_transect_lat_indices, sar_transect_lon_indices);
    sentinel.TransectSlantRanges = table("open ocean", sar_transect_slant_range, 'VariableNames', {'Field', 'Value'});

    if sarPlotsON
        figure; plotFunctions.sarTransectOnSARImageWithBuoyLocation(sar_data(2000:4000,3000:6000), sar_latGrid(2000:4000, 3000:6000), sar_lonGrid(2000:4000, 3000:6000), sar_transect_data, sar_transect_latGrid, sar_transect_lonGrid,latitude_of_interest,longitude_of_interest, "SAR transect plotted on original SAR image with buoy location");
    end
end
%% ==================== ERA5 DATA ====================
% Get Parameters of the complete
filepath = "/Users/tris/Documents/MSc/Data/Cape Point/CapePoint_ERA5-2Dws_20221014.nc";

% Store file info in MATLAB struct
% ncinfo_era5 = ncinfo(filepath);

era5 = ERA5(filepath);

% Get attributes
era5_longitudes = era5.getLongitude;
era5_latitudes = era5.getLatitude;
era5_time = era5.getTime;

era5_freq_bins = era5.FrequencyBins;
era5_direction_bins_original_degrees = era5.DirectionBins;

[era5_d2fd, era5_lat, era5_lon, era5_datetime] = era5.getSlicedWaveSpectrumD2Fd(sar_transect_center_latitude,sar_transect_center_longitude,sar_start_datetime);
era5_d2fd = era5_d2fd;

era5_direction_bins_degrees = era5_direction_bins_original_degrees;

if plotERA5
    % Plot the 2D spectra
    figure('Position', [100, 100, 800, 300]);
    subplot(1,2,1); plotFunctions.waveSpectrumAdjustedToSAR(era5_d2fd, era5_freq_bins, era5_direction_bins_degrees, 0.3, "Adjusted to SAR scene geometry 2-D Wave Spectrum, E(f,theta)"); % ADD AZ AND RANGE AXIS LABELS
    subplot(1,2,2);plotFunctions.waveSpectrum2D(0,era5_d2fd, era5_freq_bins, era5_direction_bins_original_degrees, 0.3, "2-D Wave Spectrum, E(f,theta)");
end

%% ==================== H&H: First Guess Wave Number Spectrum ====================

[Hs_init,Tm_init,direction_init,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(era5_direction_bins_degrees.*(pi/180),era5_freq_bins,era5_d2fd);
disp(['Initial Hs: ', num2str(Hs_init)]);
disp(['Inital Dir: ', num2str(direction_init)]);

% ADJUST DIRECTION BINS
disp('*==== Start angle adjustment ====*');
% We need to rotate the wave spectrum so that it is aligned with the angle
% at which the SAR data has been taken.
adjusted_direction_bins_degrees = handh.getDirectionsInSARGeometry(era5_direction_bins_degrees,sar_azimuth_to_north_angle);

[Hs_adj,Tm_adj,direction_adj,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(adjusted_direction_bins_degrees.*(pi/180),era5_freq_bins,era5_d2fd);
disp(['Angle Adjusted Hs: ', num2str(Hs_adj)]);
disp(['Angle Adjusted Dir: ', num2str(direction_adj)]);

figure('Position', [100, 100, 800, 300]);
subplot(1,2,1); plotFunctions.waveSpectrum2D(1,era5_d2fd,era5_freq_bins, era5_direction_bins_original_degrees, 0.3, "2-D Wave Spectrum, E(f,theta)");
subplot(1,2,2); plotFunctions.waveSpectrum2D(1,era5_d2fd,era5_freq_bins, adjusted_direction_bins_degrees, 0.3, "Adjusted 2-D Wave Spectrum, E(f,theta)");

% INTERPOLATE/REGRID
disp('*==== Start Interpolation/Regridding ====*');
interp_size = 128;

[S_kxky, Kx,Ky] = wavenumspec(plotsON, era5_d2fd, adjusted_direction_bins_degrees, era5_freq_bins, interp_size, sar_dx_range);

%compute Hs in kx ky space
kappa = sqrt(Kx.^2 + Ky.^2);
S_2Dxy=S_kxky./kappa;
% Compute Hs in kx-ky space
Hs_kxky = 4 * sqrt(trapz(Ky(:,1), trapz(Kx(1,:), S_2Dxy, 2)));
disp(['In kx-ky equispaced, Hs: ', num2str(Hs_kxky)]);

% Plot Cartesian spectrum
figure('Position', [0, 0, 800, 300]);
    subplot(1,2,1);
    contour(Kx,Ky, kappa); title("Kappa = Kx^2 + Ky^2");
    subplot(1,2,2);
    contour(Kx, Ky, S_kxky); title("Wave Num Spec on SAR grid", "E(k)")
    xlim([min(Kx(:)) max(Kx(:))]);
    ylim([min(Ky(:)) max(Ky(:))]); colorbar;

%1.3 First Guess Wave Number Spectra Variables
%These final variables are redundant but helpful for now in the process of understanding things.
first_guess_wave_number_spectrum = S_kxky;
first_guess_kx_range = Kx;
first_guess_ky_azimuth = Ky;
first_guess_k = kappa;
gravity = 9.81;
first_guess_omega = sqrt(gravity .* abs(first_guess_k)); % [Engen, 1995]
%% ==================== H&H: Generate First Guess SAR Spectrum ====================
disp('*==== Start SAR Spec Generation ====*');
disp(['Nonlinearity order: ', num2str(nonlinearity_order)]);

[first_guess_sar_spectrum, TS_k, Tv_k, sar_beta, xi_sqr] = handh.generateSARSpectrumFromWaveNumberSpectrum(sentinel, plotSARSpecGen, 1, nonlinearity_order, sar_sub_transect_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k,first_guess_wave_number_spectrum);

%% ==================== H&H: Observed SAR Spectrum ====================

disp('*==== Start SAR Spec Observed ====*');
cutoff = 450; % This is dependent on th
size_of_filter_window = 5; % The number of elements = dim_filt + 1
width_of_gaussian_lobe = 2; % width_of_gaussian_lobe: Controls the spread of the Gaussian. A larger width means the Gaussian will be wider and smoother.
[observed_sar_spectrum] = observedSARSpectrum(sar_transect_data, sar_transect_size, sar_sub_transect_size, sar_dx_range, first_guess_k, size_of_filter_window, width_of_gaussian_lobe, cutoff,first_guess_sar_spectrum,first_guess_kx_range, first_guess_ky_azimuth);

if plotsON
    figure('Position', [0, 0, 800, 300]);
    subplot(1,2,1); plotFunctions.generalSpectrumPlots(0,first_guess_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Simulated SAR Spectrum, P^S_k");
    subplot(1,2,2); plotFunctions().generalSpectrumPlots(0,observed_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth,  ["SAR Spectrum, P^S_{obs}"," with Butterworth cutoff and Gaussian Filter"]);
end

%% ==================== H&H: Step Zero ====================
transect_number = 1; % This will mainly be used when looping through all transects in the ice part of the pipeline
applyStepZero = 1;
if  applyStepZero
    kx = first_guess_kx_range; % This is required for integration
    ky = first_guess_ky_azimuth; % This is required for integration
 

    P_obs = observed_sar_spectrum;
    
    cosmo = sentinel;
    plotsON = 0;

    % START OF ALGO
    first_iteration_check = true;
    
    % [temp1, temp2] = find(P_obs==max(P_obs(:)));
    % k0 = first_guess_k(temp1,temp2);
    % modk = first_guess_k;
    % anello  = 1 - 1./(1+(k0./modk).^7); % MIGHT NEED TO ADJUST POWER (5 in this case) FOR DIFFERENT CASES
    anello = 1;
    
    rotation_angles = -10:1:10;
    wave_number_mag_scales = 0.9:0.1:1.1;
    energy_scales = 0.5:0.5:10;
    % rotation_angles = 0;
    % wave_number_mag_scales = 1;
    % energy_scales =1;

    J_step0 = zeros(length(rotation_angles), length(wave_number_mag_scales), length(energy_scales));
    
    % First interpolate the wave spectrum to a new wave number grid that is
    % scaled to the max possible value as per the wave_number_mag_scales
    % array.
    % max_scale = max(wave_number_mag_scales); 

    [S_tmp_0, kx_tmp,ky_tmp] = wavenumspecStepzero(plotsON, era5_d2fd, adjusted_direction_bins_degrees, era5_freq_bins, interp_size, sar_dx_range);% Convert to Cartesian (kx, ky) coordinates
   
    % max_scale = 2;
    % kx_tmp=linspace(-max_scale*max(kx(:)),max_scale*max(kx(:)),2*128) .* ones(2*128);
    % ky_tmp=linspace(-max_scale*max(ky(:)),max_scale*max(ky(:)),128*2)' .* ones(128*2);
    % 
    % Reset the temporary spectrum
    % S_tmp_0 = interp2(kx,ky,first_guess_wave_number_spectrum,kx_tmp,ky_tmp); % Reset S_tmp to original spectrum
    % S_tmp_0(isnan(S_tmp_0))=0;      
    
    for i_rot = 1:length(rotation_angles)

        % [Eq.77 HH1991]
        % theta = adjusted_direction_bins_degrees + rotation_angles(i_rot);
        % [S_tmp, rot_kx,rot_ky] = wavenumspec(plotsON, era5_d2fd, theta, era5_freq_bins, interp_size, sar_dx_range);

        for i_mag = 1:length(wave_number_mag_scales)
         
            
            for i_energy = 1:length(energy_scales)
                % METHODS THAT DID NOT WORK
                % 1. Using the [Eq.78 HH1991] plus this calculation for
                % S_tmp
                % S_tmp = interp2(kx, ky, first_guess_wave_number_spectrum,...
                %   new_kx, new_ky, 'linear', 0); % Using this results in
                % some weird outputs... good for a making a note in the
                % report to say that it does not work
    
                % 2. Calculating new ky and new ky
                % theta = -1 * rotation_angles(i_rot); % The imrotate performs counterclockwise rotation.
                % new_kx = imrotate(kx,theta,"bilinear","crop") .* wave_number_mag_scales(i_mag);
                % new_ky = imrotate(ky,theta,"bilinear","crop") .* wave_number_mag_scales(i_mag);
                % The problem with this is that it creates a matrix which
                % is not usable for the rest of the calculations. IT also
                % does not do exactly what we want, also better to do this
                % rotation from first principles rather than using built-in
                % functions.
                
                % 3. Rotation using IMROtate. This treats the spectrum like
                % an image and causes weird crops in the spectrum and
                % miscorrelation with the spectrum. 
                % theta = -1 * rotation_angles(i_rot); % The imrotate performs counterclockwise rotation.
                % S_tmp = imrotate(first_guess_wave_number_spectrum,theta,"bilinear","crop");
                
                % Reset the wave number spec
                S_tmp = S_tmp_0;

                % [Eq.78 HH1991]
                B = wave_number_mag_scales(i_mag);
                new_kx = B .* (kx_tmp .* cosd(rotation_angles(i_rot)) - ky_tmp .* sind(rotation_angles(i_rot)));
                new_ky = B .* (kx_tmp .* sind(rotation_angles(i_rot)) + ky_tmp .* cosd(rotation_angles(i_rot)));
                
                % [Eq.78 HH1991]
                A = energy_scales(i_energy);
                S_tmp = energy_scales(i_energy) .* S_tmp;
                
                F = scatteredInterpolant(new_kx(:), new_ky(:), S_tmp(:), 'linear', 'linear');
                S_tmp_n = F(kx, ky);
                S_tmp_n(isnan(S_tmp_n)) = 0;
               
                % Calculate the new SAR spectrum
                plotsON = 0;
                [P_tmp, TS_k, Tv_k, sar_beta, xi_sqr] = handh.generateSARSpectrumFromWaveNumberSpectrum(sentinel, plotsON, transect_number, nonlinearity_order, sar_sub_transect_size, kx, ky, first_guess_omega, first_guess_k, S_tmp_n); % [Eq.65 HH1991]
                P_tmp(isnan(P_tmp)) =0;

                kx_vector = kx(1,:);
                ky_vector = ky(:,1);
                % J_step0(i_rot, i_mag, i_energy) = trapz(kx_vector, trapz(ky_vector, (P_tmp - P_obs).^2 .* anello,1),2) ./ sqrt(trapz(kx_vector, trapz(ky_vector, (P_obs).^2 .* anello,1),2)) ./ sqrt(trapz(kx_vector, trapz(ky_vector, (P_tmp).^2  .* anello,1),2));
                % 
                B = 0.01 * max(first_guess_wave_number_spectrum(:))^2; % Below [Eq.76 HH1991]
                mu = 0.1 * max(P_obs(:))^2; % [Eq.76 HH1991]
                term1_63 = (P_tmp - P_obs).^2 .* P_obs;
                term1_63(isnan(term1_63))=0;
                term2_63 = mu .* ( (S_tmp_n - first_guess_wave_number_spectrum) ./ (B + first_guess_wave_number_spectrum) ).^2;
                term2_63(isnan(term2_63))=0;
                J_step0(i_rot, i_mag, i_energy) = trapz( kx_vector,trapz(ky_vector,term1_63,1),2) + trapz( kx_vector,trapz(ky_vector,term2_63,1),2);% dimensions = Y x X = 1 x 2;
      
                if first_iteration_check
                    Jmax_step0 = J_step0(i_rot, i_mag, i_energy);
                    first_iteration_check = false;
                end
                    
    
                if J_step0(i_rot, i_mag, i_energy) <= Jmax_step0
                    disp("iteration "+i_rot +"," + i_mag+","+i_energy)
                    Jmax_step0    = J_step0(i_rot, i_mag, i_energy);
                    F_new      = S_tmp_n;
                    P_new      = P_tmp;
                    energy_0     = energy_scales(i_energy);
                    mag_0 = wave_number_mag_scales(i_mag);
                    rot_angle_0  = rotation_angles(i_rot);
                    energy_0_i = i_energy;
                    TS_k_0 = TS_k;
                    Tv_k_0 = Tv_k; 
                    sar_beta_0 = sar_beta;
                    xi_sqr_0 = xi_sqr;
                end
            
            end

        end
    end
    
    S_new = F_new;
    energy_scale = energy_0;
    rot_angle = rot_angle_0;
    mag_scale = mag_0; 
    TS_k = TS_k_0;
    Tv_k = Tv_k_0;
    sar_beta = sar_beta_0;
    xi_sqr = xi_sqr_0;

    % Plot the outputs
    figure('Position', [0, 0, 1600, 800]);
        subplot(2,4,1); plotFunctions.generalSpectrumPlots(0,observed_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Observed SAR spectrum");
        subplot(2,4,2); plotFunctions.generalSpectrumPlots(0,first_guess_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "First guess SAR spectrum");
        subplot(2,4,3); plotFunctions.generalSpectrumPlots(0,P_new, kx, ky, "P(k)");
        subplot(2,4,4);
        % ix
        % countJ=J(i_rot,)
            % pcolor(J_step0(:));
            % clear temp; temp(:,:) = J_step0(:,1,:); surf(temp); view(2);
            % title("Minimisation Values"); ylabel("Wave Number Scaling Factor"); xlabel("Rotation angle (degrees)"); 
            % c = colorbar(); c.Label.String = 'Minimisation magnitude';
        subplot(2,4,6); plotFunctions.generalSpectrumPlots(0,first_guess_wave_number_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "First guess wave number spectrum");
        subplot(2,4,7); plotFunctions.generalSpectrumPlots(0,S_new, kx, ky, "E(k)");
    
    % Update the variables to the new values
    first_guess_wave_number_spectrum = S_new;
    first_guess_kx_range = kx;
    first_guess_ky_azimuth = ky;
    % first_guess_omega = sqrt(9.81 * first_guess_k);
end
%% ==================== H&H: Inversion ====================
disp('*==== Start Inversion ====*');
[J_eq_62,P_best_eq62,F_best_eq62,J_eq_63,P_best_eq63,F_best_eq63,J_eq_69,P_best_eq69,F_best_eq69] = inversionHH(inversion_iterations, nonlinearity_order, sentinel, handh, 0, observed_sar_spectrum, first_guess_wave_number_spectrum, transect_number, sar_sub_transect_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k );
%TO DO: ADD MARKER TO SHOW WHICH ITERATION WAS CHOSEN
if inversionPlots
    figure('Position', [0, 0, 1200, 900]);
    subplot(3,3,1); plotFunctions.generalSpectrumPlots(0,F_best_eq62, first_guess_kx_range, first_guess_ky_azimuth, ["Final Wave Spectrum, F_k","From Eq.62 Inversion"]);
    subplot(3,3,2); plotFunctions.generalSpectrumPlots(0,P_best_eq62, first_guess_kx_range, first_guess_ky_azimuth, ["Final SAR Spectrum, P^S_{k}"," From Eq.62 Inversion"]);
    subplot(3,3,3); 
        plot(J_eq_62); xticks((1:1:inversion_iterations)); xlabel("Number of iterations");
    subplot(3,3,4); plotFunctions.generalSpectrumPlots(0,F_best_eq63, first_guess_kx_range, first_guess_ky_azimuth, ["Final Wave Spectrum, F_k"," From Eq.63 Inversion"]);
    subplot(3,3,5); plotFunctions.generalSpectrumPlots(0,P_best_eq63, first_guess_kx_range, first_guess_ky_azimuth, ["Final SAR Spectrum, P^S_{k}"," From Eq.63 Inversion"]);
    subplot(3,3,6);
        plot(J_eq_63); xticks((1:1:inversion_iterations)); xlabel("Number of iterations");
    subplot(3,3,7); plotFunctions.generalSpectrumPlots(0,F_best_eq69, first_guess_kx_range, first_guess_ky_azimuth, ["Final Wave Spectrum, F_k"," From Eq.69 Inversion"]);    
    subplot(3,3,8); plotFunctions.generalSpectrumPlots(0,P_best_eq69, first_guess_kx_range, first_guess_ky_azimuth, ["Final SAR Spectrum, P^S_{k}","From Eq.69 Inversion"]);    
    subplot(3,3,9); 
        plot(J_eq_69); xticks((1:1:inversion_iterations)); xlabel("Number of iterations");
end

chosen_final_spectrum = F_best_eq62;
chosen_final_SAR_spectrum = P_best_eq62;

%% ==================== RESULTS ====================
disp('*==== OUTPUTS ====*');
figure('Position', [100, 100, 1200, 300]);
    subplot(1,4,1); plotFunctions.generalSpectrumPlots(0,first_guess_wave_number_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "First guess wave number spectrum");
    subplot(1,4,2); plotFunctions.generalSpectrumPlots(0,chosen_final_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Final Wave Spectrum, F_k: From Eq.62 Inversion");
    subplot(1,4,3); plotFunctions.generalSpectrumPlots(0,chosen_final_SAR_spectrum, first_guess_kx_range, first_guess_ky_azimuth,"Final SAR Spectrum: From Eq.62 Inversion");
    subplot(1,4,4); plotFunctions.generalSpectrumPlots(0,observed_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Observed SAR spectrum");

% Convert E(k) - E(k,theta)
kappa = sqrt(Kx.^2 + Ky.^2);
% Compute Hs in kx-ky space
Hs_kxky = 4 * sqrt(trapz(Ky(:,1), trapz(Kx(1,:), chosen_final_spectrum./kappa, 2)));
disp(['In SAR kx-ky equispaced, Hs: ', num2str(Hs_kxky)]);

% Recovered variables
% dimension of this grid can be selected as preferred (in principle can
% have any dimension)
n = interp_size;
k_r=logspace(log10(min(abs(kappa(:)))),log10(max(Kx(:))),n);
theta_r=linspace(0,2*pi,n);

for i=1:length(k_r)
    for j=1:length(theta_r)
        kx_r(i,j)=k_r(i)*sin(theta_r(j));
        ky_r(i,j)=k_r(i)*cos(theta_r(j));
    end
end

% Interpolate onto new grid
F = scatteredInterpolant(Kx(:), Ky(:), chosen_final_spectrum(:), 'natural', 'none');
S_r = F(kx_r, ky_r);
S_r(isnan(S_r))=0;

kappa=sqrt(kx_r.^2+ky_r.^2);
S_rxy=S_r./kappa;
% Compute Hs in kx-ky space
% Hs_kxky = 4 * sqrt(trapz(kx_r(:,1), trapz(ky_r(1,:), S_rxy, 2)));
% disp(['In kx-ky polar, Hs: ', num2str(Hs_kxky)]);

Hs_m0 = 4 * sqrt(trapz(k_r, trapz(theta_r, S_r, 2), 1));
disp(['(after regrid) In k theta space, Hs: ', num2str(Hs_m0)]);

% INVERSE_TRANSFORM_SPECTRUM Converts a wave spectrum from wavenumber (k) to frequency (f)
f_r = sqrt(k_r * g) / (2 * pi);
dk_df = (8 * pi^2 / g) * f_r';
S_f_r = S_r .* dk_df;

% Plot recovered 2D spectrum
figure('Position', [0, 0, 1200, 300]);
    subplot(1,3,1);
    contour(Kx, Ky, chosen_final_spectrum); title("Final chosen Spectrum", "on SAR kx,ky")
    xlim([min(Kx(:)) max(Kx(:))]);
    ylim([min(Ky(:)) max(Ky(:))]);
    subplot(1,3,2);
    contour(kx_r, ky_r, S_r); title("Final spectrum regridded to new kx_r, ky_r")
    xlim([min(Kx(:)) max(Kx(:))]);
    ylim([min(Ky(:)) max(Ky(:))]);
    subplot(1,3,3);
    plotLibrary().waveSpectrum2D(1, abs(S_f_r'), f_r, theta_r * 180 / pi, 1, "Recovered 2D wave spectrum");

% Compute recovered wave height
S_f_r(isnan(S_r)) = 0;
Hs_m0 = 4 * sqrt(trapz(f_r, trapz(theta_r, S_f_r, 2), 1));
disp(['Recovered Hs: ', num2str(Hs_m0)]);

% Plot recovered vs original 1D spectrum
S_1Dr = trapz(theta_r, S_f_r, 2);

S_2D = era5_d2fd;
theta = deg2rad(adjusted_direction_bins_degrees);  % Convert degrees to radians
f = era5_freq_bins;
S_1D = trapz(theta, S_2D', 2);
disp(['Inital ERA5 Hs: ', num2str(Hs_init)]);
Hs_kxky = 4 * sqrt(trapz(Ky(:,1), trapz(Kx(1,:), S_2Dxy, 2)));
disp(['Initial SAR kx-ky, Hs: ', num2str(Hs_kxky)]);

% figure;
%     plot(f_r, S_1Dr, 'k', 'DisplayName', 'Recovered');
%     hold on;
%     plot(f, S_1D, 'r', 'DisplayName', 'Original');
%     legend;
%     title('Recovered vs Original 1D Spectrum');
%     set(gca,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
figure;
    plot(f_r, S_1Dr, 'k', 'DisplayName', 'Recovered');
    hold on;
    plot(f, S_1D, 'r', 'DisplayName', 'Original');
    legend;
    title('Recovered vs Original 1D Spectrum');

%% ==================== COMPARISON ====================
disp('*==== Comparison ====*');
[final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveLibrary.calculateSpectrumCharacteristics(theta_r,f_r,S_f_r');
disp(['Final Hs: ',  num2str(final_Hs_derived)]);
disp(['Final Tm: ',  num2str(final_Tm_derived)]);
disp(['Final Dir: ',  num2str(final_direction_derived)]);
disp(['Final Dir - SAR: ',  num2str(final_direction_derived-abs(sar_azimuth_to_north_angle))]);

wavelength_Tm_final = ((final_Tm_derived)^2 * 9.81)/(2*pi);
disp(['Final wavelength from Tm: ', num2str(wavelength_Tm_final)])

[~,max_index] = max(S_1Dr);
final_Tp = 1/f_r(max_index);
disp(['Final Tp: ', num2str(final_Tp)]);
wavelength_Tp_final = ((final_Tp)^2 * 9.81)/(2*pi);
disp(['Final wavelength from Tp: ', num2str(wavelength_Tp_final)]);

disp('--------------------------------------');
disp("ERA5 date"); disp(era5_datetime)
disp(['ERA5 lon: ',  num2str(era5_lon)]);
disp(['ERA5 lat: ',  num2str(era5_lat)]);
[era5_Hs_derived,era5_Tm_derived,era5_direction_derived,~] = waveLibrary.calculateSpectrumCharacteristics(era5_direction_bins_original_degrees.*(pi/180),era5_freq_bins,era5_d2fd);

disp(['Inital ERA5 Hs: ', num2str(era5_Hs_derived)]);
disp(['Inital ERA5 Tm: ', num2str(era5_Tm_derived)]);
disp(['Inital ERA5 Dir: ', num2str(era5_direction_derived)]);


disp('--------------------------------------');
filepath_single = "/Users/tris/Documents/MSc/Data/Cape Point/CapePoint_ERA5-single-2Dws_20241014.nc";
[era5_datetime_single, era5_lon_single, era5_lat_single, era5_direction,era5_Tm,era5_Hs, era5_mdts,era5_mpts,era5_shts, era5_mdww, era5_mpww, era5_shww] = era5.getSingleValues(filepath_single, sar_transect_center_latitude,sar_transect_center_longitude,sar_start_datetime);

disp("Singles ERA5 date"); disp(era5_datetime_single)
disp(['Singles ERA5 Hs: ', num2str(era5_Hs)]);
disp(['Singles ERA5 Tm: ', num2str(era5_Tm)]);
disp(['Singles ERA5 mdts: ', num2str(era5_mdts)]);
disp(['Singles ERA5 mpts: ', num2str(era5_mpts)]);
disp(['Singles ERA5 shts: ', num2str(era5_shts)]);
disp(['Singles ERA5 mdww: ', num2str(era5_mdww)]);
disp(['Singles ERA5 mpww: ', num2str(era5_mpww)]);
disp(['Singles ERA5 shww: ', num2str(era5_shww)]);

wavelength_mww = ((era5_mpww)^2 * 9.81)/(2*pi);

disp(['Singles ERA5 wavelength from mpww: ', num2str(wavelength_mww)])
%% ==================== FUNCTIONS ====================
function [S_kxky, Kx,Ky] = wavenumspec(plotsON, wavespec, dirBins, f_bins, interp_size, sar_dx_range)
    f = f_bins;
    n = interp_size;
    dx = sar_dx_range;
    g = 9.81;
    S_2D = wavespec;
    numDirBins = length(dirBins);
    
    % Compute wavenumber (k) using deep-water dispersion relation
    k = (2 * pi * f).^2 / g;
    theta = deg2rad(dirBins);  % Convert degrees to radians
    
    % TRANSFORM_SPECTRUM Converts a wave spectrum from frequency (f) to wavenumber (k)
    % Compute the Jacobian of the transformation df/dk:
    df_dk = 1 ./ (4 * pi) .* sqrt(g ./ k);
    % S(k, theta) = S(f, theta) * |df/dk|
    S_k_theta = S_2D .* df_dk;
    
    Hs_m0 = 4 * sqrt(trapz(k, trapz(theta, S_k_theta', 2), 1));
    [Hs_derived,Tm_derived,direction_derived,~] = ...
            waveLibrary().calculateSpectrumCharacteristics(theta,k,S_k_theta);
    
    disp(['In k theta space Hs: ', num2str(Hs_m0)]);
    % Plot transformed spectrum
    if plotsON
        figure;
        plotLibrary().waveSpectrum2D(1, abs(S_k_theta), k, theta * 180 / pi, 0.3, "Wave Spectrum in (k, theta)");
    end

    % Convert to Cartesian (kx, ky) coordinates
    [kk,tt] = meshgrid(k,theta);
    y = kk(:) .* cos(tt(:));
    x =  kk(:) .* sin(tt(:));
    z = S_k_theta(:);
    % Define Cartesian grid
    dk = 2 * pi / (n * dx);
    kx = linspace(-pi/dx,pi/dx,n);
    kx = kx(kx~=0);
    ky = kx;
    [Kx, Ky] = meshgrid(kx,ky);
    F = scatteredInterpolant(x, y, z, 'linear', 'none');
    S_kxky=F(Kx, Ky);
    
    % if plotsON
    %     figure;
    %     contour(Kx, Ky, S_kxky);
    %     xlim([min(Kx(:)) max(Kx(:))]);
    %     ylim([min(Kx(:)) max(Kx(:))]);
    %     title('Wave Number Spectrum on SAR grid');
    %     colorbar;
    % end
    
    
end

function [S_tmp0, Kx_tmp,Ky_tmp] = wavenumspecStepzero(plotsON, wavespectrum, dirBins, f_bins, interp_size, sar_dx_range)
    f = f_bins;
    n = interp_size;
    dx = sar_dx_range;
    g = 9.81;
    S_2D = wavespectrum;
    numDirBins = length(dirBins);
    
    % Compute wavenumber (k) using deep-water dispersion relation
    k = (2 * pi * f).^2 / g;
    theta = deg2rad(dirBins);  % Convert degrees to radians
    
    % TRANSFORM_SPECTRUM Converts a wave spectrum from frequency (f) to wavenumber (k)
    % Compute the Jacobian of the transformation df/dk:
    df_dk = 1 ./ (4 * pi) .* sqrt(g ./ k);
    % S(k, theta) = S(f, theta) * |df/dk|
    S_k_theta = S_2D .* df_dk;
    
    Hs_m0 = 4 * sqrt(trapz(k, trapz(theta, S_k_theta', 2), 1));
    [Hs_derived,Tm_derived,direction_derived,~] = ...
            waveLibrary().calculateSpectrumCharacteristics(theta,k,S_k_theta);
    
    disp(['In k theta space Hs: ', num2str(Hs_m0)]);
    % Plot transformed spectrum
    if plotsON
        figure;
        plotLibrary().waveSpectrum2D(1, abs(S_k_theta), k, theta * 180 / pi, 0.3, "Wave Spectrum in (k, theta)");
    end

    % Convert to Cartesian (kx, ky) coordinates
    [kk,tt] = meshgrid(k,theta);
    y = kk(:) .* cos(tt(:));
    x =  kk(:) .* sin(tt(:));
    z = S_k_theta(:);
    % Define Cartesian grid
    dk = 2 * pi / (n * dx);
    kx = 10*linspace(-pi/dx,pi/dx,n);
    kx = kx(kx~=0);
    ky = kx;
    [Kx_tmp, Ky_tmp] = meshgrid(kx,ky);
    F = scatteredInterpolant(x, y, z, 'linear', 'none');
    S_tmp0=F(Kx_tmp, Ky_tmp);
    
    % if plotsON
    %     figure;
    %     contour(Kx, Ky, S_kxky);
    %     xlim([min(Kx(:)) max(Kx(:))]);
    %     ylim([min(Kx(:)) max(Kx(:))]);
    %     title('Wave Number Spectrum on SAR grid');
    %     colorbar;
    % end
    
    
end