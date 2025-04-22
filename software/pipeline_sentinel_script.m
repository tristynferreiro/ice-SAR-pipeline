%% ==================== INITIALISATION ====================
clear; close all; clc;
waveLibrary = waveLibrary;
handh = hAndH1991Library;
transectLibrary = transectLibrary;
plotLibrary = plotLibrary;

g = 9.81;

showScenePlots = 0;
plotsON = 0;
plotSARSpecGen = 0;
inversionPlots = 0;

nonlinearity_order = 3;
applyStepZero = 1;
inversion_iterations = 50;

wave_buoy_lat = -34.204;
wave_buoy_lon = 18.28666944;
latitude_of_interest = wave_buoy_lat;
longitude_of_interest = wave_buoy_lon;
sar_imagette_size = 512; % This is the standard size used.
imagette_size = 128;

%% ==================== SENTINEL DATA ====================
% Import from file
filepath = "/Users/tris/Documents/MSc/data/Cape Point Sentinel/14 Oct 2024/subset_0_of_S1A_IW_GRDH_1SDV_20241014T173427_T173452_tnr_Cal_msk.nc";

% Import metadata information for ease of figuring out where everything is
% stored (this is just a dev step)
% ncInfo_SAR = ncinfo(filepath); 

% Create a sentinel 1A satellite object to store all important information.
sentinel = Sentinel1A(filepath);

% Import the data.
sar_latGrid = sentinel.getLatitudeGrid;
sar_lonGrid = sentinel.getLongitudeGrid;
sar_dy_azimuth = sentinel.AzimuthResolution;
sar_dx_range = sentinel.RangeResolution;
sar_start_datetime = sentinel.AcquisitionStartDatetime; 

sar_data = sentinel.getCalibratedSARImage("VV");
sar_data = (sar_data-mean(sar_data(:)))./mean(sar_data(:)); % normalise [Monaldo] = relative intensity

% TO DO: Investigate the affect of using the original vs using the SNAP
% values. Does it make that big of a difference?
useOriginalS1AHeading = false;
if useOriginalS1AHeading
    sar_az_to_trueN_angle = sentinel.PlatformHeading;
else 
    sar_az_to_trueN_angle = sentinel.azimuthToNorthAngleConversion;
end

% Define TRANSECT
% * CHOOSE WHICH TRANSECT TO USE (for testing)
chosenTransect = 1;
if chosenTransect == 1
    latitude_of_interest = wave_buoy_lat;
    longitude_of_interest = wave_buoy_lon;

    clear sar_transect_data; clear sar_transect_lonGrid; clear sar_transect_latGrid; clear sar_transect_center_longitude; clear sar_transect_lon_indices; clear sar_transect_center_latitude; clear sar_transect_lat_indices; 

    [sar_transect_data,sar_transect_lonGrid,sar_transect_latGrid,sar_transect_center_longitude, sar_transect_lon_indices,sar_transect_center_latitude,sar_transect_lat_indices] = transectLibrary.transectFromSARImage(latitude_of_interest,longitude_of_interest,sar_latGrid,sar_lonGrid,sar_data, sar_imagette_size);
    %sar_transect_data = (sar_transect_data-mean(sar_transect_data(:)))./mean(sar_transect_data(:)); % normalise [Monaldo] = relative intensity

    % Get the incidence angle at the center of the SAR transect.
    sar_transect_incidence_angle_degrees = sentinel.getIncidenceAngle(sar_transect_lat_indices, sar_transect_lon_indices);
    sentinel.TransectIncidenceAngles = table("open ocean", sar_transect_incidence_angle_degrees, 'VariableNames', {'Field', 'Value'});
    
    % Get the slant range at the center of the SAR transect
    % sar_transect_slant_range = sentinel.getSlantRange(sar_transect_lat_end_index(1),sar_transect_lat_start_index(1), sar_transect_lon_start_index(1),sar_transect_lon_end_index(1));
    sar_transect_slant_range = sentinel.getSlantRange(sar_transect_lat_indices, sar_transect_lon_indices);
    sentinel.TransectSlantRanges = table("open ocean", sar_transect_slant_range, 'VariableNames', {'Field', 'Value'});
end
%% ==================== ERA5 DATA ====================
% Get Parameters of the complete
filepath = "/Users/tris/Documents/MSc/data/Cape Point ERA5/CapePoint_ERA5_2Dws_20241014.nc";

% Store file info in MATLAB struct
% ncinfo_era5 = ncinfo(filepath);

era5 = ERA5(filepath);

% Get attributes
era5_f = era5.FrequencyBins;
era5_dir_degrees = era5.DirectionBins;

[era5_2D, era5_lat, era5_lon, era5_datetime] = era5.getSlicedWaveSpectrumD2Fd(sar_transect_center_latitude,sar_transect_center_longitude,sar_start_datetime);

[Hs,Tm,Dirm,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(era5_dir_degrees.*(pi/180),era5_f,era5_2D);
disp(['Initial Hs: ', num2str(Hs)]);
disp(['Initial Tm: ', num2str(Tm)]);
disp(['Inital Dir: ', num2str(Dirm)]);

% Calculate 1D Spectrum
S_1D = trapz(era5_dir_degrees.*(pi/180), era5_2D, 1);

% Get singles values for filtering
if month(era5_datetime) < 10 %Check if it is a single digit
    date_str = num2str(year(era5_datetime))+"0"+num2str(month(era5_datetime))+""+num2str(day(era5_datetime));

else
    date_str = num2str(year(era5_datetime))+""+num2str(month(era5_datetime))+""+num2str(day(era5_datetime));

end
filepath_single = "/Users/tris/Documents/MSc/data/Cape Point ERA5/CapePoint_ERA5_single_2Dws_"+date_str+".nc";

[era5_datetime_single, era5_lon_single, era5_lat_single, era5_direction,era5_Tm,era5_Hs, era5_mdts,era5_mpts,era5_shts, era5_mdww, era5_mpww, era5_shww] = era5.getSingleValues(filepath_single, era5_lat,era5_lon,sar_start_datetime);

wavelength_mpts = ((era5_mpts)^2 * 9.81)/(2*pi);

%% ============================ Plot Scene ============================================
if showScenePlots
    % Plot SAR Transect and matching ERA5 2D spectrum
    figure('Position', [0, 0, 800, 300]);
        ax(1) = subplot(1,2,1);
        plotLibrary.sarTransect(sar_transect_data, sar_transect_latGrid, sar_transect_lonGrid, ["SAR Transect",string(sar_start_datetime)])
        colormap(ax(1),"gray");
        ax(2) = subplot(1,2,2);
        plotLibrary.waveSpectrum2D(1,era5_2D, era5_f, era5_dir_degrees, 0.3, ["ERA5 2D Wave Spectrum, E(f,theta)",string(era5_datetime)]);
        colormap(ax(2),"parula");

    % Plot SAR scene
    figure; plotLibrary.sarTransectOnSARImageWithBuoyLocationWithERA5(sar_data, sar_latGrid, sar_lonGrid, sar_transect_data, sar_transect_latGrid, sar_transect_lonGrid,latitude_of_interest,longitude_of_interest,era5_lat,era5_lon, ["SAR transect plotted on original SAR image with buoy location",string(sar_start_datetime)]);
end 
clear sar_data sar_latGrid sar_lonGrid;
%% ==================== H&H: First Guess Wave Number Spectrum ====================
% ADJUST DIRECTION BINS
disp('*==== Start angle adjustment ====*');
% We need to rotate the wave spectrum so that it is aligned with the angle
% at which the SAR data has been taken.
adjusted_dir_degrees = handh.getDirectionsInSARGeometry(era5_dir_degrees,sar_az_to_trueN_angle);

[Hs_adj,Tm_adj,Dirm_adj,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(adjusted_dir_degrees.*(pi/180),era5_f,era5_2D);
disp(['Angle Adjusted Hs: ', num2str(Hs_adj)]);
disp(['Angle Adjusted Tm: ', num2str(Tm_adj)]);
disp(['Angle Adjusted Dir: ', num2str(Dirm_adj)]);

if plotsON
    figure('Position', [100, 100, 800, 300]);
        subplot(1,2,1); plotLibrary.waveSpectrum2D(1,era5_2D,era5_f, era5_dir_degrees, 0.3, "ERA5 2D Wave Spectrum, E(f,theta)");
        subplot(1,2,2); plotLibrary.waveSpectrum2D(1,era5_2D,era5_f, adjusted_dir_degrees, 0.3, "Rotated ERA5 2D Wave Spectrum, E(f,theta)");
end

% INTERPOLATE/REGRID
disp('*==== Start Interpolation/Regridding ====*');
interp_size = 128;

[S_kxky, kx , ky] = wavenumspec(plotsON, era5_2D, adjusted_dir_degrees, era5_f, interp_size, sar_dx_range);

% Compute Hs in kx ky space
kappa = sqrt(kx.^2 + ky.^2);
S_2Dxy = S_kxky./kappa;
Hs_kxky = 4 * sqrt(trapz(ky(:,1), trapz(kx(1,:), S_2Dxy, 2)));
disp(['In kx-ky equispaced, Hs: ', num2str(Hs_kxky)]);

% Plot Cartesian spectrum
if plotsON
    figure('Position', [0, 0, 800, 300]);
        subplot(1,2,1);
        contour(kx,ky, kappa); title("Kappa = Kx^2 + Ky^2");
        subplot(1,2,2);
        contour(kx, ky, S_kxky); title("Wave Num Spec on SAR grid", "E(k)")
        xlim([min(kx(:)) max(kx(:))]);
        ylim([min(ky(:)) max(ky(:))]); colorbar;
end

%1.3 First Guess Wave Number Spectra Variables
%These final variables are redundant but helpful for now in the process of understanding things.
first_guess_wave_number_spectrum = S_kxky;
first_guess_kx_range = kx;
first_guess_ky_azimuth = ky;
first_guess_k = kappa;
first_guess_omega = sqrt(g .* first_guess_k); % [Engen, 1995]
%% =============H&H: Generate First Guess (Simulated) SAR Spectrum ====================
disp('*==== Start SAR Spec Generation ====*');
disp(['Nonlinearity order: ', num2str(nonlinearity_order)]);

[first_guess_sar_spectrum, ~, ~, ~, ~] = handh.generateSARSpectrumFromWaveNumberSpectrum(sentinel, plotSARSpecGen, 1, nonlinearity_order, imagette_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k,first_guess_wave_number_spectrum);

%% ==================== H&H: Observed SAR Spectrum ====================
disp('*==== Start SAR Spec Observed ====*');
cutoff_wavelength = wavelength_mpts + 10; % Specify the length of the longest expected wave in the region
size_of_filter_window = 3; % The number of elements = dim_filt + 1
width_of_gaussian_lobe = 3; % width_of_gaussian_lobe: Controls the spread of the Gaussian. A larger width means the Gaussian will be wider and smoother.
[observed_sar_spectrum] = observedSARSpectrum(1,sar_transect_data, sar_imagette_size, imagette_size, sar_dx_range, first_guess_k, size_of_filter_window, width_of_gaussian_lobe, cutoff_wavelength,first_guess_kx_range ,first_guess_ky_azimuth );

if plotsON
    figure('Position', [0, 0, 1200, 300]);
    subplot(1,3,1); plotLibrary.generalSpectrumPlots(0,first_guess_wave_number_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "First guess Wave Number Spectrum, E(k)","?"); 
    subplot(1,3,2); plotLibrary.generalSpectrumPlots(0,first_guess_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Simulated SAR Spectrum, P^S_k","?");
    subplot(1,3,3); plotLibrary().generalSpectrumPlots(0,observed_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth,  ["SAR Spectrum, P^S_{obs}"," with Butterworth cutoff and Gaussian Filter"],"?");
end

%% ==================== H&H: Step Zero ====================
transect_number = 1; % This will mainly be used when looping through all transects in the ice part of the pipeline

if  applyStepZero
    P_obs = observed_sar_spectrum;

    % START OF ALGO
    first_iteration_check = true;
    
    rotation_angles = -5:1:5;
    wave_number_mag_scales = 0.9:0.1:1.1;
    energy_scales = 0.9:0.5:1.5;

    J_step0 = zeros(length(rotation_angles), length(wave_number_mag_scales), length(energy_scales));
    
    % Extrapolate so that the grid space is larger for adjustment
    [S_tmp_0, kx_tmp,ky_tmp] = wavenumspecStepzero(plotsON, era5_2D, adjusted_dir_degrees, era5_f, interp_size, sar_dx_range);% Convert to Cartesian (kx, ky) coordinates  
    
    for i_rot = 1:length(rotation_angles)
        disp(num2str(i_rot));

        for i_mag = 1:length(wave_number_mag_scales)
         
            for i_energy = 1:length(energy_scales)

                % Reset the wave number spec
                S_tmp = S_tmp_0;

                % [Eq.78 HH1991]
                step0_B = wave_number_mag_scales(i_mag);
                step0_kx = step0_B .* (kx_tmp .* cosd(rotation_angles(i_rot)) - ky_tmp .* sind(rotation_angles(i_rot)));
                step0_ky = step0_B .* (kx_tmp .* sind(rotation_angles(i_rot)) + ky_tmp .* cosd(rotation_angles(i_rot)));
                
                % [Eq.78 HH1991]
                A = energy_scales(i_energy);
                S_tmp = energy_scales(i_energy) .* S_tmp;
                
                F = scatteredInterpolant(step0_kx(:), step0_ky(:), S_tmp(:), 'linear', 'linear');
                S_tmp_n = F(kx, ky);
                S_tmp_n(isnan(S_tmp_n)) = 0;
               
                % Calculate the new SAR spectrum
                [P_tmp, ~, ~, ~, ~] = handh.generateSARSpectrumFromWaveNumberSpectrum(sentinel, plotsON, transect_number, nonlinearity_order, imagette_size, kx, ky, first_guess_omega, first_guess_k, S_tmp_n); % [Eq.65 HH1991]
                P_tmp(isnan(P_tmp)) =0;

                kx_vector = kx(1,:);
                ky_vector = ky(:,1);
                
                % H&H Eq.63 cost-function
                B = 0.01 * max(first_guess_wave_number_spectrum(:))^2; % Below [Eq.76 HH1991]
                mu = 0.1 * max(P_obs(:))^2; % [Eq.76 HH1991]
                term1_63 = (P_tmp - P_obs).^2 .* P_obs;
                term1_63(isnan(term1_63))=0;
                term2_63 = mu .* ( (S_tmp_n - first_guess_wave_number_spectrum) ./ (B + first_guess_wave_number_spectrum) ).^2;
                term2_63(isnan(term2_63))=0;

                % FOR REPORT, PLOT THE TERMS FOR DEMONSTRATION REASONS
                % figure; 
                % subplot(1,2,1); pcolor(term1_63);
                % subplot(1,2,2); pcolor(term2_63);
                
                J_step0(i_rot, i_mag, i_energy) = trapz( kx_vector,trapz(ky_vector,term1_63,1),2) + trapz( kx_vector,trapz(ky_vector,term2_63,1),2);% dimensions = Y x X = 1 x 2;
               
                % Giacomo version
                % J_step0(i_rot, i_mag, i_energy) = trapz(kx_vector, trapz(ky_vector, (P_tmp - P_obs).^2 .* anello,1),2) ./ sqrt(trapz(kx_vector, trapz(ky_vector, (P_obs).^2 .* anello,1),2)) ./ sqrt(trapz(kx_vector, trapz(ky_vector, (P_tmp).^2  .* anello,1),2));
                
                if first_iteration_check
                    Jmax_step0 = J_step0(i_rot, i_mag, i_energy);
                    first_iteration_check = false;
                end
                    
                if J_step0(i_rot, i_mag, i_energy) <= Jmax_step0
                    disp("iteration "+i_rot +"," + i_mag+","+i_energy)
                    Jmax_step0      = J_step0(i_rot, i_mag, i_energy);
                    F_new           = S_tmp_n;
                    P_new           = P_tmp;
                    energy_0        = energy_scales(i_energy);
                    mag_0           = wave_number_mag_scales(i_mag);
                    rot_angle_0     = rotation_angles(i_rot);
                end
            
            end

        end
    end
    
    disp("Selected: rot= "+num2str(rot_angle_0)+ ...
        "; wave#mag = "+num2str(mag_0) +  ...
        "; energy = "+num2str(energy_0));

    % Plot the outputs
    figure('Position', [0, 0, 1600, 800]);
        subplot(2,4,1); plotLibrary.generalSpectrumPlots(0,observed_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Observed SAR spectrum","?");
        subplot(2,4,2); plotLibrary.generalSpectrumPlots(0,first_guess_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Original SAR spectrum","?");
        subplot(2,4,3); plotLibrary.generalSpectrumPlots(0,P_new, kx, ky, "Step Zero SAR Spectrum, P(k)","?");
        subplot(2,4,4);
        % ix
        % countJ=J(i_rot,)
            % pcolor(J_step0(:));
            % clear temp; temp(:,:) = J_step0(:,1,:); surf(temp); view(2);
            % title("Minimisation Values"); ylabel("Wave Number Scaling Factor"); xlabel("Rotation angle (degrees)"); 
            % c = colorbar(); c.Label.String = 'Minimisation magnitude';
        subplot(2,4,6); plotLibrary.generalSpectrumPlots(0,first_guess_wave_number_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Original Wave Number Spectrum","?");
        subplot(2,4,7); plotLibrary.generalSpectrumPlots(0,F_new, kx, ky, "Step Zero Wave Number Spectrum, E(k)","?");
    
    % Update the variables to the new values
    first_guess_wave_number_spectrum = F_new;
end
%% ==================== H&H: Inversion ====================
disp('*==== Start Inversion ====*');
[J_eq_62,P_best_eq62,F_best_eq62,J_eq_63,P_best_eq63,F_best_eq63,J_eq_69,P_best_eq69,F_best_eq69] = inversionHH(inversion_iterations, nonlinearity_order, sentinel, handh, 0, observed_sar_spectrum, first_guess_wave_number_spectrum, transect_number, imagette_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k );
%TO DO: ADD MARKER TO SHOW WHICH ITERATION WAS CHOSEN
if inversionPlots
    figure('Position', [0, 0, 1200, 900]);
    subplot(3,3,1); plotLibrary.generalSpectrumPlots(0,F_best_eq62, first_guess_kx_range, first_guess_ky_azimuth, ["Final Wave Spectrum, F_k","From Eq.62 Inversion"],"?");
    subplot(3,3,2); plotLibrary.generalSpectrumPlots(0,P_best_eq62, first_guess_kx_range, first_guess_ky_azimuth, ["Final SAR Spectrum, P^S_{k}"," From Eq.62 Inversion"],"?");
    subplot(3,3,3); 
        plot(J_eq_62); xticks((1:1:inversion_iterations)); xlabel("Number of iterations");
    subplot(3,3,4); plotLibrary.generalSpectrumPlots(0,F_best_eq63, first_guess_kx_range, first_guess_ky_azimuth, ["Final Wave Spectrum, F_k"," From Eq.63 Inversion"],"?");
    subplot(3,3,5); plotLibrary.generalSpectrumPlots(0,P_best_eq63, first_guess_kx_range, first_guess_ky_azimuth, ["Final SAR Spectrum, P^S_{k}"," From Eq.63 Inversion"],"?");
    subplot(3,3,6);
        plot(J_eq_63); xticks((1:1:inversion_iterations)); xlabel("Number of iterations");
    subplot(3,3,7); plotLibrary.generalSpectrumPlots(0,F_best_eq69, first_guess_kx_range, first_guess_ky_azimuth, ["Final Wave Spectrum, F_k"," From Eq.69 Inversion"],"?");    
    subplot(3,3,8); plotLibrary.generalSpectrumPlots(0,P_best_eq69, first_guess_kx_range, first_guess_ky_azimuth, ["Final SAR Spectrum, P^S_{k}","From Eq.69 Inversion"],"?");    
    subplot(3,3,9); 
        plot(J_eq_69); xticks((1:1:inversion_iterations)); xlabel("Number of iterations");
end

chosen_final_spectrum = F_best_eq63;
chosen_final_SAR_spectrum = P_best_eq63;

%% ==================== RESULTS ====================
disp('*==== OUTPUTS ====*');
% Convert E(k) - E(k,theta)
kappa = sqrt(kx.^2 + ky.^2);
% Compute Hs in kx-ky space
Hs_kxky = 4 * sqrt(trapz(ky(:,1), trapz(kx(1,:), chosen_final_spectrum./kappa, 2)));
disp(['In SAR kx-ky equispaced, Hs: ', num2str(Hs_kxky)]);

% Recovered variables
% dimension of this grid can be selected as preferred (in principle can
% have any dimension)
n = interp_size;
k_r=logspace(log10(min(abs(kappa(:)))),log10(max(kappa(:))),n);
theta_r=linspace(0,2*pi,n);
% remove SAR az to North angle
theta_r = (handh.getDirectionsInSARGeometry(theta_r*(180/pi),sar_az_to_trueN_angle))*pi/180;

for i=1:length(k_r)
    for j=1:length(theta_r)
        kx_r(i,j)=k_r(i)*sin(theta_r(j));
        ky_r(i,j)=k_r(i)*cos(theta_r(j));
    end
end

% Interpolate onto new grid
F = scatteredInterpolant(kx(:), ky(:), chosen_final_spectrum(:), 'natural', 'none');
S_r = F(kx_r, ky_r);
S_r(isnan(S_r))=0;

kappa=sqrt(kx_r.^2+ky_r.^2);

% Compute Hs in kx-ky space
S_rxy=S_r./kappa;
Hs_kxky = 4 * sqrt(trapz(kx_r(:,1), trapz(ky_r(1,:), S_rxy)));
disp(['In kx-ky polar, Hs: ', num2str(Hs_kxky)]);
clear S_rxy;

Hs_m0 = 4 * sqrt(trapz(k_r, trapz(theta_r, S_r, 2), 1));
disp(['(after regrid) In k theta space, Hs: ', num2str(Hs_m0)]);

% INVERSE_TRANSFORM_SPECTRUM Converts a wave spectrum from wavenumber (k) to frequency (f)
f_r = sqrt(k_r * g) / (2 * pi);
dk_df = (8 * pi^2 / g) * f_r';
S_f_r = S_r .* dk_df;

% Plot recovered 2D spectrum
figure('Position', [0, 0, 1200, 300]);
    subplot(1,3,1);
    contour(kx, ky, chosen_final_spectrum); title("Final chosen Spectrum", "on SAR kx,ky")
    xlim([min(kx(:)) max(kx(:))]);
    ylim([min(ky(:)) max(ky(:))]);
    subplot(1,3,2);
    contour(kx_r, ky_r, S_r); title("Final spectrum regridded to new kx_r, ky_r")
    xlim([min(kx(:)) max(kx(:))]);
    ylim([min(ky(:)) max(ky(:))]);
    subplot(1,3,3);
    plotLibrary().waveSpectrum2D(1, abs(S_f_r'), f_r, theta_r * 180 / pi, 1, "Recovered 2D wave spectrum");

% Compute recovered wave height
S_f_r(isnan(S_r)) = 0;
Hs_m0 = 4 * sqrt(trapz(f_r, trapz(theta_r, S_f_r, 2), 1));
disp(['Recovered Hs: ', num2str(Hs_m0)]);

% Calculate recovered 1D spectrum
S_1Dr = trapz(theta_r, S_f_r, 2);
%% ==================== FINAL PLOT  ====================
figure('Position', [0, 0, 600, 1200]);
    subplot(4,2,1); plotLibrary.waveSpectrum2D(1,era5_2D, era5_f, era5_dir_degrees, 0.3, ["ERA5 2D Wave Spectrum, E(f,theta)",string(era5_datetime)]);
    subplot(4,2,2); plotLibrary.generalSpectrumPlots(0,observed_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "Observed SAR spectrum","?");
    
    subplot(4,2,3); plotLibrary.generalSpectrumPlots(0,first_guess_wave_number_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "First Guess Wave Number Spectrum","?");
    subplot(4,2,4); plotLibrary.generalSpectrumPlots(0,first_guess_sar_spectrum, first_guess_kx_range, first_guess_ky_azimuth, "First Guess SAR spectrum","?");
    
    subplot(4,2,5); plotLibrary.generalSpectrumPlots(0,chosen_final_spectrum, first_guess_kx_range, first_guess_ky_azimuth, ["Final Wave Number Spectrum, F_k","Eq.63 Inversion"],"?");
    subplot(4,2,6); plotLibrary.generalSpectrumPlots(0,chosen_final_SAR_spectrum, first_guess_kx_range, first_guess_ky_azimuth,["Final SAR Spectrum, F_k","Eq.63 Inversion"],"?");
   
    subplot(4,2,7); plotLibrary().waveSpectrum2D(1, abs(S_f_r'), f_r, theta_r * 180 / pi, 1, ["Final 2D wave spectrum","E(f,theta)"]);

%% ==================== COMPARISON ====================
disp('*==== Comparison ====*');
[final_Hs_derived,final_Tm_derived,final_direction_derived,final_total_variance_or_energy_derived] = waveLibrary.calculateSpectrumCharacteristics(theta_r,f_r,S_f_r');
disp(['Final Hs: ',  num2str(final_Hs_derived)]);
disp(['Final Tm: ',  num2str(final_Tm_derived)]);
disp(['Final Dir: ',  num2str(final_direction_derived +180)]);

[~,max_index] = max(S_1Dr);
final_Tp = 1/f_r(max_index);
disp(['Final Tp: ', num2str(final_Tp)]);
wavelength_Tp_final = ((final_Tp)^2 * 9.81)/(2*pi);
disp(['Final wavelength from Tp: ', num2str(wavelength_Tp_final)]);

disp('--------------------------------------');
disp("ERA5 date"); disp(era5_datetime)
disp(['ERA5 lon: ',  num2str(era5_lon)]);
disp(['ERA5 lat: ',  num2str(era5_lat)]);
[era5_Hs_derived,era5_Tm_derived,era5_direction_derived,~] = waveLibrary.calculateSpectrumCharacteristics(era5_dir_degrees.*(pi/180),era5_f,era5_2D);

disp(['Inital ERA5 Hs: ', num2str(era5_Hs_derived)]);
disp(['Inital ERA5 Tm: ', num2str(era5_Tm_derived)]);
disp(['Inital ERA5 Dir: ', num2str(era5_direction_derived +180)]);

[~,max_index] = max(S_1D);
Tp = 1/era5_f(max_index);
disp(['ERA5 Tp: ', num2str(Tp)]);
wavelength_Tp = ((Tp)^2 * 9.81)/(2*pi);
disp(['ERA5 wavelength from Tp: ', num2str(wavelength_Tp)]);

disp('--------------------------------------');

disp("Singles ERA5 date"); disp(era5_datetime_single)
disp(['Singles ERA5 Hs: ', num2str(era5_Hs)]);
disp(['Singles ERA5 Tm: ', num2str(era5_Tm)]);
disp(['Singles ERA5 Dir: ', num2str(era5_direction)]);
disp(['Singles ERA5 mdts: ', num2str(era5_mdts)]);
disp(['Singles ERA5 mpts: ', num2str(era5_mpts)]);
disp(['Singles ERA5 shts: ', num2str(era5_shts)]);
disp(['Singles ERA5 mdww: ', num2str(era5_mdww)]);
disp(['Singles ERA5 mpww: ', num2str(era5_mpww)]);
disp(['Singles ERA5 shww: ', num2str(era5_shww)]);

wavelength_mww = ((era5_mpww)^2 * 9.81)/(2*pi);

disp(['Singles ERA5 wavelength from mpww: ', num2str(wavelength_mww)])
disp(['Singles ERA5 wavelength from mpts: ', num2str(wavelength_mpts)])
  
disp('--------------------------------------');
minutes = minute(era5_datetime);
if minutes<10
    minutes_str = "0"+num2str(minute(era5_datetime));
else
    minutes_str = num2str(minute(era5_datetime));
end
time_str = num2str(hour(era5_datetime)) +""+minutes_str;
date_time_str = date_str + ""+ time_str;
% Cape Point Buoy Data
filepath_buoy = "/Users/tris/Documents/MSc/data/Cape Point Buoy/CapePointBuoy_"+date_time_str+".mat";
buoy_data = load(filepath_buoy);

disp("Buoy"); disp(buoy_data.date)
disp(['Buoy Hs: ', num2str(buoy_data.Hs)]);
disp(['Buoy Tp: ', num2str(buoy_data.Tp)]);
disp(['Buoy Dir: ', num2str(buoy_data.Dir)]);

if minutes_str=="00"
    minutes_str_2="30";
else
    minutes_str_2="00";
end
time_str_2 = num2str(hour(era5_datetime)) +""+minutes_str_2;
date_time_str = date_str + "" + time_str_2;
% Cape Point Buoy Data
filepath_buoy = "/Users/tris/Documents/MSc/data/Cape Point Buoy/CapePointBuoy_"+date_time_str+".mat";
buoy_data_2 = load(filepath_buoy);


disp("Buoy"); disp(buoy_data_2.date)
disp(['Buoy Hs: ', num2str(buoy_data_2.Hs)]);
disp(['Buoy Tp: ', num2str(buoy_data_2.Tp)]);
disp(['Buoy Dir: ', num2str(buoy_data_2.Dir)]);


%% 
% figure('Position', [0, 0, 800, 300]);
figure;
    % subplot(1,2,1);
    plot(f_r, S_1Dr, 'k',"LineWidth",1, 'DisplayName', 'Recovered');
    hold on;
    plot(era5_f, S_1D, 'r',"LineWidth",1, 'DisplayName', 'ERA5');
    % plot(buoy_data.frequencies,buoy_data.energy,"Color","#0072BD","LineWidth",1, 'DisplayName', "Buoy");% [" + time_str +"]");

    plot(buoy_data_2.frequencies,buoy_data_2.energy,"Color",'#0072BD', 'DisplayName', "Buoy");%[" + time_str_2 +"]");
    legend;
    % title('1D Wave Spectra');  xlabel("Frequency [Hz]"); ylabel("E(f) [m^{2}/Hz]")
    
    % subplot(1,2,2);
    % plot(f_r, S_1Dr, 'r', 'DisplayName', 'Recovered');
    % hold on;
    % plot(era5_f, S_1D, 'k', 'DisplayName', 'ERA5');
    %  hold on;
    % plot(buoy_data.frequencies,buoy_data.energy,'b', 'DisplayName', "Buoy [" + time_str +"]");
    % hold on;
    % plot(buoy_data_2.frequencies,buoy_data_2.energy,'m', 'DisplayName', "Buoy [" + time_str_2 +"]");
    % legend;
    % title('Original, Recovered, Buoy 1D Spectrum');  xlabel("Frequency [Hz]"); ylabel("Energy [m^2/Hz]");
    % set(gca,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
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
    kx = 3*linspace(-pi/dx,pi/dx,3*n);
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