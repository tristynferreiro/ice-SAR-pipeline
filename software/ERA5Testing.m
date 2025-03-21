%% ==================== INITIALISATION ====================
% Import libraries
waveLibrary = waveLibrary;
plotLibrary = plotLibrary;

% ERA5 DATA
filepath = "/Users/tris/Documents/MSc/Data/Cape Point/CapePoint_ERA5-2Dws_20221014.nc";

era5 = ERA5(filepath);
testMultipleKTheory(era5);

%% 
% Get attributes
era5_longitudes = era5.getLongitude;
era5_latitudes = era5.getLatitude;
era5_time = era5.getTime;
% era5_d2fd_all_5dims = era5.getAllWaveSpectrumD2FD;

era5_freq_bins = era5.FrequencyBins;

era5_direction_bins_degrees = era5.DirectionBins;

sar_start_datetime = datetime("14-Oct-2024 17:34:32");
sar_transect_center_latitude = -34.204; %-32.204 has nice spectrum for Report!!!
sar_transect_center_longitude = 18.28666944;
[era5_d2fd, era5_lat, era5_lon, era5_datetime] = era5.getSlicedWaveSpectrumD2Fd(sar_transect_center_latitude,sar_transect_center_longitude,sar_start_datetime);

f = era5_freq_bins;
theta_degrees = era5_direction_bins_degrees;
S_2D = era5_d2fd;
interp_size = 128;
%% ==================== Testing ====================
    
    % Wave Number Spectrum

    % Interpolation
    testInterpolation(f,theta_degrees,S_2D,interp_size);
        
    % Calculate 1D Spectrum and compare before and after interpolation
    calculate1DSpectra(f, theta_degrees, S_2D, interp_f, interp_theta_degrees, interp_S_2D);
    
    % Wave Number Spectrum with Interpolation
    testWaveNumberSpectrumWithInterpolation(f,theta_degrees,S_2D,interp_size);
    
    % Azimuth angle adjustment
    angle = 15.*pi/180;
    testRotation(f,theta_degrees.*(pi/180),S_2D,angle, [Hs_derived,Tm_derived,mean_direction_derived_degrees], 1);
    
    % Azimuth angle adjustment with interpolation
    % 1. Rotate the spectrum and then interpolate, compare this to the original interpolation.
    % 2. Interpolate the spectrum and then rotate, compare this to other angle adjustment.
   
 
%% FUNCTIONS
function f = frequencyBins(numberOfFrequencyBins, Type)
    if(Type == "ERA5")
        if numberOfFrequencyBins<=30
            % Setup the frequency array as per the ERA5 documentation.
            f = zeros(numberOfFrequencyBins,1);
            f(1) = 0.03453;
            for i=2:length(f)
                f(i) = f(i-1)*1.1;
            end
        else 
            % Need to preserve the statistical property regardless of the
            % size of the array. i.e. it still needs to mirror the case as
            % if there were only 30 bins (as we expect in ERA5 data).
            max_f = 0.5478; % Calculated from the 30 bin case
            % Setup the frequency array as per the ERA5 documentation.
            f = zeros(numberOfFrequencyBins,1);
            f(1) = 0.03453;
            r = (max_f / f(1))^(1 / (numberOfFrequencyBins - 1));

            for i=2:length(f)
                f(i) = f(i-1)*r;
            end
        end
    elseif (Type == "linear")
        % Generic array
         f = linspace(0.01,1,numberOfFrequencyBins)';
    end
end

function [theta_degrees, theta_rad] = directionBins(numberOfDirectionBins, mean_wave_direction_degrees, Type)

    if(Type == "ERA5")
        max_degree_bin = 352.5; % ERA5 maximum
        min_degree_bin = 7.5;
        
        % bin_spacing = (max_degree_bin-min_degree_bin)/numberOfDirectionBins; % Calculate the grid spacing
        % 
        % 
        % theta_degrees = linspace(mean_wave_direction_degrees-(bin_spacing*numberOfDirectionBins/2), ...
        %                                     mean_wave_direction_degrees+(bin_spacing*numberOfDirectionBins/2) ...
        %                                     ,numberOfDirectionBins);

        theta_degrees = linspace(min_degree_bin,max_degree_bin,numberOfDirectionBins);
        theta_rad = theta_degrees.* (pi/180);

    elseif (Type == "full circle")
        max_degree_bin = 360; % ERA5 maximum
        min_degree_bin = 0;
        
        bin_spacing = (max_degree_bin-min_degree_bin)/numberOfDirectionBins; % Calculate the grid spacing
        
        
        theta_degrees = linspace(mean_wave_direction_degrees-(bin_spacing*numberOfDirectionBins/2), ...
                                            mean_wave_direction_degrees+(bin_spacing*numberOfDirectionBins/2) ...
                                            ,numberOfDirectionBins);
        theta_rad = theta_degrees.* (pi/180);
    end

end

function testRotation(f,theta_rad,S_2D,rotation_angle_rad, SimInputParams, showPlots)
    %TESTROTATION Apply the rotation angle and then remove it, compare the 
    % spectrum parameters at each point and output.
        % S_2D              = E(f, theta) where theta is in radians
        % f                 = frequency array
        % theta_rad         = direction array in radians
        % rotationAngle     = angle to rotate the spectum by in radians
        % SimInputParams    = [Hs,T0,meanWaveDirectionDegrees]

    adjusted_theta_rad = theta_rad + rotation_angle_rad;

    %COMPARE spectra with angle adjustment
    if showPlots
        figure('Position', [0, 0, 1200, 300]);
        subplot(1,3,1); plotLibrary().waveSpectrumAdjustedToSAR(S_2D, f, theta_rad.*(180/pi), 0.3, "Original 2-D Wave Spectrum, E(f,theta)");
        subplot(1,3,2); plotLibrary().waveSpectrumAdjustedToSAR(S_2D, f, adjusted_theta_rad.*(180/pi), 0.3, "Adjusted to SAR scene geometry 2-D Wave Spectrum, E(f,theta)");
        subplot(1,3,3);
            plotLibrary().waveSpectrumAdjustedToSAR(S_2D, f, adjusted_theta_rad.*(180/pi), 0.3, "");
            hold on;
            plotLibrary().waveSpectrumAdjustedToSAR(S_2D, f, theta_rad.*(180/pi), 1, "Comparison of original vs adjusted spectrum");
    end
       
    
    % COMPARE the parameter values
    [Hs_derived,Tm_derived,direction_derived,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(adjusted_theta_rad,f,S_2D);

    readjusted_theta_rad = adjusted_theta_rad - rotation_angle_rad;
    [Hs_derived_readj,Tm_derived_readj,direction_derived_readj,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(readjusted_theta_rad,f,S_2D);
    
    % Display Test results
    disp("--------------------------------------------")
    disp("         Test 2D Spectrum Rotation          ")
    disp("--------------------------------------------")
    disp("Param | Input | Rotated 2D | Un-rotated 2D ")
    disp("--------------------------------------------")
    disp([' Hs   |   ', num2str(SimInputParams(1)), '   |     ', num2str(Hs_derived), '      |     ',num2str(Hs_derived_readj)])
    disp([' Tm   |   ',num2str(SimInputParams(2)),'    |    ', num2str(Tm_derived),'  |   ', num2str(Tm_derived_readj)]);
    disp([' Dir_m |   ', num2str(SimInputParams(3)), '  |     ', num2str(direction_derived), '     |     ', num2str(direction_derived_readj)])

end

function calculate1DSpectra(f, theta_degrees, S_2D, interp_f, interp_theta_degrees, interp_S_2D)
         % Calc parameters of 1D to check that it matches the input
        [S_1D, f] = waveLibrary().calc1DWaveSpec(f, theta_degrees.*(pi/180), S_2D);
        [Hs_derived_1D,Tm_derived_1D] = waveLibrary().calculate1DSpectrumCharacteristics(f,S_1D);
    
        [interp_S_1D, interp_f] = waveLibrary().calc1DWaveSpec(interp_f, interp_theta_degrees.*(pi/180), interp_S_2D);
        [Hs_derived_1D_interp,Tm_derived_1D_interp] = waveLibrary().calculate1DSpectrumCharacteristics(interp_f,interp_S_1D);
    
        % Display Test results
        disp("------------------------------")
        disp("      Test 1D Spectrum        ")
        disp("------------------------------")
        disp("Param | Original 1D | Interp 1D")
        disp("------------------------------")
        disp([' Hs   |   ', num2str(Hs_derived_1D), '   |     ', num2str(Hs_derived_1D_interp)])
        disp([' Tm   |   ', num2str(Tm_derived_1D), '  |     -', num2str(Tm_derived_1D_interp)]);
       
        figure('Position', [0, 0, 800, 300]);
            subplot(1,2,1); plot(f,S_1D); 
                xlabel("frequency"); ylabel("m^2/Hz"); title("Original 1D Spectrum");
            subplot(1,2,2); plot(interp_f,interp_S_1D);
                xlabel("frequency"); ylabel("m^2/Hz"); title("Interpolated 1D Spectrum");

    end

function [Hs_k,Tm_k,direction_k] = testCalculationOfSpectralParamFromEk(kx,ky,E_kx_ky)

    m0_k = trapz(kx,trapz(ky,E_kx_ky,1),2);
    Hs_k = 4 * sqrt(m0_k);
    disp(Hs_k)

    % intermediate = trapz(direction_bins_rad,(1./frequency_bins) .* E_kx_ky,1);
    % intermediate(isnan(intermediate))=0;
    % m_neg1 = trapz(frequency_bins,intermediate,2); %[Eq.8]
    % 
    % intermediate = trapz(direction_bins_rad,(frequency_bins) .* E_f_theta,1);
    % intermediate(isnan(intermediate))=0;
    % m1 = trapz(frequency_bins,intermediate,2); %[Eq.8]
    % 
    % % Calculate the mean wave period (Tm) [Eq.5]
    % % Tm = m_neg1/m0; % [Eq.5]
    % Tm = (m1/m0)^(-1); % [Holthuijsen Eq.5]

    Tm_k = 0;

    direction_k = 0;
end

function [Hs,Tm,mean_direction,Hs_test, Tm_test, mean_direction_test] = testKTheory(plotsON,frequency_bins, direction_bins_degrees,S_2D,azimuth_angle)
    % figure; contour(frequency_bins,direction_bins_degrees,S_2D);
    if plotsON
        figure; plotLibrary().waveSpectrum2D(1, S_2D,frequency_bins, ...
                direction_bins_degrees, 0.3, "Original S2D");
    end


    [Hs,Tm, mean_direction] = waveLibrary().calculateSpectrumCharacteristics(direction_bins_degrees.*(pi/180),frequency_bins,S_2D);

    w = (2*pi).*frequency_bins;
    gravity = 9.81;
    k = (w).^2 ./ gravity;
    
    [E_kx_ky, Jacobian_kx_ky, kx_matrix, ky_matrix] = waveLibrary().waveNumberSpectrum(S_2D,w,k,direction_bins_degrees.*(pi/180), mean_direction*pi/180, azimuth_angle);
    
    if plotsON
        figure; contour(kx_matrix,ky_matrix,E_kx_ky)
    end
    
    kx_matrix_lin = linspace(min(kx_matrix(:)), max(kx_matrix(:)), 128);
    ky_matrix_lin = (linspace(min(ky_matrix(:)), max(ky_matrix(:)), 128))';
    k_test = sqrt(kx_matrix_lin.^2 + ky_matrix_lin.^2);

    % [xq,yq] = meshgrid(kx_matrix_test, ky_matrix_test);
    % kx_matrix_test = xq;
    % ky_matrix_test = yq;

    E_kx_ky_lin = griddata(kx_matrix,ky_matrix,E_kx_ky,kx_matrix_lin,ky_matrix_lin);
    
    if plotsON
    figure; contour(kx_matrix_lin,ky_matrix_lin,E_kx_ky_lin);
    end
    

    E_kx_ky_test = griddata(kx_matrix_lin,ky_matrix_lin,E_kx_ky_lin,kx_matrix,ky_matrix);
    if plotsON
        figure; contour(kx_matrix,ky_matrix,E_kx_ky_test);
    
        figure;
            contour(kx_matrix,ky_matrix,E_kx_ky_test); hold on;contour(kx_matrix,ky_matrix,E_kx_ky)
    end
    
    [S_2D_test,Jacobian_f_theta] = waveLibrary().waveSpectrum(E_kx_ky_test, w, k);
    S_2D_test(isnan(S_2D_test)) = 0;

    if plotsON
        figure; contour(frequency_bins,direction_bins_degrees,S_2D_test);
        figure; plotLibrary().waveSpectrum2D(1, S_2D_test,frequency_bins, ...
                direction_bins_degrees, 0.3, "Test S2d");
        %

        figure; plotLibrary().waveSpectrum2D(1, S_2D_test,frequency_bins, ...
                direction_bins_degrees, 0.3, "Test S2D");
        hold on; plotLibrary().waveSpectrum2D(1, S_2D,frequency_bins, ...
            direction_bins_degrees, 0.3, "Original S2D");
    end
    
    [Hs_test,Tm_test, mean_direction_test] = waveLibrary().calculateSpectrumCharacteristics(direction_bins_degrees.*(pi/180),frequency_bins,S_2D_test);


    % Display Test results
    % disp("------------------------------")
    % disp("      Test K Theory       ")
    % disp("------------------------------")
    % disp("Param | Original  | Test")
    % disp("------------------------------")
    % disp([' Hs   |   ', num2str(Hs), '   |     ', num2str(Hs_test)])
    % disp([' Tm   |   ', num2str(Tm), '  |     ', num2str(Tm_test)]);
    % disp([' Dir_m   |   ', num2str(mean_direction), '  |     ', num2str(mean_direction_test)]);
    % 
    % % h_diff = (Hs-Hs_test)/Hs_test;
    % % t_diff = (Tm-Tm_test)/Tm_test;
    % % dir_diff = (mean_direction-mean_direction_test)/(352.5-7.5);
    % % disp(h_diff);
    % % disp(t_diff)
    % % disp(dir_diff)
    % % 
    
    
    
    % % Giacomo line
    % sar_dy_azimuth =10;
    % dk = 2 * pi/(128 * sar_dy_azimuth);
    % kmin = -pi/sar_dy_azimuth + dk;
    % kmax = pi/sar_dy_azimuth;
    % k_giacomo = kmin:dk:kmax;
    % figure; plot(k_giacomo);
    % 
    % kx_matrix_giacomo = k_giacomo .* sind(theta);
    % ky_matrix_giacomo = k_giacomo' .* cosd(theta);
    % 
    % [xq,yq] = meshgrid(kx_matrix_giacomo, ky_matrix_giacomo);
    % kx_matrix_giacomo = xq;
    % ky_matrix_giacomo = yq;
    % 
    % E_kx_ky_giacomo = griddata(kx_matrix,ky_matrix,E_kx_ky,kx_matrix_giacomo,ky_matrix_giacomo);
    % figure; contour(kx_matrix_giacomo,ky_matrix_giacomo,E_kx_ky_giacomo)
    % 
    % E_kx_ky_test = griddata(kx_matrix_giacomo,ky_matrix_giacomo,E_kx_ky_giacomo,kx_matrix,ky_matrix);
    % figure; contour(kx_matrix,ky_matrix,E_kx_ky_test);
    % 
    % [S_2D_test,Jacobian_f_theta] = waveLibrary().waveSpectrum(E_kx_ky_test, w, k);
    % S_2D_test(isnan(S_2D_test)) = 0;
    % figure; contour(frequency_bins,direction_bins_degrees,S_2D_test);
    % figure; plotLibrary().waveSpectrum2D(1, S_2D_test,frequency_bins, ...
    %         direction_bins_degrees, 0.3, "Test S2d");
    % %
    % figure; plotLibrary().waveSpectrum2D(1, S_2D_test,frequency_bins, ...
    %         direction_bins_degrees, 0.3, "Test S2D");
    % hold on; plotLibrary().waveSpectrum2D(1, S_2D,frequency_bins, ...
    %         direction_bins_degrees, 0.3, "Original S2D");
    % 
    % figure; contour(frequency_bins,direction_bins_degrees,S_2D_test);
    % hold on;
    % contour(frequency_bins,direction_bins_degrees,S_2D);
    % %
    % [Hs_test,Tm_test, mean_direction_test] = waveLibrary().calculateSpectrumCharacteristics(direction_bins_degrees.*(pi/180),frequency_bins,S_2D_test);
    % 
    % 
    % % Display Test results
    % disp("------------------------------")
    % disp("      Test K Theory       ")
    % disp("------------------------------")
    % disp("Param | Original  | Test")
    % disp("------------------------------")
    % disp([' Hs   |   ', num2str(Hs), '   |     ', num2str(Hs_test)])
    % disp([' Tm   |   ', num2str(Tm), '  |     -', num2str(Tm_test)]);
    % disp([' Dir_m   |   ', num2str(mean_direction), '  |     -', num2str(mean_direction_test)]);
   
end

function testInterpolation(f,theta_degrees,S_2D,interp_size)
         % Calculate parameters of original spectrum
        [Hs_derived,Tm_derived,mean_direction_derived_degrees,~] = ...
            waveLibrary().calculateSpectrumCharacteristics(theta_degrees.*(pi/180),f,S_2D);
        
        % Interpolate the spectrum
        [interp_S_2D, interp_f, interp_theta_degrees] = waveLibrary().interpolateWaveSpectrum(interp_size, f,theta_degrees,S_2D,"ERA5");
        
        % testKTheory(interp_f, interp_theta_degrees,interp_S_2D,15)
        
        
        % COMPARE the parameter values
        [Hs_derived_interp,Tm_derived_interp,direction_derived_interp,~] = ...
            waveLibrary().calculateSpectrumCharacteristics(interp_theta_degrees.*(pi/180),interp_f,interp_S_2D);
        
        % Display Test results
        disp("----------------------------------------------")
        disp(" Testing Interpolation of ERA5 data ")
        disp("----------------------------------------------")
        disp("Param | Original spectrum | Interp       ")
        disp("----------------------------------------------")
        disp([' Hs   |      ', num2str(Hs_derived), '       |     ', num2str(Hs_derived_interp)])
        disp([' Tm   |      ', num2str(Tm_derived),'       |      ', num2str(Tm_derived_interp)]);
        disp([' Dir_m |     ', num2str(mean_direction_derived_degrees), '       |     ', num2str(direction_derived_interp)])
        
        figure('Position', [0, 0, 800, 300]);
            subplot(1,2,1)
            contour(interp_f,interp_theta_degrees,interp_S_2D)
            colorbar;
            hold on; contour(f,theta_degrees,S_2D);
            subplot(1,2,2)
            plotLibrary().waveSpectrum2D(1, interp_S_2D, interp_f, interp_theta_degrees, 0.5, "")
            colorbar;
            hold on; plotLibrary().waveSpectrum2D(1, S_2D, f, theta_degrees, 0.5, "")
    
end

function testWaveNumberSpectrumWithInterpolation(f,theta_degrees,S_2D,interp_size)
        w = f.*(2*pi);
        k = (w).^2 ./ 9.81;
        
        [Hs_derived,Tm_derived,mean_direction_derived_degrees,~] = ...
            waveLibrary().calculateSpectrumCharacteristics(theta_degrees.*(pi/180),f,S_2D);

        [E_kx_ky, Jacobian_kx_ky, kx_matrix, ky_matrix] = waveLibrary().waveNumberSpectrum(S_2D,w,k,theta_degrees.*(pi/180));
        [wave_spectrum,Jacobian_f_theta] = waveLibrary().waveSpectrum(E_kx_ky, w, k);
        

        % INTERPOLATE
        [interp_S_2D, interp_f, interp_theta_degrees] = waveLibrary().interpolateWaveSpectrum(interp_size, f,theta_degrees,S_2D,"ERA5");
        [Hs_derived_interp,Tm_derived_interp,direction_derived_interp,~] = ...
            waveLibrary().calculateSpectrumCharacteristics(interp_theta_degrees.*(pi/180),interp_f,interp_S_2D);

        interp_w = interp_f.*(2*pi);
        interp_k = (interp_w).^2 ./ 9.81;
        [interp_E_kx_ky, Jacobian_kx_ky, interp_kx_matrix,interp_ky_matrix] = waveLibrary().waveNumberSpectrum(interp_S_2D,interp_w,interp_k,interp_theta_degrees.*(pi/180));
        [interp_wave_spectrum,Jacobian_f_theta] = waveLibrary().waveSpectrum(interp_E_kx_ky, interp_w, interp_k);
    
        figure('Position', [0, 0, 1200, 1200]);
            subplot(3,3,1); plotLibrary().waveSpectrum2D(1, S_2D,f, ...
                theta_degrees, 0.3, "Original ERA5 2D spectrum E(f,theta)");
            subplot(3,3,2); plotLibrary().waveNumberSpectrum(E_kx_ky, ...
                kx_matrix, ky_matrix, "Original Wave number spectrum, E(kx,ky)");
                xlim([-0.15 0.15]);ylim([-0.15 0.15]);
                grid on;
            subplot(3,3,3); 
                plotLibrary().waveSpectrum2D(1, wave_spectrum, f, ...
                    theta_degrees, 0.3, "Original Wave Spectrum after conversion " + ...
                    "E(Kx,ky) --> E(f,theta)");
            
            
            subplot(3,3,4); plotLibrary().waveSpectrum2D(1, interp_S_2D,interp_f, ...
                interp_theta_degrees, 0.3, "Interpolated 2D wave spectrum E(f,theta)");
            subplot(3,3,5); plotLibrary().waveNumberSpectrum(interp_E_kx_ky, ...
                interp_kx_matrix, interp_ky_matrix, "Interpolated Wave number spectrum, E(kx,ky)");
                xlim([-0.15 0.15]);ylim([-0.15 0.15]);
                grid on;
            subplot(3,3,6); 
                plotLibrary().waveSpectrum2D(1, interp_wave_spectrum, interp_f, ...
                    interp_theta_degrees, 0.3, "Interpolated Wave Spectrum after conversion " + ...
                    "E(Kx,ky) --> E(f,theta)");
                
            subplot(3,3,7); 
                plotLibrary().waveSpectrum2D(1, interp_S_2D,interp_f, ...
                interp_theta_degrees, 0.3, "Overlay E(f,theta)");
                hold on;
                plotLibrary().waveSpectrum2D(1, S_2D,f, ...
                theta_degrees, 0.3, "Overlay E(f,theta)");
            subplot(3,3,8); 
                plotLibrary().waveNumberSpectrum(interp_E_kx_ky, ...
                interp_kx_matrix, interp_ky_matrix, "Overlay E(Kx,ky)");
                xlim([-0.15 0.15]);ylim([-0.15 0.15]);grid on;
                hold on;
                plotLibrary().waveNumberSpectrum(E_kx_ky, ...
                kx_matrix, ky_matrix, "Overlay E(Kx,ky)");
                xlim([-0.15 0.15]);ylim([-0.15 0.15]);grid on;
            subplot(3,3,9); 
                plotLibrary().waveSpectrum2D(1, interp_wave_spectrum, interp_f, ...
                    interp_theta_degrees, 0.3, "Overlay E(Kx,ky) --> E(f,theta)");
                hold on;
                plotLibrary().waveSpectrum2D(1, wave_spectrum, f, ...
                    theta_degrees, 0.3, "Overlay E(Kx,ky) --> E(f,theta)");
    
        % COMPARE the parameter values
        [Hs_derived_interp_Ek,Tm_derived_interp_Ek,direction_derived_interp_Ek,~] = ...
            waveLibrary().calculateSpectrumCharacteristics(interp_theta_degrees.*(pi/180),interp_f,interp_wave_spectrum);
    
        % Display Test results
        disp("------------------------------------------------------------------")
        disp("              Testing Interpolation of ERA5 data                  ")
        disp("------------------------------------------------------------------")
        disp("Param | Original spectrum |     Interp    |   After Ek->Ef,theta  ")
        disp("------------------------------------------------------------------")
        disp([' Hs   |      ', num2str(Hs_derived), '       |     ', num2str(Hs_derived_interp), '    |     ', num2str(Hs_derived_interp_Ek)])
        disp([' Tm   |      ', num2str(Tm_derived),'       |      ', num2str(Tm_derived_interp), '   |     ', num2str(Tm_derived_interp_Ek)]);
        disp([' Dir_m |     ', num2str(mean_direction_derived_degrees), '       |     ', num2str(direction_derived_interp), '   |     ', num2str(direction_derived_interp_Ek)])
        
end

function testMultipleKTheory(era5)
    %%TESTMULTIPLE want to run multiple era5 spectra through the 'k theory'
    %%test to see what the errors are each time to see if there is any
    %%consistency.


    % Get attributes
    era5_longitudes = era5.getLongitude;
    era5_latitudes = era5.getLatitude;
    era5_time = era5.getTime;
    % era5_d2fd_all_5dims = era5.getAllWaveSpectrumD2FD;

    f = era5.FrequencyBins;

    theta_degrees = era5.DirectionBins;
   
    for t =1:length(era5_time) 
        chosen_time = era5_time(t);
        sar_transect_center_latitude = -34.204; %-32.204 has nice spectrum for Report!!!
        sar_transect_center_longitude = 18.28666944;
        [era5_d2fd, era5_lat, era5_lon, era5_datetime] = era5.getSlicedWaveSpectrumD2Fd(sar_transect_center_latitude,sar_transect_center_longitude,chosen_time);

        S_2D = era5_d2fd;
        S_2D_all(t,:,:) = era5_d2fd;
        interp_size = 128;

        % Interpolate the spectrum
        [interp_S_2D, interp_f, interp_theta_degrees] = waveLibrary().interpolateWaveSpectrum(interp_size, f,theta_degrees,S_2D,"ERA5");
        interp_S(t,:,:) = interp_S_2D;
        interp_f_all(t,:) = interp_f;
        interp_theta_degrees_all(t,:) =interp_theta_degrees;

        [Hs_original(t),Tm_original(t),mean_direction_original(t),Hs_test(t), Tm_test(t), mean_direction_test(t)] = testKTheory(0,interp_f, interp_theta_degrees, interp_S_2D, 15);
    
    end

    figure('Position', [0, 0, 1200, 300]);
            subplot(1,3,1); plot(Hs_original); hold on; plot(Hs_test);
            legend("Original","Test"); title("Hs");
            subplot(1,3,2); plot(Tm_original); hold on; plot(Tm_test);
            legend("Original","Test"); title("Tm");
            subplot(1,3,3); plot(mean_direction_original); hold on; plot(mean_direction_test);
            legend("Original","Test"); title("Mean Wave Direction");
end