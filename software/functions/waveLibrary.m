function waves = waveLibrary
% Wave Function Library 
%   This contains all of the functions relating to waves

waves.calc1DWaveSpec = @calc1DWaveSpec;
waves.waveNumberSpectrum = @waveNumberSpectrum; % E(f, direction) -> E(kx,ky)
waves.waveSpectrum = @waveSpectrum; % E(kx,ky) -> E(f, direction)
waves.lineariseKxandKyMatrices = @lineariseKxandKyMatrices;


waves.interpolateWaveNumberSpectrum = @interpolateWaveNumberSpectrum;
waves.calculate1DSpectrumCharacteristics = @calculate1DSpectrumCharacteristics;
waves.calculateSpectrumCharacteristics = @calculate2DSpectrumCharacteristics;

waves.interpolateWaveSpectrum = @interpolateWaveSpectrum;


waves.simulateWaveSpectrum = @simulateWaveSpectrum;
waves.generateDirectionalDistribution = @generateDirectionalDistribution;
end

function [S_1D, f] = calc1DWaveSpec(f, theta_rad, S_2D)
    S_1D = trapz(theta_rad,S_2D,1);
end

function [E_kx_ky, Jacobian, kx_matrix, ky_matrix] = waveNumberSpectrum(wave_spectrum, omega ,k , direction_bins_rad, mean_wave_dir, azimuth_angle)
%WaveNumberSpectrum convert a wave spectrum E(frequency, direction) to k domain, E(kx,ky)
%   According to Holthuijsen: "Remote-sensing and numerical wave models can 
%   estimate the full two-dimensional spectrum, usually the wave-number 
%   spectrum E (k x , k y )" [pg.52 Holthuijsen] 

% Wave Number Spectrum E(kx,ky)
    % First convert to omega from frequency
    era5_d2wd = wave_spectrum ./ (2 * pi); % conver E(f, theta) to E(w, theta)

    % Wave phase speed (c_w) for deep water
    % era5_cw_wave_speed = gravity ./ (era5_omega); % c_0, [Eq.5.4.24 Holthuijsen]
    era5_cw_wave_speed = omega ./ k; % c, [Eq.3.5.36 Holthuijsen] same as above
    
    % Group wave speed (c_g)
    era5_n_dispersion_scaling_factor = 0.5; % 0.5 is used for deep water, [Eq.5.4.32 Holthuijsen]
    era5_cg_group_wave_speed = era5_n_dispersion_scaling_factor * era5_cw_wave_speed; % [Eq.5.4.32 Holthuijsen]
    
    % Calculate the wave number spectrum E(kx,ky)
    Jacobian = ((era5_cw_wave_speed .* era5_cg_group_wave_speed) ./ omega ) ; % [Eq.3.5.36 Holthuijsen]
    E_kx_ky = Jacobian .* era5_d2wd; % [Eq.3.5.36 Holthuijsen]
    E_kx_ky(isnan(E_kx_ky)) = 0;
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
    % direction = direction_bins_rad(105);
   
    % kx_matrix = k .* sind(mean_wave_dir+abs(azimuth_angle)); % [rad/m] [Eq.3.5.19b Holthuijsen]
    % ky_matrix = k .* cosd(mean_wave_dir+abs(azimuth_angle)); % [rad/m] [Eq.3.5.19b Holthuijsen]
    kx_matrix = k .* sin(abs(direction_bins_rad)); % [rad/m] [Eq.3.5.19b Holthuijsen]
    ky_matrix = k .* cos(abs(direction_bins_rad)); % [rad/m] [Eq.3.5.19b Holthuijsen]
    
    
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

function [wave_spectrum, Jacobian] = waveSpectrum(wave_number_spectrum, omega, k)
%WaveSpectrum convert E(kx,ky) -> E(frequency, theta)

    % Calculate wave phase speed 
    cw = omega./k;
    % Calculate group wave phase speed 
    cg = 0.5 .* cw;
    
    % Calculate the E(w,theta) spectrum
    Jacobian = omega ./ (cw .* cg);
    d2wd = wave_number_spectrum .* Jacobian;
    % Calculate the E(f, theta) spectrum
    wave_spectrum = d2wd .* (2 * pi);
    
    % f = omega ./ (2*pi);
    % theta_degrees = asind(kx_matrix./k);
        
    % kx_matrix = k .* sind(direction_bins_degrees); % [rad/m] [Eq.3.5.19b Holthuijsen]
    % ky_matrix = k .* cosd(direction_bins_degrees); % [rad/m] [Eq.3.5.19b Holthuijsen]

    % figure;
    % contour(omega./(2*pi),direction_bins_degrees, wave_spectrum');
    % xlabel("Frequency"); ylabel('Direction [degrees]');title("Final 2-D Wave Spectrum, E(f,theta)", "calculated from E(k)"); c = colorbar();c.Label.String = '[m^2 s / degree]"';
end

function [linspace_kx_matrix,linspace_ky_matrix,regridded_E_f_theta] = lineariseKxandKyMatrices(kx_matrix,ky_matrix,E_f_theta)

    matrix_size = size(kx_matrix,1); % Should be square

    linspace_kx_matrix = linspace(min(abs(kx_matrix(:))),max(abs(kx_matrix(:))),matrix_size/2);
    linspace_kx_matrix = [-1*fliplr(linspace_kx_matrix) linspace_kx_matrix]; % Necessary to remove k=0 as this does not exist.
    
    matrix_size = size(ky_matrix,1); % Should be square
    linspace_ky_matrix = linspace(min(abs(ky_matrix(:))),max(abs(ky_matrix(:))),matrix_size/2);
    linspace_ky_matrix = [-1*fliplr(linspace_ky_matrix) linspace_ky_matrix]';% Necessary to remove k=0 as this does not exist.
    
    % linspace_ky_matrix = ones(matrix_size,matrix_size).*linspace_ky_matrix;

    [linspace_kx_matrix, linspace_ky_matrix] = meshgrid(linspace_kx_matrix, linspace_ky_matrix);

    % Perform interpolation
    regridded_E_f_theta = griddata(kx_matrix, ky_matrix, E_f_theta, linspace_kx_matrix, linspace_ky_matrix, 'linear');
    regridded_E_f_theta(isnan(regridded_E_f_theta))=0;
    % regridded_E_f_theta = E_f_theta;

end

function [Hs_m0,Tm] = calculate1DSpectrumCharacteristics(f,S_1D)
    %CALCULATE1DSPECTRUMCHARACTERICS 
    
    m0 = trapz(f,S_1D,2);
    Hs_m0 = 4 * sqrt(m0);
    m1 = trapz(f,S_1D .* f.^1,2);
    Tm = (m1/m0)^(-1); % [Holthuijsen Eq.5] 
end

function [Hs_m0,Tm,mean_wave_direction_degrees,total_variance_or_energy] = calculate2DSpectrumCharacteristics(direction_bins_rad,frequency_bins,E_f_theta)
    %%calculate2DSpectrumCharacteristics Calculate the integral values of the wave
    %%spectrum
        % frequency_bins         = original frequency array; [1xN row vector]
        % direction_bins         = original direction array in RADIANS; [Mx1 column vector]
        % E_f_theta              = original 2D wave spectrum on the f, theta grid;
        %                           [MxN matrix]

% Significant wave height

    % Calculate the 0th order moments (m_n) [Eq.8]
    m0 = trapz(direction_bins_rad, trapz(frequency_bins ,frequency_bins.^0 .* E_f_theta,2),1); %[Eq.8], dimension 2,3 = direc,freq

    % Calculate the significant wave height (Hs) [Eq.4]
    Hs_m0 = 4*sqrt(m0);

    % Calculate mean wave height (H)
    % H = sqrt(pi/8) * Hs_m0; %[Holthuijsen Eq.4.2.6]
    
% Mean Wave period
    % Calculate the -1st order moments (m_n)
    intermediate = trapz(direction_bins_rad,(1./frequency_bins) .* E_f_theta,1);
    intermediate(isnan(intermediate))=0;
    m_neg1 = trapz(frequency_bins,intermediate,2); %[Eq.8]

    intermediate = trapz(direction_bins_rad,(frequency_bins) .* E_f_theta,1);
    intermediate(isnan(intermediate))=0;
    m1 = trapz(frequency_bins,intermediate,2); %[Eq.8]

    % Calculate the mean wave period (Tm) [Eq.5]
    % Tm = m_neg1/m0; % [Eq.5]
    Tm = (m1/m0)^(-1); % [Holthuijsen Eq.5]
    Tz = (m0/m_neg1)^(1/2); %Tz --Average zero-crossing period [s]
    % T0 = 1./w0; %T0  --Modal (peak) period (T0 = 2 pi /w0) [s]
    
% Mean Wave direction
    % Calculate the angular moments [Eq.9&10]
    a1 = trapz(frequency_bins,trapz(direction_bins_rad,cos(direction_bins_rad).*E_f_theta,1),2); %[Eq.9]
    b1 = trapz(frequency_bins,trapz(direction_bins_rad,sin(direction_bins_rad).*E_f_theta,1),2); %[Eq.10]
   
    % Calculate the mean Direction of the wave
    mean_wave_direction = atan(b1./a1); %[Eq.6]
    mean_wave_direction_degrees = mean_wave_direction *180/pi;

    if mean_wave_direction_degrees < 0
        mean_wave_direction_degrees =  180 + mean_wave_direction_degrees; % To get the direction as "coming from"
    else
        mean_wave_direction_degrees = mean_wave_direction_degrees; % To get the direction as "coming from"
    end
    
   
% Spectral Denisty
    % Total variance/energy of the waves = spectral density integrated over all
    % frequencies and wave numbers at specified time of interest
    total_variance_or_energy = trapz(trapz(E_f_theta,1),2);

end

%%
function [interp_E_f_theta, interp_f, interp_theta_degrees] = interpolateWaveSpectrum(interp_size, f,theta_degrees,E_f_theta,ModelName)
    %%INTERPOLATEWAVESPECTRUM Interpolates the original 2D wave spectrum to
    %%the new desired size
        % interp_size   = size of new interpolation grid
        % f             = original frequency array; [1xN row vector]
        % theta         = original direction array; [Mx1 column vector]
        % E_f_theta     = original 2D wave spectrum on the f, theta grid;
        %                 [MxN matrix
        % ModelName     = name of the source of the data i.e. ERA5

    if ModelName == "ERA5"
       f_min = 0.03453; % AS per the ERA5 documentation
       % Need to preserve the statistical property regardless of the
       % size of the array. i.e. it still needs to mirror the case as
       % if there were only 30 bins (as we expect in ERA5 data).
       f_max = max(f); % Calculated from the 30 bin case

       % Setup the frequency array as per the ERA5 documentation.
        interp_f = zeros(1,interp_size);
        interp_f(1) = f_min;
        a = (f_max / f_min)^(1 / (interp_size - 1));
    
        for i=2:length(interp_f)
            interp_f(i) = interp_f(i-1)*a;
        end


        max_degree_bin = max(theta_degrees);%352.5; % ERA5 maximum from docs
        min_degree_bin = min(theta_degrees);%7.5; % ERA5 min from docs
        
        interp_theta_degrees = linspace(min_degree_bin,max_degree_bin,interp_size)';

        interp_E_f_theta = griddata(f,theta_degrees, E_f_theta, interp_f, interp_theta_degrees,'cubic');
        interp_E_f_theta(isnan(interp_E_f_theta)) = 0;
        interp_E_f_theta(interp_E_f_theta<0) = 0;

    end
end
% function [interp_E_kx_ky,interp_kx_matrix, interp_ky_matrix,interp_k_matrix] = interpolateWaveNumberSpectrum(sar_transect_size,sar_range_resolution, sar_azimuth_resolution, wave_number_spectrum, kx_matrix, ky_matrix)
% %interpolateWaveNumberSpectrum create a 128x128 wave number spectrum
% %   The original spectrum is smaller than the SAR image; interpolation is
% %   used to increase the size.
% 
%     % % Resize the kx and ky matrices.
%     % extrap_kx = linspace((min(min(abs(era5_kx_matrix)))),(max(max(abs(era5_kx_matrix)))),sar_sub_transect_size/2);
%     % extrap_kx = [-1*fliplr(extrap_kx) extrap_kx];
%     % interp_kx = linspace(min(kx_matrix(:)),max(kx_matrix(:)),sar_sub_transect_size/2);
% 
%     % NEW OPTION
%     % sar_dx_range=5;
%     % sar_dy_azimuth=5;
%     % kx_grid = linspace(-sar_dx_range/pi, sar_dx_range/pi, 128);
%     % ky_grid = linspace(-sar_dy_azimuth/pi, sar_dy_azimuth/pi, 128);
% 
% 
%     min_kx = min(min(abs(kx_matrix)));
%     max_kx = max(max(abs(kx_matrix)));
%     interp_kx = [-1*linspace(max_kx,min_kx,sar_transect_size/2) linspace(min_kx,max_kx,sar_transect_size/2)];
% 
%     min_ky = min(min(abs(ky_matrix)));
%     max_ky = max(max(abs(ky_matrix)));
%     interp_ky = [-1*linspace(max_ky,min_ky,sar_transect_size/2) linspace(min_ky,max_ky,sar_transect_size/2)];
% 
%     % OLD WAY
%     % interp_kx = (2*pi/(sar_range_resolution*sar_transect_size)) * ([-1*sar_transect_size/2:1:-1 1:1:sar_transect_size/2]);
%     % interp_ky = (2*pi/(sar_azimuth_resolution*sar_transect_size)) * ([-1*sar_transect_size/2:1:-1 1:1:sar_transect_size/2]);
%     % extrap_ky = linspace((min(min(abs(era5_ky_matrix)))),(max(max(abs(era5_ky_matrix)))),sar_sub_transect_size/2);
%     % extrap_ky = [-1*fliplr(extrap_ky) extrap_ky];
% 
%     %% GIACOMO's WAY
%     % Flat interpolation Method
%     kx_flat = kx_matrix(:);
%     ky_flat = ky_matrix(:);
%     S_flat = wave_number_spectrum(:);
% 
%     F = scatteredInterpolant(kx_flat , ky_flat , S_flat, 'natural','none');
%     [interp_kx_matrix, interp_ky_matrix] = meshgrid(interp_kx,interp_ky);
%     interp_E_kx_ky = F(interp_kx_matrix,interp_ky_matrix);
%     interp_E_kx_ky(isnan(interp_E_kx_ky))=0;
% 
%     interp_k_matrix = sqrt(interp_ky_matrix.^2 + interp_kx_matrix.^2);
% 
%     % % Non-flat Interpolation
%     % [KX, KY] = meshgrid(kx, ky);
%     % [KX_new, KY_new] = meshgrid(new_kx, new_ky);
%     % S_k_interpolated = interp2(KX, KY, wave_number_spectrum, KX_new, KY_new, 'linear');
% 
% 
%     % [interp_k_matrix, indx] = sort(interp_k_matrix);
% 
%     % The code below is a sanity check to understand the effects of the interpolation on the spectrum. The peak wave period of the original (old) and interpolated (new) are calculated and compared.
%     % % Check the quality of the extrapolation (internal sanity check)
% 
%     % [test_index1_new, test_index2_new] = find(interp_E_kx_ky==(max(max(interp_E_kx_ky))));
%     % test_k_value_new = sqrt(interp_kx_matrix(test_index1_new, test_index2_new).^2 + interp_ky_matrix(test_index1_new, test_index2_new).^2);
%     % test_omega_value_new = sqrt(gravity * test_k_value_new);
%     % test_period_new = 2*pi/test_omega_value_new
%     % [test_index1_old, test_index2_old] = find(era5_wave_number_spectrum==(max(max(era5_wave_number_spectrum))));
%     % test_k_value_old = sqrt(era5_kx_matrix(test_index1_old, test_index2_old).^2 + era5_ky_matrix(test_index1_old, test_index2_old).^2);
%     % test_omega_value_old = sqrt(gravity .* test_k_value_old);
%     % test_period_old = 2*pi/test_omega_value_old
%     % % TODO: COMPARE TO THE PEAK WAVE PERIOD
% 
% 
%     % Marcello
% 
% end


% function [interp_E_kx_ky,interp_kx_matrix, interp_ky_matrix] = interpolateWaveNumberSpectrum(sar_transect_size,sar_range_resolution, sar_azimuth_resolution, wave_number_spectrum, kx_matrix, ky_matrix)
% %interpolateWaveNumberSpectrum create a 128x128 wave number spectrum
% %   The original spectrum is smaller than the SAR image; interpolation is
% %   used to increase the size.
% 
%     % % Resize the kx and ky matrices.
%     % extrap_kx = linspace((min(min(abs(era5_kx_matrix)))),(max(max(abs(era5_kx_matrix)))),sar_sub_transect_size/2);
%     % extrap_kx = [-1*fliplr(extrap_kx) extrap_kx];
%     interp_kx = (2*pi/(sar_range_resolution*sar_transect_size)) * ([-1*sar_transect_size/2:1:-1 1:1:sar_transect_size/2]);
% 
%     interp_ky = (2*pi/(sar_azimuth_resolution*sar_transect_size)) * ([-1*sar_transect_size/2:1:-1 1:1:sar_transect_size/2]);
%     % extrap_ky = linspace((min(min(abs(era5_ky_matrix)))),(max(max(abs(era5_ky_matrix)))),sar_sub_transect_size/2);
%     % extrap_ky = [-1*fliplr(extrap_ky) extrap_ky];
% 
%     [interp_kx_matrix, interp_ky_matrix] = meshgrid(interp_kx,interp_ky);
%     F = scatteredInterpolant(kx_matrix(:) , ky_matrix(:) , wave_number_spectrum(:), 'natural');
%     interp_E_kx_ky = F(interp_kx_matrix,interp_ky_matrix);
% 
%     % The code below is a sanity check to understand the effects of the interpolation on the spectrum. The peak wave period of the original (old) and interpolated (new) are calculated and compared.
%     % % Check the quality of the extrapolation (internal sanity check)
% 
%     % [test_index1_new, test_index2_new] = find(interp_E_kx_ky==(max(max(interp_E_kx_ky))));
%     % test_k_value_new = sqrt(interp_kx_matrix(test_index1_new, test_index2_new).^2 + interp_ky_matrix(test_index1_new, test_index2_new).^2);
%     % test_omega_value_new = sqrt(gravity * test_k_value_new);
%     % test_period_new = 2*pi/test_omega_value_new
%     % [test_index1_old, test_index2_old] = find(era5_wave_number_spectrum==(max(max(era5_wave_number_spectrum))));
%     % test_k_value_old = sqrt(era5_kx_matrix(test_index1_old, test_index2_old).^2 + era5_ky_matrix(test_index1_old, test_index2_old).^2);
%     % test_omega_value_old = sqrt(gravity .* test_k_value_old);
%     % test_period_old = 2*pi/test_omega_value_old
%     % % TODO: COMPARE TO THE PEAK WAVE PERIOD
% end

%%
function [S_1D_f, S_2D_cosm] = simulateWaveSpectrum(Hs,gamma,f, T0,  mean_wave_direction_degrees, direction_bins_degrees, PlotFlag,spectrumType)
    %SIMULATEWAVESPECTRUM generates a simulated 2D wave spectrum E(f,theta)
        % 

        % Use:  S = simulateWaveSpectrum(Hs,Tm,gamma,w0)
        %
        % Inputs:
        %   Hs = significant wave height 
        %   Tm = mean period
        %   gamma = Peakedness factor; usually between 1<gamma<5
        %   w = Column vector of wave frequencies [rad/s]
        %   w0 = modal (peak) frequency [rad/s]
        %   Tsig = significant wave period [m/s]
        %   PlotFlag	- 1 to plot the spectrum, 0 for no plot
        %
        % Output:
        %   S = matrix vector of wave spectrum values [m^2 s]; evaluated at W(k) 

        % Constants
        g = 9.81; % Gravitational acceleration [m/s²]
        fp = 1 / T0; % Peak frequency [Hz]
        f0 = fp / 0.877; % PM scaling correction
        
        % Compute 1D wave spectrum
        if strcmp(spectrumType, 'PM')
            % Pierson-Moskowitz Spectrum
            alpha = 8.1e-3;
            beta = 0.74;
            S_1D_f = (alpha * g^2) ./ ((2 * pi * f).^5) .* exp(-beta * (f0^4) ./ (f.^4));
        elseif strcmp(spectrumType, 'JONSWAP')
            % JONSWAP Spectrum
            alpha = 5e-3 * Hs^2 * T0^(-4); % Scale parameter
            sigma = zeros(size(f));
            sigma(f <= fp) = 0.07;
            sigma(f > fp) = 0.09;
            G = exp(-((f - fp).^2) ./ (2 * sigma.^2 * fp^2)); % Peakedness factor
            S_1D_f = alpha * g^2 .* f.^(-5) .* exp(-5/4 * (fp ./ f).^4) .* gamma.^G;
        else
            error('Invalid spectrum type. Use ''PM'' or ''JONSWAP''.');
        end
        
        % Compute total variance (m0) from 1D spectrum
        target_m0 = (Hs / 4)^2;
        computed_m0 = trapz(f, S_1D_f);
        
        % **Scaling fix**: Ensure the computed variance matches the target variance
        scaling_factor = target_m0 / computed_m0;
        S_1D_f = S_1D_f * scaling_factor;

        S_1D_omega = S_1D_f ./ (2*pi);
        
        % CREATE 2D MODEL
        % direction_bins = direction_bins_degrees .* pi/180;
        % mean_wave_direction = mean_wave_direction_degrees * pi/180;
        % [D_cosm] = generateDirectionalDistribution(2, mean_wave_direction,direction_bins);
        % %Test if D_cosm is in rads or degrees
        % test = trapz(direction_bins,D_cosm);
        % S_2D_cosm = S_1D_f .* D_cosm;

        % Convert direction to radians
        theta_r = direction_bins_degrees .* (pi / 180);
        theta_0 = mean_wave_direction_degrees * pi / 180;

        % Compute directional spreading function using a cosine power model
        D_cosm = cos((theta_r - theta_0) ./ 2).^2;
        D_cosm = D_cosm ./ trapz(theta_r, D_cosm); % Normalize

        % Compute 2D spectrum
        S_2D_cosm = S_1D_f .* D_cosm;

end


function [D_cosm] = generateDirectionalDistribution(width_of_distribution, mean_wave_direction, direction_bins)
    %GENERATEDIRECTIONALDISTRIBUTION based on Holthuijsen pg.165
    %
    % s -- controls the width of the distribution and is a value between 1
    %       and 16 (Holthuijsen)
    % D(meandirection) = D(meandirection,f) -- [1/angle = 1/degree or
    %       1/rad]
    % m -- width of the distribution; default = 2; [Holthuijesen Eq.6.3.23]

    % D(meandirection;f) gives the normalised distribution of the wave
    % energy density over directions at ONE frequency, whereas E(f,theta)
    % gives the non-normalised distribution over both frequency and
    % direction.
    
    % direction_bins = direction_bins_degrees.* pi/180;
    

    % Cos^m(meandirection) model - limits direction to -90<meandirection<90
    % Limits the wave propogation to a sector of 90 degrees on either side
    % of the mean wave direction.
    %   Most widely used version is cos^2(meandirection) model by [Pierson,
    %   et al. 1952] % default = 2 [Holthuijesen Eq.6.3.23]
    %   But generalised to cos^m(meandirection) model

    m = width_of_distribution; %[Holthuijesen Eq.6.3.24]
    A1 = gamma((0.5 * m) +1)/(gamma((0.5*m) +0.5)*sqrt(pi)); % Normalisation coefficient [Holthuijesen Eq.6.3.25]
    for i = 1: length(direction_bins)
        if abs(direction_bins(i) - mean_wave_direction) <= pi/2
            D_cosm(i) = A1 * cos(direction_bins(i)-mean_wave_direction)^m; % [1/rads] [Holthuijesen Eq.6.3.24]
        else
            D_cosm(i) = 0;
        end
    end

    % % sigma_theta = Directional spreading of the waves [Holthuijesen Eq.6.3.22]
    % for i=1:length(freqArray)
    %     if freqArray(i) < f_peak
    %         sigma_theta(i) = 26.9 .* (freqArray(i)/f_peak).^(-1.05);
    %     else
    %         sigma_theta(i) = 26.9 .* (freqArray(i)/f_peak).^(0.68);
    %     end
    % end
    % 
    % s = (2./ (sigma_theta.*pi./180).^2) -1; % [Holthuijesen Eq.6.3.26] CHECK THIS!!!
    % 
    % % Cos^2s(0.5 meandirection) model 
    % % Does not have the same limit as Cos^m(meandirection) model, instead allows
    % % waves to propagate against the wind (i.e. meandirection larger than 90
    % % degrees).
    % A2 = gamma(s+1)./(gamma(s+0.5).*(2*sqrt(pi))); % [Holthuijesen Eq.6.3.25]
    % % For -180<theta<=180 
    % for i = 1: length(s)
    %     for j = 1: length(direction_bins)
    %         direction_difference = abs(direction_bins(j) - mean_wave_direction);
    %         if direction_difference <= 180 && direction_difference>-180
    %             D_cos2s(i,j,:) = abs(A2 .* cosd(0.5 .* direction_bins(j)).^(2*s(i))); % [Holthuijesen Eq.6.3.25]
    %         else
    %             D_cos2s(i,j,:) = 0;
    %         end
    %     end
    % end


    % WRAPPED NORMAL DISTRIBUTION

end