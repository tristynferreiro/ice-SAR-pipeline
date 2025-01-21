function handh = hAndH1991Functions
% H&H,1991: All functions used in H&H, 1991
%{

ALL ANGLES ARE IN DEGREES, all trig functions are in degrees
ALL k_neg functions are the function when -k is used.

Variable Definitions:
    zeta_k = wave amplitude of the wave spectrum (the fourier coefficient)

    TR_k = RAR MTF
    TR_k_neg = TR(-k) = additive inverse of TR_k
    Tt_k = tilt MTF
    Th_k = hydrodynamic MTF
    k_l = component of the incident wave number vector in the radar look
            direction 
    mu = damping factor, introduced to describe the response of the short
            waves to the long wave modulation, μ
    Y_r+ 1i*Y_i = complex feedback factor, representing the long wave
            modulation of the wind input to the short waves
    k = wave number vector
    omega = gravity wave frequency, ω=√gk
    F_k = F(k) OR E(k) = ocean wave variance spectrum
    F_k_neg = F(-k) = additive negerse of F_k

    xi = azimuthal displacement of the backscattering element, ξ
    beta = β, velocity bunching parameter
    v = satellite_orbital_velocity OR upsilon, υ OR 
                               = range component of the long wave orbital 
                                 velocity 
                               = time average over the period during which 
                                 the scattering element is viewed by the SAR.
                                 BUT in H&H, first order is equal to the 
                                 INSTANTANEOUS orbital velocity in the 
                                 center of the viewing window 
    sar_platform_velocity = U
    sar_slant_range = rho, ρ OR R OR slant range from satellite
    Tv_k = Tυ_k = range velocity (υ) transfer function
    sar_incidence_angle_degrees = theta, θ in degrees
    TS_k = net SAR imaging MTF
    
    fv_r = f^v(r) = orbital velocity covariance function
    r = the x and y coords vector

%}

handh.klUsingSARlook = @klUsingSARlook;
% The Frozen Surface Contribution / RAR
handh.gravityWaveFrequency = @defineGravityWaveFrequency;
handh.rarMTF = @rarMTF;
handh.tiltMTF = @tiltMTF;
handh.hydrodynamicMTF = @hydrodynamicMTF;
handh.hydrodynamicMTFwithComplexFeedback = @hydrodynamicMTFwithComplexFeedback;
handh.rarImageVarianceSpectrum = @rarImageVarianceSpectrum;

% Motion Effects of SAR
handh.azimuthalDisplacement = @azimuthalDisplacement;
handh.defineBeta = @defineBeta;
handh.rangeVelocityMTF = @rangeVelocityMTF;
handh.sarImageAmplitudeSpectrum = @sarImageAmplitudeSpectrum;
handh.velocityBunchingMTF = @velocityBunchingMTF;
handh.sarImagingMTF = @sarImagingMTF;
handh.sarImageVarianceSpectrumLinearMappingTransform = @sarImageVarianceSpectrumLinearMappingTransform; %SAR image variance spectrum

% General Nonlinear Mapping Relation
handh.orbitalVelocityCovarianceFunction = @orbitalVelocityCovarianceFunction;
handh.meanSquareAzimuthalDisplacement = @meanSquareAzimuthalDisplacement;
handh.sarImageVarianceSpectrumQuasilinearMappingTransform = @sarImageVarianceSpectrumQuasilinearMappingTransform;

handh.rarImageIntensityAutocovariance = @rarImageIntensityAutocovariance;
handh.rarImageIntensityCovariance = @rarImageIntensityCovariance;
handh.spectralExpansion2nTerm = @spectralExpansion2nTerm;
handh.spectralExpansion2nMinus1Term = @spectralExpansion2nMinus1Term;
handh.spectralExpansion2nMinus2Term = @spectralExpansion2nMinus2Term;
handh.sarImageVarianceSpectrumNonlinearMappingTransform = @sarImageVarianceSpectrumNonlinearMappingTransform;

% Generate SAR Spectrum
handh.generateSARSpectrumFromWaveNumberSpectrum = @generateSARSpectrumFromWaveNumberSpectrum;
end

% From paper
function omega = defineGravityWaveFrequency(g,k)
    omega = sqrt(g*k);
end 

function [k_l] = klUsingSARlook(sar_look_direction,k_range)
    % This is defined on the bottom of page 10,715 after [Eq.6, H&H 1991] 
    if(isequal(sar_look_direction,"right"))
        k_l = -1*k_range;
        % k_y = -1*k_y;
    elseif(isequal(sar_look_direction,"left"))s
        k_l = k_range;
        % k_y = k_y;
    end
end

% -------------------------------------------------------------------
%% SAR Imaging of Ocean Waves: Frozen Surface Contribution (pg.10,715 - 10,716)

function TR_k = rarMTF(Tt_k, Th_k)
    %[Eq.4, H&H 1991]
    TR_k = Tt_k + Th_k;
end

function Tt_k = tiltMTF(sar_polarisation, sar_incidence_angle_degrees, k_l)
    %[Eq.5, H&H 1991]
    if(strcmp(sar_polarisation, 'VV') && all(sar_incidence_angle_degrees <= 60, 'all'))
        Tt_k = 4i .* k_l .* cotd(sar_incidence_angle_degrees) ./(1 + sind(sar_incidence_angle_degrees).^2); 
    elseif (strcmp(sar_polarisation, 'HH') && all(sar_incidence_angle_degrees <= 60, 'all'))
        Tt_k = 8i.* k_l ./ (sind(2 .* sar_incidence_angle_degrees)); 
    end
    Tt_k = abs(Tt_k);
end

function Th_k = hydrodynamicMTF(omega, mu, k, k_range,Y_r,Y_i)
    %[Eq.6, H&H 1991] 
    % We do not use the complex feedback factor: Y_r+1i*Y_i term, thus it
    % is generally set to 0
    
    complex_feedback_term = Y_r+1i*Y_i;
    Th_k = (omega - (1i * mu)) ./ (omega.^2 + mu.^2) .* 4.5 .* k .* omega .* ((k_range.^2./k.^2)+complex_feedback_term); %[Eq.6, H&H 1991] 
    Th_k = abs(Th_k);
end

function PR_k = rarImageVarianceSpectrum(TR_k, F_k, TR_k_neg, F_k_neg)
    %[Eq.13, H&H 1991]
    PR_k = 0.5 .* (abs(TR_k).^2.*F_k + abs(TR_k_neg).^2.*F_k_neg);
end

% function IR_k = imageModulationIntensity(TR_k,zeta_k,Tr_k_neg,zeta_k_neg)
%     %[Eq.10, H&H 1991]
%     IR_k = TR_k*zeta_k+conj(Tr_k_neg,zeta_k_neg);
% end

% -------------------------------------------------------------------
%% SAR Imaging of Ocean Waves: Motion Effects (pg.10,716 - 10,717)
function Tv_k = rangeVelocityMTF(omega, sar_incidence_angle_degrees, k_range, k)
    %[Eq.17, H&H 1991]
    Tv_k = -1 .* omega .* ((sind(sar_incidence_angle_degrees) .* (k_range ./ abs(k))) + (1i .* cosd(sar_incidence_angle_degrees))); 
    Tv_k = abs(Tv_k);
end

function beta = defineBeta(sar_slant_range, sar_platform_velocity)
    %[Eq.15, H&H 1991], velocity bunching parameter
    beta = sar_slant_range / sar_platform_velocity;
end

function Tvb_k = velocityBunchingMTF(sar_beta, k_azimuth, Tv_k)
    %[Eq.24, H&H 1991]
    Tvb_k = -1 .* sar_beta .* k_azimuth .* Tv_k;% [Eq.24, H&H 1991]

    Tvb_k = abs(Tvb_k);
end

function TS_k = sarImagingMTF(TR_k, Tvb_k)
    %[Eq.27, H&H 1991]
    TS_k = TR_k + Tvb_k;
end

function PS_k = sarImageVarianceSpectrumLinearMappingTransform(TS_k, F_k, TS_k_neg, F_k_neg)
    %[Eq.26, H&H 1991]
    PS_k = ((abs(TS_k).^2 .* F_k./2) + (abs(TS_k_neg).^2 .* F_k_neg./2));
end

% -------------------------------------------------------------------
%% General nonlinear mapping relation pg.5-6

function fv_r = orbitalVelocityCovarianceFunction(F_k, Tv_k)
    %[Eq.43, H&H 1991]
    fv_r = abs(ifftshift(ifft2(F_k .* abs(Tv_k).^2)));
end

function xi_sqr = meanSquareAzimuthalDisplacement(sar_beta, fv_r)
    %[Eq.44, H&H 1991]
    % The integral term in Eq.44 is equivalent to fv(0) in Eq.43.
    xi_sqr = (sar_beta).^2 .* abs(fv_r(1,1)); % [Eq.44, H&H 1991]
end

function [PS_ql_k,azimuthal_cutoff_factor] = sarImageVarianceSpectrumQuasilinearMappingTransform(k_azimuth, xi_sqr, PS_1)
    % [Eq.56, H&H 1991]
        % PS_1 = Linear Mapping Transform SAR Spectrum % [Eq.55, H&H 1991]
        % xi_sqr = mean square azimuthal displacement 
        % sar_transect_size is typically = 128

    % Apply filter to filter out high frequencies
    azimuthal_cutoff_factor = exp(-1 .* k_azimuth.^2 .* xi_sqr);
    PS_ql_k = azimuthal_cutoff_factor .* PS_1; % [Eq.56, H&H 1991]
end

function fR_r = rarImageIntensityAutocovariance(PR_k)
    %[Eq.47, H&H 1991]
    fR_r = 0.5 * abs( ifftshift( ifft2(PR_k)) );
end

function  fRv_r = rarImageIntensityCovariance(F_k, TR_k, Tv_k, F_k_neg, TR_k_neg, Tv_k_neg)
    %[Eq.48, H&H 1991]
    % Covariance between the RAR image intensity, I(x), and the orbital
    % velocity, v(x)
    fRv_r = 0.5 .* abs(ifftshift(ifft2(F_k .* TR_k .* conj(Tv_k) + F_k_neg .* conj(TR_k_neg) .* Tv_k_neg)));
end

function PS_2n = spectralExpansion2nTerm(n, fv_r)
    %[Eq.51, H&H 1991]
    PS_2n = (2*pi)^(-2) .* ifftshift(ifft2((fv_r.^n)./factorial(n)));
end

function PS_2n_minus_1 = spectralExpansion2nMinus1Term(n,fv_r, fRv_r, fRv_neg_r)
    %[Eq.52, H&H 1991]
    PS_2n_minus_1 = (2*pi)^(-2) .* ifftshift(ifft2( ...
        (1i .* (fRv_r - fRv_neg_r) .* fv_r.^(n-1)) ./ factorial(n-1) ...
        ));
end

function PS_2n_minus_2 = spectralExpansion2nMinus2Term(n, fv_r, fRv_r, fRv_neg_r, fRv_0, fR_r)
    %[Eq.53, H&H 1991]
    if (n == 1)
        % second_factorial_term = 0;
        PS_2n_minus_2 = (2*pi)^(-2) .* ifftshift(ifft2( ...
        (1/factorial(n-1)) .* fR_r .* fv_r.^(n-1) ...
        ));
    else
        second_factorial_term = 1/factorial(n-2);
        PS_2n_minus_2 = (2*pi)^(-2) .* ifftshift(ifft2( ...
        (1/factorial(n-1)) .* fR_r .* fv_r.^(n-1) + ...
        (second_factorial_term) .* (fRv_r - fRv_0) .* (fRv_neg_r - fRv_0) .* fR_r.^(n-2) ...
        ));
    end
    
end

function PS_nl = sarImageVarianceSpectrumNonlinearMappingTransform(plotsON, nonlinearity_order, PS_ql, k_azimuth, sar_beta, fv_r, fRv_r, fR_r, xi_sqr, k_range, sar_dazimuth, sar_sub_transect_size)
    % [Eq.50, H&H 1991]
        % PS_ql = Quasilinear Mapping Transform Spectrum, this is PS_1 (the
        % first term of the nonlinearity mapping transform)\
        % sar_transect_size is typically = 128
    ps_k = 0;       

    % Calculate the spectrum for all nonlinear orders > 1
    for nonlinearity = 2:nonlinearity_order
        % Spectral Expansion Term 1
        PS_2n = abs(spectralExpansion2nTerm(nonlinearity, fv_r)); % [Eq.51, H&H 1991]
        coefficient_2n = (k_azimuth .* sar_beta).^(2*nonlinearity);
    
        % Spectral Expansion Term 2
        PS_2n_minus_1 = abs(spectralExpansion2nMinus1Term(nonlinearity,fv_r, fRv_r, rot90(fRv_r,2))); % [Eq.52, H&H 1991]
        coefficient_2n_minus_1 = ((k_azimuth .* sar_beta).^(2*nonlinearity-1));
    
        % Spectral Expansion Term 3
        PS_2n_minus_2 = abs(spectralExpansion2nMinus2Term(nonlinearity, fv_r, fRv_r, rot90(fRv_r,2),  fRv_r(1,1), fR_r)); % [Eq.53, H&H 1991]
        coefficient_2n_minus_2 = (k_azimuth .* sar_beta).^(2*nonlinearity-2);
        
        sum = ((coefficient_2n .* PS_2n) + (coefficient_2n_minus_1 .* PS_2n_minus_1) + (coefficient_2n_minus_2 .* PS_2n_minus_2));
        
        % % Alternate approach for the code:
        % sum = zeros(sar_sub_transect_size);
        % for m = (2*nonlinearity-2):2*nonlinearity
        %     coefficient = (first_guess_ky_azimuth .* sar_beta).^m;
        %     if m == 2*nonlinearity-2
        %         PS_k = coefficient .* PS_2n_minus_2;
        %     elseif m == 2*nonlinearity-1
        %         PS_k = PS_k + coefficient .* PS_2n_minus_1;
        %     elseif m == 2*nonlinearity
        %         PS_k = PS_k + coefficient .* PS_2n;
        %     end
        % end

        if plotsON
            figure('Position', [100, 100, 1600, 300]);
                title("Nonlinear order = %d",nonlinearity);
            subplot(1,4,1); 
                plotFunctions().waveNumberSpectrum(PS_2n, k_range, k_azimuth, ("Spectral Expansion Term 1, P^S_{2n} & n = " + nonlinearity));
                % plotFunctions().waveNumberSpectrum(abs(PS_2n), k_range, k_azimuth, ("Spectral Expansion Term 1, P^S_{2n} & n = " + nonlinearity));
                % xlim([-0.08 0.08]); ylim([-0.08 0.08]);
            subplot(1,4,2);
                plotFunctions().waveNumberSpectrum(PS_2n_minus_1, k_range, k_azimuth, ("Spectral Expansion Term 2, P^S_{2n-1} & n = " + nonlinearity));
                % plotFunctions().waveNumberSpectrum(abs(PS_2n_minus_1), k_range, k_azimuth, ("Spectral Expansion Term 2, P^S_{2n-1} & n = " + nonlinearity));
                % xlim([-0.08 0.08]); ylim([-0.08 0.08]);
            
            subplot(1,4,3);
                plotFunctions().waveNumberSpectrum(PS_2n_minus_2, k_range, k_azimuth, ("Spectral Expansion Term 3, P^S_{2n-2} & n = " + nonlinearity));
                % plotFunctions().waveNumberSpectrum(abs(PS_2n_minus_2), k_range, k_azimuth, ("Spectral Expansion Term 3, P^S_{2n-2} & n = " + nonlinearity));
                % xlim([-0.08 0.08]); ylim([-0.08 0.08]);
                
            subplot(1,4,4);
                plotFunctions().waveNumberSpectrum(sum, k_range, k_azimuth, ("Sum of the nonlinear terms, P^S_{sum of terms} & n = " + nonlinearity));
                % plotFunctions().waveNumberSpectrum(abs(sum), k_range, k_azimuth, ("Sum of the nonlinear terms, P^S_{sum of terms} & n = " + nonlinearity));
        end

         % Update ps_k
        ps_k = ps_k + sum;
    end
    
    % Non-linear spectrum calculation [Eq.50, H&H 1991]
    ps_nl_only = ps_k .* exp(-xi_sqr * k_azimuth.^2) .* (sar_dazimuth / (sar_sub_transect_size*2*pi))^2; % [Eq.45, H&H 1991]; % This factor is for converting from 2pi to spatial domain (dx / (n* 2*pi))^2 it is the FFT and dk term in hHH (SEE NOTES BELOW)
    PS_nl = abs(PS_ql + ps_nl_only);

    if plotsON
        figure('Position', [100, 100, 1200, 300]);
        subplot(1,3,1); plotFunctions().waveNumberSpectrum(ps_nl_only, k_range, k_azimuth, ("Nonlinear Spectral Terms & nonlinearity order = "+nonlinearity));
        
        subplot(1,3,2); plotFunctions().waveNumberSpectrum(PS_ql, k_range, k_azimuth, ("Quasilinear Spectrum, P^S_{ql} & nonlinearity order = "+nonlinearity));
       
        subplot(1,3,3); plotFunctions().waveNumberSpectrum(PS_nl, k_range, k_azimuth, ("Generated SAR Spectrum, P^S_{k} & nonlinearity order = "+nonlinearity));
    end
end


% -------------------------------------------------------------------
%% Generate the first guess SAR Spectrum
% Use an input wave number spectrum to generate the equivalent SAR spectrum

function [PS_k,TS_k, Tv_k, sar_beta,xi_sqr] = generateSARSpectrumFromWaveNumberSpectrum(SatelliteObject, plotsON, transect_number, nonlinearity_order, sar_transect_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k,first_guess_wave_number_spectrum)
    % sar_transect_size typically = 128

    %% SAR metadata attributes
    sar_polarisation = SatelliteObject.Polarisation;
    
    % Get the plarform velocity
    sar_platform_velocity = SatelliteObject.SatelliteVelocity; %m/s from Giacomo's output file.
    
    % Get the incidence angle at the center of the SAR transect.
    sar_incidence_angle_degrees = SatelliteObject.TransectIncidenceAngles(transect_number,2).Value;
    
    % Get the slant range at the scenter of the SAR transect
    sar_transect_slant_range = SatelliteObject.TransectSlantRanges(transect_number,2).Value;

    sar_look_direction = SatelliteObject.LookDirection;

    sar_azimuth_resolution = SatelliteObject.AzimuthResolution;

    % disp("Successfully read in SAR parameters.");
    %% Frozen Surface Contribution pg.3-4 
    % [pg43, Hasselmann, 1990] - "...in all cases as mu = 0.5 1/sec and gamma = 0.
    % This is consistent with field and laboratory measurements (cf. Keller and Wright, 1975)"
    mu = 0.5;
    
    % k in radar look direction variable - convert k in the range component to
    % H&H definition
    [first_guess_kx_range] = klUsingSARlook(sar_look_direction,first_guess_kx_range);

    % Tilt MTF: The change in incidence angle is due to the change in the slope of the wave face.
    Tt_k = tiltMTF(sar_polarisation, sar_incidence_angle_degrees,first_guess_kx_range); % [Eq.5, H&H 1991]
    
    % Hydrodynamic MTF: The interactions between short and long waves, 
    % which modulate the energy and wave number of the short, Bragg scattering, ripple waves
    Th_k = hydrodynamicMTF(first_guess_omega, mu, first_guess_k, first_guess_kx_range,0,0); %[Eq.6, H&H 1991] 
    
    % RAR MTF
    TR_k = rarMTF(Tt_k,Th_k); % [Eq.4, H&H 1991]
    % isequal(TR_k(:,65:end), fliplr(TR_k(:,1:64))) % check for symmetry
    
    % Plots: RAR MTF Development
    % The RAR MTF plot should look something like the H&H Figure 3.
    if plotsON
        figure('Position', [0, 0, 1000, 300]);
        subplot(1,3,1); plotFunctions().generalSpectrumPlots(0,Tt_k, first_guess_kx_range, first_guess_ky_azimuth, "Tilt MTF");

        subplot(1,3,2); plotFunctions().generalSpectrumPlots(0,Th_k, first_guess_kx_range, first_guess_ky_azimuth, "Hydrodynamic MTF");
        
        subplot(1,3,3); plotFunctions().generalSpectrumPlots(0,TR_k, first_guess_kx_range, first_guess_ky_azimuth, "RAR MTF");
    end

    % RAR Image Variance Spectrum calculation
    PR_k = rarImageVarianceSpectrum(TR_k,first_guess_wave_number_spectrum, ...
        rot90(TR_k,2),rot90(first_guess_wave_number_spectrum,2)); % [Eq.13, H&H 1991]
    
    % Plots: RAR Spectrum
    if plotsON
        figure('Position', [0, 0, 300, 300]);
        plotFunctions().generalSpectrumPlots(0,PR_k, first_guess_kx_range, first_guess_ky_azimuth, "RAR Image Variance Spectrum, P^{R}_{k}");
        % xlim([-0.08 0.08]); ylim([-0.08 0.08]);
    end
    
    % disp("Successfully completed Frozen Surface Contribution calculations.");
    %% Motion Effects:
    % "We consider now the modification of the frozen image induced by the
    % surface motion. This is normally described by two effects:
    % the azimuthal displacement, xi, of the apparent position of the
    % backscattering element in the image plane, and an azimuthal smearing
    % or broadening of the image of the (theoretically infinitesimal)
    % backscattering element" [H&H 1991]

    % Simulation of velocity bunching caused by the orbital velocity of
    % the long waves, which produces a Doppler shift in the received, return signal.

    % Range Velocity MTF
    Tv_k = rangeVelocityMTF(first_guess_omega, ...
        sar_incidence_angle_degrees, first_guess_kx_range, first_guess_k); % [Eq.17, H&H 1991]

    % Velocity Bunching Factor - this slightly differs between transects
    % and so needs to be calculated for each transect.
    sar_beta = defineBeta(sar_transect_slant_range, sar_platform_velocity); %[Eq.15, H&H 1991]

    % Velocity Bunching MTF
    Tvb_k = velocityBunchingMTF(sar_beta, first_guess_ky_azimuth, Tv_k);% [Eq.24, H&H 1991]

    % Net SAR imaging MTF
    TS_k = sarImagingMTF(TR_k, Tvb_k); %[Eq.27, H&H 1991]

    % Plots: SAR MTF Development
    if plotsON
        figure('Position', [0, 0, 1600, 300]);
        subplot(1,4,1); plotFunctions().generalSpectrumPlots(0,Tv_k, first_guess_kx_range, first_guess_ky_azimuth, "Range Velocity MTF");
        
        subplot(1,4,2); plotFunctions().generalSpectrumPlots(0,Tvb_k, first_guess_kx_range, first_guess_ky_azimuth, "Velocity Bunching MTF");
        
        subplot(1,4,3); plotFunctions().generalSpectrumPlots(0,TR_k, first_guess_kx_range, first_guess_ky_azimuth, "RAR MTF");
        
        subplot(1,4,4); plotFunctions().generalSpectrumPlots(0,TS_k, first_guess_kx_range, first_guess_ky_azimuth, "SAR imaging MTF");
    end

    % SAR Image Variance Spectrum, PS_k = PS_1 (Linear Mapping Transform)
    % This is the case in which certain conditions are met such that the 
    % SAR image is imaging roughness of the ocean surface linearly – 
    % see Radar Imaging of Ocean Waves textbook pg 117 onwards. 
    % Nonlinearity is then discussed on page 124.
    PS_k_linear = sarImageVarianceSpectrumLinearMappingTransform( ...
        TS_k,first_guess_wave_number_spectrum,rot90(TS_k,2), ...
        rot90(first_guess_wave_number_spectrum,2));%[Eq.26, H&H 1991]
    
    % Plots: Linear SAR Spectrum
    if plotsON
        figure('Position', [0, 0, 300, 300]); 
        plotFunctions().generalSpectrumPlots(0,PS_k_linear, first_guess_kx_range, first_guess_ky_azimuth, "Image Variance Spectrum, P^S_1 ");
        % xlim([-0.08 0.08]); ylim([-0.08 0.08]);

    end
    
    % disp("Successfully completed Motion Effects calculations.");
    %% General Nonlinear mapping
    % The reasoning for doing this is given in H&H on pg 5 in the paragraph before section 3.
    % "To determine the dependence of the SAR image Fourier components IS_k on the wave Fourier components in the general nonlinear case, we first apply a Fourier transform to the basic mapping relation [Eq.20, H&H 1991]... this yields the SAR image variance soectrum PS_k [Eq.30, H&H 1991]"

    % Orbital Velocity Covariance
    fv_r = (orbitalVelocityCovarianceFunction(first_guess_wave_number_spectrum,Tv_k)); % [Eq.43, H&H 1991]

    % Mean Square Azimuthal Displacement
    xi_sqr = meanSquareAzimuthalDisplacement(sar_beta, fv_r); % [Eq.44, H&H 1991]

    % Quasi-linear Approximation / Quasilinear Mapping Transform
    PS_1 = PS_k_linear; % [Eq.55, H&H 1991]
    % Add filter to filter out high frequencies
    [PS_ql,azimuthal_cutoff_factor] = sarImageVarianceSpectrumQuasilinearMappingTransform(first_guess_ky_azimuth, xi_sqr, PS_1); % [Eq.56, H&H 1991]

     % Plots: Quasilinear SAR Spectrum Development
    if plotsON
        figure('Position', [0, 0, 1200, 300]);
        subplot(1,3,1); plotFunctions().generalSpectrumPlots(0,fv_r, first_guess_kx_range, first_guess_ky_azimuth, "Orbital Velocity Covariance Function, f^v(r)");

        subplot(1,3,2); plotFunctions().generalSpectrumPlots(0,azimuthal_cutoff_factor, first_guess_kx_range, first_guess_ky_azimuth, "Azimuthal Cutoff Factor");
        
        subplot(1,3,3); plotFunctions().generalSpectrumPlots(0,PS_ql, first_guess_kx_range, first_guess_ky_azimuth, "Quasilinear Mapping Transform, P^S_{ql}");
    end


    % RAR Image Intensity Autocovariance Function
    fR_r = rarImageIntensityAutocovariance(PR_k); % [Eq.47, H&H 1991]
    % fR_r = ifft2(first_guess_wave_number_spectrum.*abs(TR_k).^2); % WE DO NOT KNOW WHY

    % Covariance Function of RAR Image and Oribital Velocity
    fRv_r = rarImageIntensityCovariance(first_guess_wave_number_spectrum, TR_k, Tv_k, rot90(first_guess_wave_number_spectrum,2), rot90(TR_k,2), rot90(Tv_k,2)); % [Eq.48, H&H 1991]

    % Plots: Autocovariance and covariance functions
    if plotsON
        figure('Position', [0, 0, 1200, 300]);
        subplot(1,3,1); plotFunctions().generalSpectrumPlots(0,fv_r, first_guess_kx_range, first_guess_ky_azimuth, "Orbital Velocity Covariance, f^v(r)");

        subplot(1,3,2); plotFunctions().generalSpectrumPlots(0,fR_r, first_guess_kx_range, first_guess_ky_azimuth, "RAR Image Intensity Autocovariance, f^{R}(r)");

        subplot(1,3,3); plotFunctions().generalSpectrumPlots(0,fRv_r, first_guess_kx_range, first_guess_ky_azimuth, "Covariance of RAR Image and Oribital Velocity, f^{Rv}(r)");
    end

    % Nonlinear Mapping Transform
    % Calculate the Nonlinear Mapping Transform SAR Spectrum using the spectral series expansion terms
    PS_nonlinear = sarImageVarianceSpectrumNonlinearMappingTransform(plotsON, nonlinearity_order, PS_ql, first_guess_ky_azimuth, sar_beta, fv_r, fRv_r, fR_r, xi_sqr, first_guess_kx_range, sar_azimuth_resolution, sar_transect_size);
    
    % dk = (sar_azimuth_resolution / (sar_transect_size*2*pi))^2; % [Eq.45, H&H 1991]
    dk = 1;

    PS_k = PS_nonlinear .* dk;

    % Plots: Generated SAR Spectrum
    if plotsON
        figure('Position', [0, 0, 1200, 300]);
        subplot(1,3,1); plotFunctions().generalSpectrumPlots(0,PS_1, first_guess_kx_range, first_guess_ky_azimuth, "Linear SAR Spectrum, P^S_linear = P^S_1");
 
        subplot(1,3,2); plotFunctions().generalSpectrumPlots(0,PS_ql, first_guess_kx_range, first_guess_ky_azimuth, "Quasilinear SAR Spectrum, P^S_ql");

        subplot(1,3,3); plotFunctions().generalSpectrumPlots(0,PS_k, first_guess_kx_range, first_guess_ky_azimuth, "Generated SAR Spectrum, P^S_k");
    end
    
    % disp("Successfully generated SAR Spectrum.");
end
% -------------------------------------------------------------------


% -------------------------------------------------------------------


function xi = azimuthalDisplacement(beta,v)
    % Azimuthal displacement of the apparent position backscattering
    % element in the image plane
    %[Eq.14, H&H 1991]
    xi = beta * v;
end

function v = orbitalVelocity(orbital_velocity_metadata)
    % This is defined below [Eq.15, H&H 1991] on page 10,716
    % The time average over the period during which the scattering element
    % is viewed by the SAR

    % BUT First order v may be set to the instantaneous orbital velocity in the
    % center of the viewing window.
    v = 0; % NEED TO FIGURE OUT HOW TO DO THIS.
end

function xi_mean = azimuthalDisplacementMean(sar_slant_range, v_all)
    % This is defined midway through the first column on page 10,717 
    % v_all = orbital velocity of the cell facet ensemble

    xi_mean = sar_slant_range * mean(v_all);
end

function [IS_k,IS_k_linear] = sarImageAmplitudeSpectrum(IR_k,Tvb_k,zeta_k,Tvb_k_neg,zeta_k_neg,TS_k,TS_k_neg)
    %[Eq.23, H&H 1991]
    IS_k = IR_k + (Tvb_k*zeta_k +conj(Tvb_k_neg * zeta_k_neg));
    %[Eq.25, H&H 1991]: linear approximation
    IS_k_linear = TS_k*zeta_k +conj(TS_k_neg * zeta_k_neg);
end















