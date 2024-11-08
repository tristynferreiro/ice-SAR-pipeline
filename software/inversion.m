function [ps_k_nl] = inversion(lin_ord, cosmo,first_guess_kx, first_guess_ky, first_guess_k, first_guess_omega, first_guess_wave_number_spectrum, sar_center_incidence_angle_degrees, sar_sub_transect_size)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %% Vars
    % Import the functions
    handh = hAndH1991Functions;

    % %% Import the SAR metadata attributes
    sar_polarisation = cosmo.Polarisation;
    
    % Get the slant range
    sar_slant_range = cosmo.SlantRange;
    % TODO: Calculate the slant range
    % c = physconst('LightSpeed');
    % slant_range_time_grid = ncread(cosmo.Filepath,"slant_range_time");
    % sar_slant_range_time_transect = slant_range_time_grid(sar_transect_lat_end_index(1):sar_transect_lat_start_index(1),sar_transect_lon_start_index(1):sar_transect_lon_end_index(1)); % center of the transect
    % sar_center_slant_range = sar_slant_range_time_transect(sar_transect_size/2,sar_transect_size/2) .* c/2; % c/2 because it is the RTT of the signal from the antenna --> earth --> antenna
    
    % Get the plarform velocity
    sar_platform_velocity = cosmo.SatelliteVelocity; %m/s from Giacomo's output file.
    % Calculate the platform velocity
    % orbital_vectors_metadata = ncinfo(filepath).Groups.Groups(1).Groups(1).Groups;

    sar_dy = cosmo.AzimuthResolution;
    sar_dx = cosmo.RangeResolution;
   
    %% Frozen Surface:
    % [pg43, Hasselmann, 1990] - "...in all cases as mu = 0.5 1/sec and gamma = 0.
    % This is consistent with field and laboratory measurements (cf. Keller and Wright, 1975)"
    mu = 0.5;
    
    % k in radar look direction variable
    [sar_look,handh_kl] = handh.klUsingSARlook(cosmo.LookDirection,first_guess_kx);
    
    % Tilt MTF
    Tt_k = handh.tiltMTF(sar_polarisation, sar_center_incidence_angle_degrees,handh_kl); % [Eq.5, H&H 1991]

    % Hydrodynamic MTF
    Th_k = ( (first_guess_omega - 1i * mu) ./ (first_guess_omega.^2 + mu^2) ) .* 4.5 .* first_guess_k .* first_guess_omega .* (first_guess_kx.^2 ./ first_guess_k.^2);
    
    % RAR MTF calculation
    TR_k = handh.rarMTF(Tt_k,Th_k); % [Eq.4, H&H 1991]

    % RAR Image Variance Spectrum calculation
    PR_k = handh.rarImageVarianceSpectrum(TR_k,first_guess_wave_number_spectrum,rot90(TR_k,2),rot90(first_guess_wave_number_spectrum,2)); % [Eq.13, H&H 1991]

    %% Motion Effects:
    % Range Velocity MTF
    Tv_k = handh.rangeVelocityMTF(first_guess_omega, sar_center_incidence_angle_degrees, handh_kl, first_guess_k); % [Eq.17, H&H 1991]

    % Velocity Bunching MTF
    sar_beta = handh.beta(sar_slant_range, sar_platform_velocity); %[Eq.15, H&H 1991]

    Tvb_k = -1 .* sar_beta .* first_guess_ky .* Tv_k;% [Eq.24, H&H 1991]

    % Net SAR imaging MTF
    TS_k = handh.sarImagingMTF(TR_k, Tvb_k); %[Eq.27, H&H 1991]
    
    %%
    PS_k_linear = handh.sarImageVarianceSpectrumLinearMappingTransform(TS_k,first_guess_wave_number_spectrum,rot90(TS_k,2),rot90(first_guess_wave_number_spectrum,2));%[Eq.26, H&H 1991]
    
    %% General Nonlinear mapping

    % Orbital velocity covar
    fv_r = ifft2(first_guess_wave_number_spectrum .* abs(Tv_k).^2); % [Eq.43, H&H 1991]

    % RAR image intensity covar
    fR_r = 0.5 * ifft2(PR_k); % [Eq.47, H&H 1991]

    % Covar of ...
    fRv_r = 0.5 .* ifft2(first_guess_wave_number_spectrum.*TR_k.*conj(Tv_k) + rot90(first_guess_wave_number_spectrum,2) .* conj(rot90(TR_k,2)) .* rot90(Tv_k,2) ); % [Eq.48, H&H 1991]

    % Mean square az displacement
    xi_sqr = sar_beta.^2 .* fv_r(1,1); % [Eq.44, H&H 1991]

    % Quasilinear
    PS_1 = PS_k_linear; % [Eq.55, H&H 1991]

    % Add filter to filter out high frequencies
    PS_ql = exp(-1 .* first_guess_ky.^2 .* xi_sqr) .* PS_1; % [Eq.56, H&H 1991]

    %% Spectral series expansion 
    % Initialize variables (assuming they are defined elsewhere in your script)
    fv = fv_r;
    frv = fRv_r;
    frv_r = rot90(fRv_r,2);
    frv_00 = fRv_r(1,1);
    dx = sar_dx;
    n = sar_sub_transect_size;
    n21 = sar_sub_transect_size/2 -1;
    cshi = xi_sqr;
    kaz = first_guess_ky;
    kaz1 = kaz .* sar_beta;
    fr = fR_r;
    
    ps_k = 0;
    
    
    for l = 2:lin_ord
        % tmp1 calculation
        tmp1 = (fv.^l) / factorial(l);
        tmp1 = ifft2(tmp1, 'symmetric');
    
        % a1 calculation
        a1 = (kaz1.^(2*l)) .* tmp1;
    
        % tmp2 calculation
        tmp2 = 1i * (frv - frv_r) .* (fv.^(l-1)) / factorial(l-1);
        tmp2 = ifft2(tmp2, 'symmetric');
    
        % a2 calculation
        a2 = (kaz1.^(2*l-1)) .* tmp2 + a1;
    
        % tmp3 calculation
        tmp3 = fr .* (fv.^(l-1)) / factorial(l-1) + ...
            (frv - frv_00) .* (frv_r - frv_00) .* (fv.^(l-2)) / factorial(l-2);
        tmp3 = ifft2(tmp3, 'symmetric');
    
        % tmp4 calculation
        tmp4 = (kaz1.^(2*l-2)) .* tmp3 + a2;
    
        % Update ps_k
        ps_k = ps_k + tmp4;
    end
    
    % Non-linear spectrum calculation
    ps_2 = ps_k .* exp(-cshi * kaz.^2) * (dx / (n* 2*pi))^2; % This factor is for converting from 2pi to spatial domain (dx / (n* 2*pi))^2 it is the FFT and dk term in hHH (SEE NOTES BELOW)
    ps_k_nl = (PS_ql + ps_2 );

    %

end

