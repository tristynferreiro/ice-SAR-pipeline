function [S_OUT, P_OUT,J] = inverti_sea_HH91_MIO(S_INP, P_obs1, P_guess, TS_k, first_guess_ky_azimuth,first_guess_kx_range,xi_sqr, first_guess_omega,cosmo,first_guess_k,sar_center_incidence_angle_degrees, sar_sub_transect_size,HandHLibrary)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % INPUTS:
    % S_INP  --> input first guess ocean wave spectrum
    % rg_res --> range resolution
    %
    % OUTPUTS:
    % S_OUT --> output best fit ocean wave spectrum
    % P_OUT --> SAR spectrum corresponding to S_BEST
    
    % % Define constants and parameters
    % grav = 9.81; % acceleration due to gravity
    % k0 = 2 * pi / lc_obs;
    % anello = 1 - 1 / (1 + (k0 / modk)^5);
    % k0 = 2 * pi / 500;
    % butterworth = 1 / (1 + (k0 / modk)^20);
    % beta_ = beta1; % =R/V -> slant range OVER platform velocity
    % 
    % % First guess directional wave spectrum
    % S_guess = S_INP;
    % S0 = S_guess;
    % P_obs = P_obs1;
    % teta = teta1; % incidence angle at the center of image tile [rad]
    % i = 1i;
    % t = look_sep; % time interval between two looks in SAR cross spectra
    % n = nn; % size along each direction
    % dx = ddx; % pixel size in meters
    % lin_ord = lin_ord; % non linearity order in generating the simulated SAR image spectrum
    % 
    % % Load SAR imaging parameters
    % Ti = load('Ti.dat');
    % Ty = load('Ty.dat');
    % Ts_k = Ti - i * beta_ * kaz * Ty;
    % n21 = floor(n / 2) - 1;
    % xx = conj(Ts_k);
    % tmp = circshift(xx, [-n21, -n21]);
    % Ts_k_meno = fft2(tmp) / (n^2);
    % Ts_k_meno = circshift(Ts_k_meno, [n21 + 1, n21 + 1]);
    % 
    % % Initialize variables
    % S_ice = S_guess;
    % mu0 = 1e-1 * max(P_obs)^2;
    % B = 1e-2 * max(S_guess);
    % 
    % % Forward model
    % [P_guess, P_nl, cshi] = forward_insieme_range(lin_ord, beta_, look_sep, S_guess, rg_res);
    % P_guess = P_guess * butterworth;
   
    % J_guess = int_tabulated_2d(kx, ky, anello * P_obs .* (P_obs - P_guess).^2) + ...
    %     mu0 * int_tabulated_2d(kx, ky, anello * (S_ice - S_guess).^2 ./ (B + S_guess).^2);
    
    % k0 = 2 * pi / lc_obs;
    % anello = 1 - 1 / (1 + (k0 / modk)^5);

    P_obs = P_obs1;
    
    % Ring
    [temp1, temp2] = find(P_obs==max(P_obs(:)));
    k0 = first_guess_k(temp1,temp2);
    modk = first_guess_k;
    anello  = 1 - 1./(1+(k0./modk).^5); % MIGHT NEED TO ADJUST POWER (5 in this case) FOR DIFFERENT CASES
    
    % Butterworth
    % butterworth  = 1 - 1./(1+(k0./modk).^20); 
    butterworth =  1 - 1./(1+(k0./(2*modk))^20);% Removing all the contributions from waves 2 x peak wave length 

    kx = first_guess_kx_range(1,:);
    ky = first_guess_ky_azimuth(:,1);
    
    kaz = first_guess_ky_azimuth;
    cshi = xi_sqr;
    n = sar_sub_transect_size;
    omega = first_guess_omega;
    % t = look_sep; % time interval between two looks in SAR cross spectra
    t=0;

    S_guess = S_INP;
    S0 = S_guess;
    
    P_ice = P_obs;


    S_ice = S_guess;


    B = 1e-2 * max(S_guess);
    mu0 = 1e-1 * max(P_obs(:)).^2;
    mu = mu0; % assumption in H&H

    term1 = (P_guess - P_obs).^2 .*  P_obs .* anello;
    term2 = anello .* ( (S_ice - S_guess) ./ (B + S_guess) ).^2;
    
    J_guess_0 = trapz(kx, trapz(ky, term1,1),2) + trapz(kx, trapz(ky, mu .* term2,1),2); % MAYBE just check this trapz()
    % J_guess = trapz(kx, trapz(ky, term1,1),2); % MAYBE just check this trapz()
    
    eps_guess = trapz(kx, trapz(ky, (P_guess - P_obs).^2 .* anello ,1),2) ./ sqrt(trapz(kx, trapz(ky, (P_obs).^2 .* anello ,1),2)) ./ sqrt(trapz(kx, trapz(ky, (P_guess).^2 .* anello ,1),2));
    eps_guess0 = eps_guess;
    eps(1) = eps_guess;

    P_ice = P_guess;
    
    J(1) = J_guess_0;
    Jmax = J(1);

    S_best = 0;
    P_best =0;

    % Iterative inversion
    for iterazione = 1:20
        B = 1e-2 * max(S_guess);
        
        mu_k = mu0 ./ (B + S_guess).^2;
        mu_k = mu_k .* anello;
        mu_k_minus = rot90(mu_k,2);

        % PLOT mu_k IT IS THE FILTER!!!!
        
        % tmp = circshift(mu_k, [-n21, -n21]);
        % tmp = rot90(mu_k,2);
        % mu_k_minus= fft2(tmp) / (n^2);
        % mu_k_minus = circshift(mu_k_meno, [n21 + 1, n21 + 1]);
        % mu_k_minus = rot90(mu_k,2);

        % W_k = 0.5 * abs(Ts_k).^2 .* exp(-kaz.^2 * cshi) .* cos(omega * t);
         % % tmp = circshift(W_k, [-n21, -n21]);
        % W_k_meno = fft2(tmp) / (n^2);
        % W_k_meno = circshift(W_k_meno, [n21 + 1, n21 + 1]);
        W_k = abs(TS_k).^2 .* exp(-kaz.^2 * cshi) .* cos(omega); % [Eq. 75, H&H 1991]
        W_k_minus = rot90(W_k,2); 
        W_k_minus = W_k_minus.*anello;
        W_k = W_k.*anello;
        
        
        % A_k = P_obs .* W_k.^2 + mu_k;
        % tmp = circshift(A_k, [-n21, -n21]);
        % A_k_meno = fft2(tmp) / (n^2);
        % A_k_meno = circshift(A_k_meno, [n21 + 1, n21 + 1]);
        A_k = P_obs .* W_k.^2 +mu_k;
        A_k_minus = rot90(A_k,2);
        A_k = A_k .* anello;
        A_k_minus = A_k_minus .* anello;

        % A_k = W_k.^2 + 2 * mu0; % [Eq. 73, H&H 1991]
        % A_k_minus = rot90(A_k,2); 

        
        B_k = P_obs .* W_k .* W_k_minus .* anello;
        % B_k = W_k .* W_k_minus .* anello;  % [Eq. 74, H&H 1991]
        
        % dP_k = (P_obs - P_ice) .* anello; 
        dP_k = (P_obs - P_ice) .* anello; % [Eq. 71, H&H 1991]


        % tmp = circshift(dF_k, [-n21, -n21]);
        % dF_k_meno = fft2(tmp) / (n^2);
        % dF_k_meno = circshift(dF_k_meno, [n21 + 1, n21 + 1]);
        dF_k = S_guess - S_ice; % [Eq. 72, H&H 1991]
        dF_k_minus = rot90(dF_k,2);
        dF_k = dF_k .* anello;
        dF_k_minus = dF_k_minus .* anello;

        % [Eq. 70, H&H 1991]
        den = (A_k .* A_k_minus - B_k.^2);
        den(den == 0) = 1e-6;
        
        % num = ( A_k_minus .* (W_k .*dP_k + mu0 .* dF_k) ) - ( B_k .* (W_k_minus .* dP_k + mu0 .* dF_minus_k) );
        num = (A_k_minus .* (W_k .* P_obs .* dP_k + mu_k .* dF_k) - B_k .* (W_k_minus .* P_obs .* dP_k + mu_k_minus .* dF_k_minus));

        % deltaF_k = (A_k_minus .* (W_k .* P_obs .* dP_k + mu_k .* dF_k) - ...
            % B_k .* (W_k_meno .* P_obs .* dP_k + mu_k_meno .* dF_k_meno)) ./ den;
        % deltaF_k = deltaF_k .* anello;
        deltaF_k = num ./ den .* anello; % [Eq. 70, H&H 1991]
        
        S_tmp = S_ice + deltaF_k;
        S_tmp(S_tmp < 0) = 0;

        % [P_tmp, P_nl, cshi] = forward_insieme_range(lin_ord, beta_, look_sep, S_tmp, rg_res);
        % P_tmp = P_tmp * butterworth; % CHECK IF WE NEED THIS
        P_tmp = HandHLibrary.generateSARSpectrumFromWaveNumberSpectrum(cosmo, 0, 3, sar_sub_transect_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k,S_tmp);
        P_tmp = P_tmp .* butterworth;
        
        % J(iterazione) = int_tabulated_2d(kx, ky, anello * P_obs .* (P_obs - P_tmp).^2) + ...
        % int_tabulated_2d(kx, ky, anello * mu_k .* (S_tmp - S_guess).^2);
        term1_tmp = (P_tmp - P_obs).^2 .*  P_obs .* anello;
        % term2_tmp = ( (S_tmp - S_guess) ./ (B + S_guess) ).^2 .* anello;
        term2_tmp = anello .* mu_k .*(S_tmp - S_guess).^2;
        J(iterazione+1) = trapz(kx, trapz(ky, term1_tmp,1),2) + trapz(kx, trapz(ky, mu .* term2_tmp,1),2); % MAYBE just check this trapz() 
        % J(iterazione) = trapz(kx, trapz(ky, mu .* term2_tmp,1),2);
        % J(iterazione) = trapz(kx, trapz(ky, term1_tmp,1),2); % MAYBE just check this trapz()
        eps(iterazione+1) = trapz(kx, trapz(ky, (P_guess - P_obs).^2 .* anello ,1),2) ./ sqrt(trapz(kx, trapz(ky, (P_obs).^2 .* anello ,1),2)) ./ sqrt(trapz(kx, trapz(ky, (P_tmp).^2 .* anello ,1),2));

        if eps(iterazione) < eps_guess0
            display(iterazione)
            Jmax = J(iterazione);
            epsmax = eps(iterazione);
            
            S_best = S_tmp;
            P_best = P_tmp;
        end
        S_guess = S_ice;
        S_ice = S_tmp;
        P_ice = P_tmp;
    end
    
    % Output results
    S_OUT = S_best;
    P_OUT = P_best;
    %  S_OUT = 0;
    % P_OUT = 0;
end
S_new = first_guess_wave_number_spectrum;

[S_OUT, P_OUT, J] = inverti_sea_HH91_MIO(S_new,observed_sar_spectrum,first_guess_sar_spectrum, TS_k, first_guess_ky_azimuth ,first_guess_kx_range,xi_sqr, first_guess_omega, cosmo,first_guess_k,sar_incidence_angle_degrees,sar_sub_transect_size,handh);