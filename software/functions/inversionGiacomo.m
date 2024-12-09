function [F_best_giacomo_eps, P_best_giacomo_eps, eps, F_best_giacomo_J, P_best_giacomo_J, J] = inversionGiacomo(inversion_iterations,nonlinearity_order, first_guess_wave_number_spectrum, P_obs, TS_k, first_guess_ky_azimuth, first_guess_kx_range,xi_sqr, first_guess_omega,cosmo, first_guess_k, transect_number, sar_sub_transect_size,HandHLibrary)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Ring:
    % This is modification to Giacomo's original approach, here instead the
    % filter characteristics are determined by the peak values in the
    % observed SAR spectrum.
    [temp1, temp2] = find(P_obs==max(P_obs(:)));
    k0 = first_guess_k(temp1,temp2);
    modk = first_guess_k;
    anello  = 1 - 1./(1+(k0./modk).^5); % MIGHT NEED TO ADJUST POWER (5 in this case) FOR DIFFERENT CASES
    anello = 1;
    
    % Butterworth
    % butterworth  = 1 - 1./(1+(k0./modk).^20); 
    % butterworth =  1 - 1./(1+(k0./(2*modk)).^20);% Removing all the contributions from waves 2 x peak wave length 
    % butterworth = 1;

    kx = first_guess_kx_range(1,:);
    % if(cosmo.LookDirection == "right")
    %     kx = fliplr(kx); % needed to make J positive rather than negative
    % end
    ky = first_guess_ky_azimuth(:,1);

    Jmax = 1e15;
    epsmax = 1e15;
    % Iterative inversion
    for n = 1:inversion_iterations

        if(n == 1)
            S_old_guess = first_guess_wave_number_spectrum;
            S_guess = first_guess_wave_number_spectrum;
            
            B = 1e-2 * max(S_old_guess);
            mu0 = 1e-1 * max(P_obs(:)).^2;
            mu_k = mu0;
        else
            B = 1e-2 * max(S_old_guess);
            mu_k = mu0 ./ (B + S_old_guess).^2;
            mu_k = mu_k .* anello;
            mu_k_minus = rot90(mu_k,2);
    
            % PLOT mu_k IT IS THE FILTER!!!!
            
            W_k = abs(TS_k).^2 .* exp(-first_guess_ky_azimuth.^2 * xi_sqr) .* cos(first_guess_omega); % [Eq. 75, H&H 1991]
            W_k_minus = rot90(W_k,2); 
            W_k_minus = W_k_minus.*anello;
            W_k = W_k.*anello;
            
            A_k = P_obs .* W_k.^2 + mu_k; % This is different to H&H:
            % A_k = W_k.^2 + 2 * mu0; % [Eq. 73, H&H 1991]

            A_k_minus = rot90(A_k,2);
            A_k = A_k .* anello;
            A_k_minus = A_k_minus .* anello;
  
            B_k = P_obs .* W_k .* W_k_minus .* anello; % This is different to H&H
            % B_k = W_k .* W_k_minus .* anello;  % [Eq. 74, H&H 1991]
            
            dP_k = (P_obs - P_guess) .* anello; % [Eq. 71, H&H 1991]

            dF_k = S_old_guess - S_guess; % [Eq. 72, H&H 1991]
            dF_k_minus = rot90(dF_k,2);
            dF_k = dF_k .* anello;
            dF_k_minus = dF_k_minus .* anello;
    
            % GIACOMO EQUATION
            denominator = (A_k .* A_k_minus - B_k.^2);
            denominator(denominator == 0) = 1e-6;
            
            % num = ( A_k_minus .* (W_k .*dP_k + mu0 .* dF_k) ) - ( B_k .* (W_k_minus .* dP_k + mu0 .* dF_minus_k) ); 
            num = (A_k_minus .* (W_k .* P_obs .* dP_k + mu_k .* dF_k) - B_k .* (W_k_minus .* P_obs .* dP_k + mu_k_minus .* dF_k_minus));
    
            % deltaF_k = (A_k_minus .* (W_k .* P_obs .* dP_k + mu_k .* dF_k) - ...
                % B_k .* (W_k_meno .* P_obs .* dP_k + mu_k_meno .* dF_k_meno)) ./ den;
            % deltaF_k = deltaF_k .* anello;
            deltaF_k = num ./ denominator .* anello; % [Eq. 70, H&H 1991]
            
            S_guess = S_old_guess + deltaF_k;
            S_guess(S_guess < 0) = 0;
        end

        P_guess = HandHLibrary.generateSARSpectrumFromWaveNumberSpectrum(cosmo, 0, transect_number, nonlinearity_order, sar_sub_transect_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k,S_guess);
        % P_guess = P_guess .* butterworth; % CHECK IF WE NEED THIS

        % J(iterazione) = int_tabulated_2d(kx, ky, anello * P_obs .* (P_obs - P_tmp).^2) + ...
        % int_tabulated_2d(kx, ky, anello * mu_k .* (S_tmp - S_guess).^2);
        term1 = (P_guess - P_obs).^2 .*  P_obs .* anello;
        term2 = mu_k .* ( (S_guess - S_old_guess) ./ (B + S_old_guess) ).^2 .* anello;
        % term2_tmp = anello .* mu_k .* (S_guess - S_old_guess).^2;

        J(n) = trapz(kx, trapz(ky, term1,1),2) + trapz(kx, trapz(ky, term2,2),1);
        eps(n) = trapz(kx, trapz(ky, (P_guess - P_obs).^2 .* anello ,1),2) ./ sqrt(trapz(kx, trapz(ky, (P_obs).^2 .* anello ,1),2)) ./ sqrt(trapz(kx, trapz(ky, (P_guess).^2 .* anello ,1),2));
        
        %% Store the best fit spectra
        if eps(n) < epsmax
            display("EPS iteration " + n)
            epsmax = eps(n);
            
            F_best_giacomo_eps = S_guess;
            P_best_giacomo_eps = P_guess;
        end
        if J(n) < Jmax
            display("J iteration " + n)
            Jmax = J(n);
            
            F_best_giacomo_J = S_guess;
            P_best_giacomo_J = P_guess;
  
        S_old_guess = S_guess;
        end
    end
end