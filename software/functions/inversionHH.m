function [J_eq_62,P_best_eq62,F_best_eq62,J_eq_63,P_best_eq63,F_best_eq63,J_eq_69,P_best_eq69,F_best_eq69] = inversionHH(inversion_iterations, nonlinearity_order, cosmo, handh, plotsON, P_obs, F_first_guess, sar_sub_transect_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k , Tv_k, sar_beta, TS_k)
%H&H1991 Inversion Contains all three inversion equations
%   

    
    kx = first_guess_kx_range(1,:); % This is required for integration
    ky = first_guess_ky_azimuth(:,1); % This is required for integration
    
    % For Eq 63 & 62 cost function
    mu = 0.1 * max(P_obs(:))^2; % [Eq.76 HH1991]
    B = 0.01 * max(F_first_guess(:))^2; % Below [Eq.76 HH1991]
    
    Jmax_eq63 = 1e10;
    Jmax_eq62  = 1e10;
    Jmax_eq69 = 1e10;
    for n = 1:inversion_iterations
        if(n == 1) % Below Eq 64
            F_n_k = F_first_guess;
        else
            % update the spectrum for the next iteration
            F_n_k = F_n_old_k + delta_F_n_old; % [Eq.66 HH1991]
        end
        
        P_n_k = handh.generateSARSpectrumFromWaveNumberSpectrum(cosmo, plotsON, nonlinearity_order, sar_sub_transect_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k, F_n_k); % [Eq.65 HH1991]
    
        % [Eq.62 HH1991]
        term1_62 = (P_n_k - P_obs).^2;
        term2_62 = mu .* ( (F_n_k - F_first_guess) ./ (B + F_first_guess) ).^2;
        J_eq_62(n) = trapz(kx,trapz(ky,term1_62,1),2) + trapz(kx,trapz(ky,term2_62,1),2);
        
        % [Eq.63 HH1991]
        term1_63 = (P_n_k - P_obs).^2 .* P_obs;
        term2_63 = mu .* ( (F_n_k - F_first_guess) ./ (B + F_first_guess) ).^2;
        J_eq_63(n) =  trapz(kx,trapz(ky,term1_63,1),2) + trapz(kx,trapz(ky,term2_63,1),2);
      
        %% Store the best fit spectra
        if(abs(J_eq_62(n)) < abs(Jmax_eq62))
            disp("Eq62: iteration "+n)
            Jmax_eq62 = J_eq_62(n);
            P_best_eq62 = P_n_k;
            F_best_eq62 = F_n_k;
        end
    
        if( abs(J_eq_63(n)) < abs(Jmax_eq63))
            disp("Eq63: iteration "+n)
            Jmax_eq63 = J_eq_63(n);
            P_best_eq63 = P_n_k;
            F_best_eq63 = F_n_k;
        end
    
        %% [Eq 71 --> 75}
        Delta_F_k = F_first_guess - F_n_k; %[Eq.72 HH1991]
        Delta_P_k = P_obs - P_n_k; %[Eq.71 HH1991]
       
        fv_r_n = handh.orbitalVelocityCovarianceFunction(F_n_k, Tv_k);
        xi_sqr_n = handh.meanSquareAzimuthalDisplacement(sar_beta, fv_r_n);
        W_k = abs(TS_k).^2 .* exp(-first_guess_ky_azimuth.^2 .* xi_sqr_n); % [Eq.75 HH1991]
        W_neg_k = rot90(W_k,2); % rotated by 180 degrees
    
        B_k = W_k .* W_neg_k; % [Eq.74 HH1991]
        A_k = W_k.^2 + 2*mu; % [Eq.73 HH1991]
        
        % [Eq.70 HH1991]
        numerator_term1 = rot90(A_k,2) .* (W_k .* Delta_P_k + mu .* Delta_F_k);
        numerator_term2 = B_k .* (W_neg_k .* Delta_P_k + mu .* rot90(Delta_F_k,2));
        denomenator = A_k .* rot90(A_k,2) - B_k.^2;
        delta_F_n = (numerator_term1 + numerator_term2)./denomenator; % [Eq.70 HH1991]
        
        % [Eq.68 HH1991] - Only used for Eq.69
        delta_P_n = 0.5 .* exp( -first_guess_ky_azimuth.^2 .* xi_sqr_n) .* (abs(TS_k).^2 .* Delta_F_k + abs(rot90(TS_k,2)).^2 .* rot90(Delta_F_k,2));
        
        % [Eq.69 HH1991]
        intergral_term1 = (delta_P_n - (P_obs - P_n_k)).^2;
        intergral_term2 = mu .* (delta_F_n - (F_first_guess - F_n_k)).^2;
        J_eq_69(n) = trapz(kx, trapz(ky,intergral_term1,2),1) + trapz(kx, trapz(ky,intergral_term2,2),1) ;  % [Eq.69 HH1991]
        
        %% Store the best fit spectra
        if abs(J_eq_69(n)) < abs(Jmax_eq69)
            disp("Eq96: iteration "+n)
            Jmax_eq69 = J_eq_69(n);
            F_best_eq69 = F_n_k;
            P_best_eq69 = P_n_k;
        end
        
        % update the spectrum for the next iteration
        F_n_old_k = F_n_k; 
        delta_F_n_old = delta_F_n;
    end
end