function [J_step0,F_new, P_new, energy_0_i] = inversionStep0(nonlinearity_order, plotsON, P_obs,first_guess_wave_number_spectrum,kx,ky,first_guess_omega,first_guess_k,sar_sub_transect_size, cosmo, handh)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    Jmax_step0   = 1e+10;
    
    % [temp1, temp2] = find(P_obs==max(P_obs(:)));
    % k0 = first_guess_k(temp1,temp2);
    % modk = first_guess_k;
    % anello  = 1 - 1./(1+(k0./modk).^7); % MIGHT NEED TO ADJUST POWER (5 in this case) FOR DIFFERENT CASES
    anello = 1;
    
    rot = -30:1:30;
    energy = 0.1:0.05:1;
    wave_number_mag = 0.1:0.1:1;
    J_step0 = zeros(length(energy), length(wave_number_mag),length(rot));
    for i_energy = 1:length(energy)
        for i_mag = 1:length(wave_number_mag)
            for i_rot = 1:length(rot)
                
                % [Eq.78 HH1991]
                % B = wave_number_mag(i_mag);
                % new_kx = B * (kx .* cosd(rot(i_rot)) - ky .* sind(rot(i_rot)));
                % new_ky = B * (kx .* sind(rot(i_rot)) + ky .* cosd(rot(i_rot)));
                % 
                % new_k = sqrt(new_kx.^2 + new_ky.^2);
                % new_omega = sqrt(gravity .* new_k);
                new_kx = kx;
                new_ky = ky;
                new_omega = first_guess_omega; 
                new_k = first_guess_k;
                
                % [Eq.77 HH1991]
                theta = -1 * rot(i_rot); % The imrotate performs counterclockwise rotation.
                S_tmp = imrotate(first_guess_wave_number_spectrum,theta,"bilinear","crop");
                % S_tmp = interp2(kx, ky, first_guess_wave_number_spectrum, new_kx, new_ky, 'linear', 0);
                S_tmp = energy(i_energy) .* S_tmp;
               
                P_tmp = handh.generateSARSpectrumFromWaveNumberSpectrum(cosmo, plotsON, nonlinearity_order, sar_sub_transect_size, new_kx, new_ky, new_omega, new_k, S_tmp); % [Eq.65 HH1991]
    
                J_step0(i_energy,i_mag,i_rot) = trapz(kx, trapz(ky, (P_tmp - P_obs).^2 .* anello,1),2) ./ sqrt(trapz(kx, trapz(ky, (P_obs).^2 .* anello,1),2)) ./ sqrt(trapz(kx, trapz(ky, (P_tmp).^2  .* anello,1),2));
    
                if J_step0(i_energy,i_mag,i_rot) < Jmax_step0
                    disp("iteration "+i_energy +"," + i_mag+","+i_rot)
                    Jmax_step0    = J_step0(i_energy,i_mag,i_rot);
                    F_new      = S_tmp;
                    P_new      = P_tmp;
                    energy_0     = energy(i_energy);
                    mag_0 = wave_number_mag(i_mag);
                    rot_angle_0  = rot(i_rot);
                    energy_0_i = i_energy;
                    new_kx_0 = new_kx;
                    new_ky_0 = new_ky;
                    new_k_0 = new_k;
                    new_omega_0 = new_omega;
                end
           
            end
        end
    end
end