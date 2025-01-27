function [J_step0,F_new, P_new, energy_0, new_kx_0,new_ky_0,rot_angle_0,mag_0,TS_k_0, Tv_k_0, sar_beta_0,xi_sqr_0] = inversionStep0(nonlinearity_order, plotsON,transect_number, P_obs,first_guess_wave_number_spectrum,kx,ky,first_guess_omega,first_guess_k,sar_sub_transect_size, cosmo, handh)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    first_iteration_check = true;
    
    % [temp1, temp2] = find(P_obs==max(P_obs(:)));
    % k0 = first_guess_k(temp1,temp2);
    % modk = first_guess_k;
    % anello  = 1 - 1./(1+(k0./modk).^7); % MIGHT NEED TO ADJUST POWER (5 in this case) FOR DIFFERENT CASES
    anello = 1;
    
    energy_scales = 0.1:1:10;
    wave_number_mag_scales = 0.1:1:10;
    rotation_angles = -30:10:30;

    J_step0 = zeros(length(energy_scales), length(wave_number_mag_scales),length(rotation_angles));
    for i_energy = 1:length(energy_scales)
        for i_mag = 1:length(wave_number_mag_scales)
            for i_rot = 1:length(rotation_angles)
                
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

                % [Eq.78 HH1991]
                B = wave_number_mag_scales(i_mag);
                new_kx = B * (kx .* cosd(rotation_angles(i_rot)) - ky .* sind(rotation_angles(i_rot)));
                new_ky = B * (kx .* sind(rotation_angles(i_rot)) + ky .* cosd(rotation_angles(i_rot)));
                new_k = sqrt(new_kx.^2 + new_ky.^2);
                gravity = 9.81;
                new_omega = sqrt(gravity .* new_k);
                
                % [Eq.77 HH1991]
                theta = -1 * rotation_angles(i_rot); % The imrotate performs counterclockwise rotation.
                S_tmp = imrotate(first_guess_wave_number_spectrum,theta,"bilinear","crop");
                
                S_tmp = energy_scales(i_energy) .* S_tmp;
               
                [P_tmp, TS_k, Tv_k, sar_beta,xi_sqr] = handh.generateSARSpectrumFromWaveNumberSpectrum(cosmo, plotsON, transect_number, nonlinearity_order, sar_sub_transect_size, new_kx, new_ky, new_omega, new_k, S_tmp); % [Eq.65 HH1991]
                
                kx_new = new_kx(1,:);
                ky_new = new_ky(:,1);

                J_step0(i_energy,i_mag,i_rot) = trapz(kx_new, trapz(ky_new, (P_tmp - P_obs).^2 .* anello,1),2) ./ sqrt(trapz(kx_new, trapz(ky_new, (P_obs).^2 .* anello,1),2)) ./ sqrt(trapz(kx_new, trapz(ky_new, (P_tmp).^2  .* anello,1),2));
                
                if first_iteration_check
                    Jmax_step0 = J_step0(i_energy,i_mag,i_rot);
                    first_iteration_check = false;
                end
                    
    
                if J_step0(i_energy,i_mag,i_rot) < Jmax_step0
                    disp("iteration "+i_energy +"," + i_mag+","+i_rot)
                    Jmax_step0    = J_step0(i_energy,i_mag,i_rot);
                    F_new      = S_tmp;
                    P_new      = P_tmp;
                    energy_0     = energy_scales(i_energy);
                    mag_0 = wave_number_mag_scales(i_mag);
                    rot_angle_0  = rotation_angles(i_rot);
                    energy_0_i = i_energy;
                    new_kx_0 = new_kx;
                    new_ky_0 = new_ky;
                    new_k_0 = new_k;
                    new_omega_0 = new_omega;
                    TS_k_0 = TS_k;
                    Tv_k_0 = Tv_k; 
                    sar_beta_0 = sar_beta;
                    xi_sqr_0 = xi_sqr;
                end
           
            end
        end
    end
end