function [S_new,P_new,mag_0,rot_angle_0,J_new] = step0(P_obs,first_guess_k,first_guess_wave_number_spectrum,kx,ky)
%STEP0 Summary of this function goes here
%   Detailed explanation goes here
[temp1, temp2] = find(P_obs==max(P_obs(:)));
k0 = first_guess_k(temp1,temp2);
modk = first_guess_k;
anello  = 1 - 1./(1+(k0./modk).^7); % MIGHT NEED TO ADJUST POWER (5 in this case) FOR DIFFERENT CASES
anello =1;

Jmax_step0   = 1e+10;

% EPS_GUESS = 1.e+10;
rot = -30:1:30;
mag = 0.3:0.1:2;
J_new = zeros(length(rot),length(mag));
for i_rot = 1:length(rot)
    for i_mag = 1:length(mag)
        % for x_shift = -3:1:3
        %     for y_shift = -5:1:5

                % A = circshift(first_guess_wave_number_spectrum, [y_shift x_shift]);

                % S_tmp = i_energy*rot(S_2D_buoy, i_rot, i_mag, /interp)
                S_tmp = mag(i_mag) .* imrotate(first_guess_wave_number_spectrum, rot(i_rot), 'bilinear', 'crop');

                % find(isnan(S_tmp))
                 %     forward_insieme_range, lin_ord, beta1, look_sep, S_tmp, P_tmp, TMP, cshi, fv_00, rg_res
                %     P_tmp = double(P_tmp)
                %     sub = where(P_tmp lt 0d)
                %     if sub[0] ne -1l then P_tmp[sub] = 0d
                 %     P_tmp = P_tmp*BUTTERWORTH
                P_tmp = handh.generateSARSpectrumFromWaveNumberSpectrum(cosmo, 0, nonlinearity_order, sar_sub_transect_size, first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, first_guess_k, S_tmp); % [Eq.65 HH1991]

                J_new(i_rot,i_mag) = trapz(kx, trapz(ky, (P_tmp - P_obs).^2 .* anello ,1),2) ./ sqrt(trapz(kx, trapz(ky, (P_obs).^2 .* anello ,1),2)) ./ sqrt(trapz(kx, trapz(ky, (P_tmp).^2 .* anello ,1),2));


                if J_new(i_rot,i_mag) < Jmax_step0
                    disp("iteration "+i_rot+","+i_mag)
                    Jmax_step0    = J_new(i_rot,i_mag);
                    S_new      = S_tmp;
                    P_new      = P_tmp;
                    % x_shift0   = x_shift;
                    % y_shift0= y_shift;
                    mag_0       = mag(i_mag);
                    rot_angle_0  = rot(i_rot);

                end
        %     end
        % end

    end
end
end


