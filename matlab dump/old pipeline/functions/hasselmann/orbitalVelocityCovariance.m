function [f_v_r] = orbitalVelocityCovariance(k,waveSpectrum,r,th,Tv_k)

% Need:
% - F(k)
% - T^V_k
% -- omega
% -- theta (Radar incidence angle)
% -- k_l (Wavenumber in look direction)
% --- Look
% -- k
% - k and r (vectors)
% -- r is the x and y coords
% function: f = @(x,y) (x^.2 + y.^2)

% x-axis = SAR flight direction (azimuthal) (heading)

%dk=[0,diff(k)']; % resize to correct size
f_v_r = trapz(waveSpectrum.*abs(Tv_k).^2.*exp(1i.*k.*r)).*dk;

end



