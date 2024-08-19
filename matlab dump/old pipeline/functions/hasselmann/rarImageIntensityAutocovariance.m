function [f_r_r] = rarImageIntensityAutocovariance(k,waveSpectrum,E_k_inv, r, Tr_k, Tr_k_inv)
    
% Need:
% - F(k)
% - T^R_k
% -- T^t_k
%   - k_l (Wavenumber in look direction)
%   -- Look
%   - theta (Radar incidence angle)
%   - Polarisation
% -- T^h_k
%   - omega
%   - mu
%   - k
%   - k_y
% -- theta (Radar incidence angle)
% -- k_l (Wavenumber in look direction)
% - k and r (vectors)
% -- r is the x and y coords
% function: f = @(x,y) (x^.2 + y.^2)

% x-axis = SAR flight direction (azimuthal) (heading)

dk=[0,diff(k)']; % resize to correct size
f_r_r = 0.5.*cumtrapz((waveSpectrum.*abs(Tr_k).^2)+(E_k_inv.*abs(Tr_k_inv).^2).*exp(1i.*k.*r)).*dk;


end