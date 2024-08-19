function [f_rv_r] = rarImageIntensityCovariance(k,waveSpectrum,E_k_inv,r,Tr_k,Tr_k_inv,Tv_k,Tv_k_inv)
    
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
% -- T^V_k
%   - Same as above
% - k and r (vectors)
% -- r is the x and y coords
% function: f = @(x,y) (x^.2 + y.^2)

% x-axis = SAR flight direction (azimuthal) (heading)

dk=[0,diff(k)']; % resize to correct size

f_rv_r = 0.5.*cumtrapz((waveSpectrum.*Tr_k.*conj(Tv_k))+(E_k_inv.*conj(Tr_k_inv).*Tv_k_inv).*exp(1i.*k.*r)).*dk;


end

function [look,incidence_near,incidence_far,num_pixels,polarisation] = getMetadata_f_rv_r(metadata)
    req_atributes = ["antenna_pointing","incidence_near","incidence_far","num_output_lines","mds1_tx_rx_polar","mds2_tx_rx_polar"];
    meta_frvr = filterAttributesNetCDF(metadata.Attributes, req_atributes);
    polarisation= ["";""];
    look = meta_frvr(1).Value;
    incidence_near = meta_frvr(2).Value;
    incidence_far = meta_frvr(3).Value;
    num_pixels = meta_frvr(4).Value;
    polarisation(1) = meta_frvr(5).Value;
    polarisation(2) = meta_frvr(6).Value;
end
