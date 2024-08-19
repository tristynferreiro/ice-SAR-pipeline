function p_s_k = imageVarianceSpectrum(waveSpectrum,waveSpectrum_inv, Ts_k,Ts_k_inv)
% Equation 26 in Hasselmann
%% Get required metadata
p_s_k = abs(Ts_k).^2.*(waveSpectrum./2) + abs(Ts_k_inv).^2.*(waveSpectrum_inv./2);
end