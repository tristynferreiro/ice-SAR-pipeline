% HandH 1991 pipeline
% This is the non-matlab live version

%% Generate the SAR Spectrum
% Import the functions
handh = hAndH1991Functions;

% Run the function
generatedSARSpectrum = handh.generateSARSpectrumFromWaveNumberSpectrum( ...
    cosmo, 1, nonlinearity_order, sar_sub_transect_size, ...
    first_guess_kx_range, first_guess_ky_azimuth, first_guess_omega, ...
    first_guess_k, first_guess_wave_number_spectrum);
