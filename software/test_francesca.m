% SCRIPT TO SIMULATE A WAVE SPECTRUM (Francesca)
clear all
close all
clc

waveLibrary = waveLibrary;
plotLibrary = plotLibrary;

% Grid definition
numberOfFrequencyBins = 24;
numberOfDirectionBins = 24;

f = linspace(0.01, 1, numberOfFrequencyBins)'; % ERA5 Frequency Grid [Hz]
Hs = 15;       % Significant wave height [m]
T0 = 2;       % Significant wave period [s] (defines peak frequency)
mean_wave_direction_degrees = 170; % Mean wave direction
gammaWaveValueInput = 1.3715; % Peakedness factor (1-5)

% Compute the directional bins centered around the mean wave direction
max_degree_bin = 352.5;
min_degree_bin = 7.5;
bin_spacing = (max_degree_bin - min_degree_bin) / numberOfDirectionBins;
direction_bins_degrees = linspace(mean_wave_direction_degrees - (bin_spacing * numberOfDirectionBins / 2), ...
                                  mean_wave_direction_degrees + (bin_spacing * numberOfDirectionBins / 2), ...
                                  numberOfDirectionBins);

% Select the wave spectrum type
spectrumType = 'PM'; % 'PM' or 'JONSWAP'

% Generate the wave spectrum
[S_1D_f, S_2D_cosm] = generateWaveSpectrum(f, direction_bins_degrees, Hs, T0, mean_wave_direction_degrees, gammaWaveValueInput, spectrumType);

% Plot 1D spectrum
figure;
plot(f, S_1D_f, 'b', 'LineWidth', 2);
title("1D Simulated Wave Spectrum, E(f)")
xlabel("Frequency, f (Hz)")
ylabel("Magnitude, E(f) [m^2/Hz]")
grid on;

% Plot 2D spectrum
figure;
plotLibrary.waveSpectrum(1, abs(S_2D_cosm)', f, direction_bins_degrees, 1, "Simulated 2D wave spectrum using JONSWAP and cos²(θ) model");

% Compute final characteristics
[final_Hs_derived, final_Tm_derived, final_direction_derived, final_total_variance_or_energy_derived] = ...
    calculateSpectrumCharacteristics(direction_bins_degrees .* pi / 180, f, S_2D_cosm, 24)




%%%%%%%%%%%%%%%%%%%functions for the library
function [S_1D_f, S_2D] = generateWaveSpectrum(f, direction_bins_degrees, Hs, T0, meanWaveDirection, gamma, spectrumType)
% GENERATEWAVESPECTRUM generates a 1D and 2D ocean wave spectrum.
%
% Inputs:
%   f - Frequency bins [Hz]
%   direction_bins_degrees - Direction bins [degrees]
%   Hs - Significant wave height [m]
%   T0 - Significant wave period [s]
%   meanWaveDirection - Mean wave direction [degrees]
%   gamma - Peakedness factor (only used for JONSWAP)
%   spectrumType - 'PM' for Pierson-Moskowitz, 'JONSWAP' for JONSWAP'
%
% Outputs:
%   S_1D_f - 1D wave spectrum [m^2s/Hz]
%   S_2D - 2D directional wave spectrum [m^2s/Hz/rad]

% Constants
g = 9.81; % Gravitational acceleration [m/s²]
fp = 1 / T0; % Peak frequency [Hz]
f0 = fp / 0.877; % PM scaling correction

% Compute 1D wave spectrum
if strcmp(spectrumType, 'PM')
    % Pierson-Moskowitz Spectrum
    alpha = 8.1e-3;
    beta = 0.74;
    S_1D_f = (alpha * g^2) ./ ((2 * pi * f).^5) .* exp(-beta * (f0^4) ./ (f.^4));
elseif strcmp(spectrumType, 'JONSWAP')
    % JONSWAP Spectrum
    alpha = 5e-3 * Hs^2 * T0^(-4); % Scale parameter
    sigma = zeros(size(f));
    sigma(f <= fp) = 0.07;
    sigma(f > fp) = 0.09;
    G = exp(-((f - fp).^2) ./ (2 * sigma.^2 * fp^2)); % Peakedness factor
    S_1D_f = alpha * g^2 .* f.^(-5) .* exp(-5/4 * (fp ./ f).^4) .* gamma.^G;
else
    error('Invalid spectrum type. Use ''PM'' or ''JONSWAP''.');
end

% Compute total variance (m0) from 1D spectrum
target_m0 = (Hs / 4)^2;
computed_m0 = trapz(f, S_1D_f);

% **Scaling fix**: Ensure the computed variance matches the target variance
scaling_factor = target_m0 / computed_m0;
S_1D_f = S_1D_f * scaling_factor;

% Convert direction to radians
theta_r = direction_bins_degrees * pi / 180;
theta_0 = meanWaveDirection * pi / 180;

% Compute directional spreading function using a cosine power model
D_cosm = cos((theta_r - theta_0) / 2).^2;
D_cosm = D_cosm / trapz(theta_r, D_cosm); % Normalize

% Compute 2D spectrum
S_2D = S_1D_f .* D_cosm;
end



function [Hs_m0,Tm,mean_wave_direction_degrees,total_variance_or_energy] = calculateSpectrumCharacteristics(direction_bins,frequency_bins,wave_spectrum,sar_transect_size)
%CalculateSpectrumCharacteristics Calculate the integral values of the wave
%spectrum
%   wave_spectrum = E(w, theta)

    
% Significant wave height

    % Calculate the 0th order moments (m_n) [Eq.8]
    m0 = trapz(frequency_bins ,trapz(direction_bins, frequency_bins.^0 .* wave_spectrum,2),1); %[Eq.8], dimension 2,3 = direc,freq

    % Calculate the significant wave height (Hs) [Eq.4]
    Hs_m0 = 4*sqrt(m0);

    % Calculate mean wave height (H)
    % H = sqrt(pi/8) * Hs_m0; %[Holthuijsen Eq.4.2.6]
    
% Mean Wave period
    % Calculate the -1st order moments (m_n)
    intermediate = trapz(direction_bins,(1./frequency_bins) .* wave_spectrum,1);
    intermediate(isnan(intermediate))=0;
    m_neg1 = trapz(frequency_bins,intermediate,2); %[Eq.8]
    
    intermediate = trapz(direction_bins,(frequency_bins) .* wave_spectrum,1);
    intermediate(isnan(intermediate))=0;
    m1 = trapz(frequency_bins,intermediate,2); %[Eq.8]
    
    % Calculate the mean wave period (Tm) [Eq.5]
    % Tm = m_neg1/m0; % [Eq.5]
    Tm = (m1/m0)^(-1); % [Holthuijsen Eq.5]
    Tz = (m0/m_neg1)^(1/2); %Tz --Average zero-crossing period [s]
    % T0 = 1./w0; %T0  --Modal (peak) period (T0 = 2 pi /w0) [s]
    Tm = 0;
    
% Mean Wave direction
    % Calculate the angular moments [Eq.9&10]
    a1 = trapz(frequency_bins,trapz(direction_bins,cos(direction_bins).*wave_spectrum,2),1); %[Eq.9]
    b1 = trapz(frequency_bins,trapz(direction_bins,sin(direction_bins).*wave_spectrum,2),1); %[Eq.10]
   
    % Calculate the mean Direction of the wave
    mean_wave_direction = atan(b1./a1); %[Eq.6]
    mean_wave_direction_degrees = mean_wave_direction *180/pi;

    if mean_wave_direction_degrees < 0
        mean_wave_direction_degrees =  180 + mean_wave_direction_degrees; % To get the direction as "coming from"
    else
        mean_wave_direction_degrees = mean_wave_direction_degrees; % To get the direction as "coming from"
    end
    
   
% Spectral Denisty
    % Total variance/energy of the waves = spectral density integrated over all
    % frequencies and wave numbers at specified time of interest
    total_variance_or_energy = trapz(trapz(wave_spectrum,2),1);

end