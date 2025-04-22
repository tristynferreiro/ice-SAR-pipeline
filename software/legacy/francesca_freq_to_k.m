%% ==================== INITIALIZATION ====================
clear; close all; clc;

plotLibrary = plotLibrary; % Load custom plotting functions

% Define physical and grid parameters
dx = 5.0;       % Pixel spacing [m]
n = 512;         % Cartesian grid size
numFreqBins = 30; % Number of frequency bins
numDirBins = 24;  % Number of directional bins
g = 9.81;        % Gravitational acceleration [m/s²]

%% ==================== SIMULATED SPECTRUM ====================
% Define Frequency Bins (Logarithmically spaced)
f(1) = 0.03453;
for i = 2:numFreqBins
    f(i) = f(i - 1) * 1.1;
end
f = f';

% Define Direction Bins
indeg = 360 / numDirBins / 2;
dirBins = linspace(indeg, 360 - indeg, numDirBins); % Direction bins [deg]

% Define Wave Spectrum Parameters
Hs = 8;          % Significant wave height [m]
T0 = 11;         % Significant wave period [s]
meanDir = 120;   % Mean wave direction [deg]
gamma = 1.3715;  % Peakedness factor (JONSWAP)

% Generate 1D & 2D Spectra
[S_1D, S_2D] = generateWaveSpectrum(f, dirBins, Hs, T0, meanDir, gamma, 'JONSWAP');
S_2D(~isfinite(S_2D)) = 0; % Remove NaN values

% Plot 2D Spectrum
figure;
plotLibrary.waveSpectrum2D(1, abs(S_2D)', f, dirBins, 1, "Simulated 2D wave spectrum");

% Compute Initial Significant Wave Height
Hs_init = 4 * sqrt(trapz(f, trapz(dirBins * pi / 180, S_2D, 2), 1));
disp(['Initial Hs: ', num2str(Hs_init)]);

%% ==================== REFINEMENT OF FREQUENCY AND DIRECTION ====================
nf_fine = 128;               % Number of refined frequency bins
numDirBins_fine = 128;       % Number of refined directional bins

% Refine frequency bins using logarithmic interpolation
f_fine = logspace(log10(min(f)), log10(max(f)), nf_fine);
tmp_fine = zeros(nf_fine, numDirBins);

for j = 1:length(dirBins)
    tmp_fine(:, j) = exp(interp1(log(f), log(S_2D(:, j)), log(f_fine), 'pchip', 'extrap'));
end

% Refine direction bins using periodic interpolation (handling 0° = 360°)
dirBins_fine = linspace(0, 360, numDirBins_fine);
S_2D_fine = zeros(nf_fine, numDirBins_fine);

for j = 1:nf_fine
    % Use interp1 with periodic boundary conditions
    S_2D_fine(j, :) = interp1([dirBins-360, dirBins, dirBins+360], ...
        repmat(tmp_fine(j, :), 1, 3), ...
        dirBins_fine, 'pchip','extrap');
end

% Plot refined spectrum
figure;
plotLibrary.waveSpectrum2D(1, abs(S_2D_fine'), f_fine, dirBins_fine, 1, "Refined 2D wave spectrum");

% Compute refined wave height using more accurate integral
Hs_m0 =  4 * sqrt(trapz(f_fine, trapz(dirBins_fine * pi / 180, S_2D_fine, 2), 1));;
disp(['After refinement Hs: ', num2str(Hs_m0)]);

%% ==================== PLOT REFINED VS ORIGINAL 1D SPECTRUM ====================
S_2D_fine(isnan(S_2D_fine)) = 0;
S_1D_fine = trapz(dirBins_fine * pi / 180, S_2D_fine, 2);

figure;
plot(f_fine, S_1D_fine, 'k', 'DisplayName', 'Refined');
hold on;
plot(f, S_1D, 'ro', 'DisplayName', 'Original');
legend;
title('Refined vs Original 1D Spectrum');
set(gca,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
% We can use the refined spectrum
fold=f;
diBinsold=dirBins;
f=f_fine;
dirBins=dirBins_fine;
S_2D=S_2D_fine;
%% ==================== CONVERT TO (k, theta) ====================
% Compute wavenumber (k) using deep-water dispersion relation
k = (2 * pi * f).^2 / g;
theta = deg2rad(dirBins);  % Convert degrees to radians

% TRANSFORM_SPECTRUM Converts a wave spectrum from frequency (f) to wavenumber (k)
% Compute the Jacobian of the transformation df/dk:
df_dk = 1 ./ (4 * pi) .* sqrt(g ./ k);
% S(k, theta) = S(f, theta) * |df/dk|
S_k_theta = S_2D .* df_dk';

Hs_m0 = 4 * sqrt(trapz(k, trapz(theta, S_k_theta, 2), 1));
disp(['In k theta space Hs: ', num2str(Hs_m0)]);
% Plot transformed spectrum
figure;
plotLibrary.waveSpectrum2D(1, abs(S_k_theta'), k, theta * 180 / pi, 0.3, "Wave Spectrum in (k, theta)");

%% ==================== CONVERT TO CARTESIAN COORDINATES (kx, ky) ====================
% Convert to Cartesian (kx, ky) coordinates

[tt, kk] = meshgrid(theta, k);
y = kk(:) .* cos(tt(:));
x =  kk(:) .* sin(tt(:));
z=S_k_theta(:);
% Define Cartesian grid
dk = 2 * pi / (n * dx);
kx=linspace(-pi/dx,pi/dx,n);
kx=kx(kx~=0);
ky=kx;
[Kx, Ky] = meshgrid(kx,ky);
F = scatteredInterpolant(x, y, z, 'linear', 'none');
S_kxky=F(Kx, Ky);

figure;
contour(Kx, Ky, S_kxky);
xlim([min(Kx(:)) max(Kx(:))]);
ylim([min(Kx(:)) max(Kx(:))]);
title('Spettro ricollocato in coordinate cartesiane');
colorbar;

%compute Hs in kx ky space
kappa=sqrt(Kx.^2+Ky.^2);
S_2Dxy=S_kxky./kappa;
% Compute Hs in kx-ky space
Hs_kxky = 4 * sqrt(trapz(Ky(:,1), trapz(Kx(1,:), S_2Dxy, 2)));
disp(['In kx-ky equispaced, Hs: ', num2str(Hs_kxky)]);
% Plot Cartesian spectrum
figure;
contour(Kx, Ky, S_kxky);
xlim([min(Kx(:)) max(Kx(:))]);
ylim([min(Ky(:)) max(Ky(:))]);

%% ==================== REINTERPOLATE TO POLAR GRID ====================
% Recovered variables
% dimension of this grid can be selected as preferred (in principle can
% have any dimension)
k_r=logspace(log10(min(kappa(:))),log10(max(kappa(:))),nf_fine);
theta_r=linspace(0,2*pi,numDirBins_fine);
% k_r = (2 * pi * fold).^2 / g;
% theta_r=diBinsold*pi/180;
for i=1:length(k_r)
    for j=1:length(theta_r)
        kx_r(i,j)=k_r(i)*sin(theta_r(j));
        ky_r(i,j)=k_r(i)*cos(theta_r(j));
    end
end

% Interpolate onto new grid
F = scatteredInterpolant(Kx(:), Ky(:), S_kxky(:), 'natural', 'none');
S_r = F(kx_r, ky_r);
S_r(isnan(S_r))=0;

kappa=sqrt(kx_r.^2+ky_r.^2);
S_rxy=S_r./kappa;
% Compute Hs in kx-ky space
Hs_kxky = 4 * sqrt(trapz(Ky(:,1), trapz(Kx(1,:), S_2Dxy, 2)));
disp(['In kx-ky polar, Hs: ', num2str(Hs_kxky)]);

% Plot reinterpolated spectrum
figure;
contour(kx_r, ky_r, S_r);
xlim([min(Kx(:)) max(Kx(:))]);
ylim([min(Ky(:)) max(Ky(:))]);

Hs_m0 = 4 * sqrt(trapz(k_r, trapz(theta_r, S_r, 2), 1));
disp(['In k theta space Hs: ', num2str(Hs_m0)]);

% INVERSE_TRANSFORM_SPECTRUM Converts a wave spectrum from wavenumber (k) to frequency (f)
f_r = sqrt(k_r * g) / (2 * pi);
dk_df = (8 * pi^2 / g) * f_r';
S_f_r = S_r .*dk_df;

% Plot recovered 2D spectrum
figure;
plotLibrary.waveSpectrum2D(1, abs(S_f_r)', f_r, theta_r * 180 / pi, 1, "Recovered 2D wave spectrum");

% Compute recovered wave height
disp('Recovered Hs');
S_f_r(isnan(S_r)) = 0;
Hs_m0 = 4 * sqrt(trapz(f_r, trapz(theta_r, S_f_r, 2), 1))

% Plot recovered vs original 1D spectrum
S_1Dr = trapz(theta_r, S_f_r, 2);
figure;
plot(f_r, S_1Dr, 'k', 'DisplayName', 'Recovered');
hold on;
plot(fold, S_1D, 'r', 'DisplayName', 'Original');
legend;
title('Recovered vs Original 1D Spectrum');
set(gca,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');

%% ==================== FUNCTIONS ====================

function [S_1D_f, S_2D] = generateWaveSpectrum(f, dirBins, Hs, T0, meanDir, gamma, type)
% Generates 1D & 2D ocean wave spectrum (JONSWAP or Pierson-Moskowitz)

g = 9.81; fp = 1 / T0;
if strcmp(type, 'PM')
    S_1D_f = (8.1e-3 * g^2) ./ ((2 * pi * f).^5) .* exp(-0.74 * (fp ./ f).^4);
else
    alpha = 5e-3 * Hs^2 * T0^(-4);
    sigma = 0.07 * (f <= fp) + 0.09 * (f > fp);
    G = exp(-((f - fp).^2) ./ (2 * sigma.^2 * fp^2));
    S_1D_f = alpha * g^2 .* f.^-5 .* exp(-5/4 * (fp ./ f).^4) .* gamma.^G;
end

% Scale to match variance
S_1D_f = S_1D_f * ((Hs / 4)^2) / trapz(f, S_1D_f);

% Compute 2D spectrum with cos²(θ) spreading
theta_r = deg2rad(dirBins);
D_cosm = cos((theta_r - deg2rad(meanDir)) / 2).^2;
D_cosm = D_cosm / trapz(theta_r, D_cosm); % Normalize
S_2D = S_1D_f .* D_cosm;
end
