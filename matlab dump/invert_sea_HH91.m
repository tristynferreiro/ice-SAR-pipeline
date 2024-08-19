% Procedure for SAR inversion in open sea using HH91 method.
% INPUTS:
% P_SAR  --> observed SAR spectrum (complex if LOOK_SEP > 0)
% S_INP  --> input first guess ocean wave spectrum
%
% OUTPUTS:
% S_OUT --> output best fit ocean wave spectrum
% P_OUT --> corresponding SAR spectrum to S_BEST

function [S_OUT, P_OUT] = invert_sea_HH91(S_INP, rg_res)

% Define paths
% Add paths for required directories
addpath('H:\ICE-ARC\Procedure CROSS\Procedure Inversioni SAR\open sea');
addpath('H:\Odden97\vento_marzo');

% Common variables
global nn ddx dk kran kaz modk omega HF_orbital_var
global look_sar band pol Bragg
global Uvento fi_rg
global lin_ord beta1 look_sep
global P_obs1 teta1 fi_wind_azi lc_obs

% Define integration domain for 2D integrals (instead of using TOTAL)
X = reshape(kran(:,1), [], 1);
kx = X;
Y = reshape(kaz(1,:), [], 1);
ky = Y;

grav = 9.81;  % acceleration due to gravity

lc_smoothed = lc_obs;

% Define wavenumber ring for inversion
k0 = 2 * pi / lc_obs;
anello = single(1 - 1 / (1 + (k0 / modk)^5));

% Define Butterworth band-pass filter already applied to the observed spectrum
k0 = 2 * pi / 500;
butterworth = single(1 / (1 + (k0 / modk)^20));

beta_ = beta1;  % = R/V -> slant range OVER platform velocity

% First guess directional wave spectrum
S_guess = S_INP;

S0 = S_guess;
P_obs = P_obs1;

% SAR imaging parameters
teta = teta1;  % incidence angle at the center of image tile [rad]

i = complex(0, 1);
t = look_sep;  % time interval between two looks in SAR cross spectra
               % t = 0 ==> SAR image co-spectrum

% Define SAR image and image spectrum variables
n = nn;       % size along each direction
dx = ddx;     % pixel size in meters
lin_ord = lin_ord;  % non-linearity order in generating the simulated SAR image spectrum

% MTF wind function (assuming it is a MATLAB function defined elsewhere)
% mtf_wind(Uvento, fi_rg, teta);

% Initialize Ti and Ty
Ti = complex(zeros(n, n));
Ty = Ti;

% Load Ti.dat and Ty.dat (assuming they are binary files)
fileID = fopen('Ti.dat', 'r');
Ti = fread(fileID, [n, n], 'double');
fclose(fileID);

fileID = fopen('Ty.dat', 'r');
Ty = fread(fileID, [n, n], 'double');
fclose(fileID);

Ts_k = Ti - i * beta_ * kaz .* Ty;
n21 = int32(n / 2) - 1;

xx = conj(Ts_k);
tmp = circshift(xx, -n21);
Ts_k_meno = fft2(tmp) / (n^2);
Ts_k_meno = circshift(Ts_k_meno, n21 + 1);

S_ice = S_guess;

% Assuming other functions like `INT_TABULATED_2D_MIO`, `forward_insieme_range`, 
% `stimo_cut_off` are defined elsewhere in MATLAB
% for example: INT_TABULATED_2D_MIO(kx, ky, anello .* (P_obs - P_guess).^2);

mu0 = 1e-1 * max(P_obs(:))^2;
B = 1e-2 * max(S_guess(:));

% Call to forward_insieme_range
% [P_guess, P_nl, cshi, fv_00] = forward_insieme_range(lin_ord, beta_, look_sep, S_guess, rg_res);

% Estimate cut-off
% [lc_ice, err_lc, c, y_fit, x_fit, lag] = stimo_cut_off(P_guess, Nc, n_spec, 'pix_size', dx);

P_guess = single(P_guess) .* butterworth;
disp(['Iteration 0: Lambda cut_off = ', num2str(lc_ice)]);

S_ice = S_guess;
J_guess = INT_TABULATED_2D_MIO(kx, ky, anello .* P_obs .* (P_obs - P_guess).^2) + ...
    mu0 * INT_TABULATED_2D_MIO(kx, ky, anello .* (S_ice - S_guess).^2 ./ (B + S_guess).^2);
J_guess_0 = J_guess;

eps_guess = INT_TABULATED_2D_MIO(kx, ky, anello .* (P_obs - P_guess).^2) / ...
    sqrt(INT_TABULATED_2D_MIO(kx, ky, anello .* P_obs.^2)) / ...
    sqrt(INT_TABULATED_2D_MIO(kx, ky, anello .* P_guess.^2));
eps_guess_0 = eps_guess;

P_ice = P_guess;

J = -999 * ones(1000, 1);
eps = J;
iterazione = 0;
J(iterazione + 1) = J_guess_0;  % MATLAB index starts from 1
Jmax = J(iterazione + 1);
eps(iterazione + 1) = eps_guess_0;
epsmax = eps(iterazione + 1);
itera_max = 0;
disp(['Iteration # ', num2str(iterazione), ' J = ', num2str(J(iterazione + 1))]);

S_best = S_guess;
P_best = P_guess;
lc_best = lc_ice;
S_guess = S_ice;

for iterazione = 1:20
    B = 1e-2 * max(S_guess(:));
    mu_k = mu0 ./ (B + S_guess).^2;
    tmp = circshift(mu_k, -n21);
    mu_k_meno = fft2(tmp) / (n^2);
    mu_k_meno = single(circshift(mu_k_meno, n21 + 1));
    mu_k = mu_k .* anello;
    mu_k_meno = mu_k_meno .* anello;

    W_k = 0.5 * abs(Ts_k).^2 .* exp(-kaz.^2 * cshi) .* cos(omega * t);
    tmp = circshift(W_k, -n21);
    W_k_meno = fft2(tmp) / (n^2);
    W_k_meno = single(circshift(W_k_meno, n21 + 1));
    W_k = W_k .* anello;
    W_k_meno = W_k_meno .* anello;

    A_k = P_obs .* W_k.^2 + mu_k;
    tmp = circshift(A_k, -n21);
    A_k_meno = fft2(tmp) / (n^2);
    A_k_meno = single(circshift(A_k_meno, n21 + 1));
    A_k = A_k .* anello;
    A_k_meno = A_k_meno .* anello;

    B_k = P_obs .* W_k .* W_k_meno .* anello;
    dP_k = (P_obs - P_ice) .* anello;
    dF_k = S_guess - S_ice;
    tmp = circshift(dF_k, -n21);
    dF_k_meno = fft2(tmp) / (n^2);
    dF_k_meno = single(circshift(dF_k_meno, n21 + 1));
    dF_k = dF_k .* anello;
    dF_k_meno = dF_k_meno .* anello;

    den = (A_k .* A_k_meno - B_k.^2);
    sub = find(den == 0);
    if ~isempty(sub)
        den(sub) = 1e-6;
    end
    deltaF_k = (A_k_meno .* (W_k .* P_obs .* dP_k + mu_k .* dF_k) - ...
        B_k .* (W_k_meno .* P_obs .* dP_k + mu_k_meno .* dF_k_meno)) ./ den;
    deltaF_k(sub) = 0;
    deltaF_k = deltaF_k .* anello;

    S_tmp = S_ice + single(deltaF_k);
    sub = find(S_tmp < 0);
    if ~isempty(sub)
        S_tmp(sub) = 0;
    end

    % Call to forward_insieme_range
    % [P_tmp, P_nl, cshi, fv_00] = forward_insieme_range(lin_ord, beta_, look_sep, S_tmp, rg_res);

    % Estimate cut-off
    % [lc_ice, err_lc, c, y_fit, x_fit, lag] = stimo_cut_off(P_tmp, Nc, n_spec, 'pix_size', dx);

    P_tmp = single(P_tmp) .* butterworth;
    disp(['Iteration ', num2str(iterazione), ': Lambda cut_off = ', num2str(lc_ice)]);

    J(iterazione + 1) = INT_TABULATED_2D_MIO(kx, ky, anello .* P_obs .* (P_obs - P_tmp).^2) + ...
        mu0 * INT_TABULATED_2D_MIO(kx, ky, anello .* (S_tmp - S_guess).^2 ./ (B + S_guess).^2);

    eps(iterazione + 1) = INT_TABULATED_2D_MIO(kx, ky, anello .* (P_obs - P_tmp).^2) / ...
        sqrt(INT_TABULATED_2D_MIO(kx, ky, anello .* P_obs.^2)) / ...
        sqrt(INT_TABULATED_2D_MIO(kx, ky, anello .* P_tmp.^2));

    J(iterazione + 1) = double(J(iterazione + 1));
    eps(iterazione + 1) = double(eps(iterazione + 1));

    if J(iterazione + 1) < Jmax
        Jmax = J(iterazione + 1);
        epsmax = eps(iterazione + 1);
        itera_max = iterazione;
        S_best = S_tmp;
        P_best = P_tmp;
        lc_best = lc_ice;
    end

    disp(['Iteration # ', num2str(iterazione), ' J = ', num2str(J(iterazione + 1)), ...
        ' Eps = ', num2str(eps(iterazione + 1))]);

    % Check stop criterion
    if abs(J(iterazione + 1) - J_guess_0) / J_guess_0 < 1e-3
        break;
    end

    S_guess = S_tmp;
    S_ice = S_tmp;
    P_ice = P_tmp;
end

disp(['Best iteration # ', num2str(itera_max), ': J = ', num2str(Jmax), ...
    ' Eps = ', num2str(epsmax), ' Lambda cut_off = ', num2str(lc_best)]);
S_OUT = S_best;
P_OUT = P_best;
end
