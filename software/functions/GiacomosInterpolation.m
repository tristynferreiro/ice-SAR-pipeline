function [A_cart1,A_cart2,A_cart3,A_cart4] = GiacomosInterpolation(dx,sar_transect_size, frequency_bins,direction_bins, gravity,wave_spectrum,PASS)
%Giacomo's Interpolation 
%   - based on leggi_wam_csk_20190321_72S_27W.pro and to_cartesian.pro
%
%   PASS = ASC or DESC
%   LOOK = LEFT or RIGHT
%   dx = range resolution (pixel spacing)
%   heading = [degrees north] scene orientation
%   sar_transect_size = size of the filtered SAR transec = 128
%

delta_k = 2 * pi / (dx * sar_transect_size);
tmp = -1 * pi / dx + (0:sar_transect_size) * delta_k;
x1 = tmp(1:sar_transect_size);
y1=x1;

k_range = repmat(x1',1,numel(y1));
k_azimuth = repmat(y1, numel(x1), 1);


if PASS == 'ASC'
    k_azimuth = -1 .* k_azimuth;
end

XX = k_range(:,1);
YY = k_azimuth(1,:);

if PASS == 'ASC'
    YY = -1 .* k_azimuth(1,:);
end

%% Vectors

% Wave number vector
kappa = (2 .* pi .* frequency_bins).^2 ./ gravity; % Deep water approx for wave number, k

% Direction vector [rads]
% direction_bins = (0:35) .* 10 .* pi ./ 180; % This was in Giacomo's code.
theta = direction_bins;
phi = direction_bins .* pi ./ 180;

%% E(f, theta) --> E(k, theta)
n_theta = numel(theta);

Swat = zeros(36,36);
Swat = wave_spectrum;
Swat(~isfinite(Swat)) = 0;

for i=1:n_theta
    % Swat_(:,i) = 2 .* Swat(:,i) .* gravity ./ (8 .* pi.^2 .* frequency_bins .* kappa);
    Swat_(i,:) = 2 .* Swat(i,:) .* gravity ./ (8 .* pi.^2 .* frequency_bins .* kappa);
end

%% Convert to cartesian coordinates
% from to_cartesian.pro

% INPUT:
% A_pol --> spectrum expressed in (k,phi) with dimension [N_k, N_phi]. phi =
%           direction in radians.
% KAPPA --> wavenumber vector [1/m]
% N_k   --> length of kappa vector
% phi    --> direction vector [rads]
% N_phi  --> length of phi vector
% N_cart --> Linear dimension of the output cartesian spectrum. Default =
%           128
% dx    --> Linear spacing of hypothetical SAR image source.

% OUTPUT:
% A_cart --> cartesian spectrum (kx,ky) with dimension [N_cart, N_cart]

A_pol = Swat_';
n_phi = numel(phi);
n_k = numel(kappa);
n_cart = sar_transect_size;

nc = n_phi;
nr = n_k;   % Number radius elements
r_k = kappa; % Polar radius

% Flatten A_pol to 1D array
z = A_pol;
s = zeros(numel(z),1,'single');
s(1:numel(s)) = z(:);

% Flatten phi
phi_out = zeros(nc * nr, 1, 'single');
for i = 1:nr
    phi_out((i-1) * numel(phi) + 1:(i*numel(phi))) = phi;
    theta_out((i-1) * numel(theta) + 1:(i*numel(theta))) = theta; 
end

theta_out = theta_out';

% k interpolation
r_out_k = zeros(nc*nr,1,'single');
count = 0;
for i = 1:nr
    count = count+1;
    r_out_k((i-1)*numel(phi)+1:(i*numel(phi))) = r_k(count);
end

dk = 2 * pi / (n_cart * dx);
% kmin = -1 * pi / (dx * dk);
% kmax = pi / (dx *dk);
kmin = -1 * pi / (dx);
kmax = pi / (dx);

% %% Chat GPT inversion
% % Convert polar to cartesian
% % Define Cartesian grid
% [x_cart, y_cart] = meshgrid(linspace(kmin, kmax, n_cart), linspace(kmin, kmax, n_cart));
% 
% % Convert polar coordinates to Cartesian
% kx = r_out_k .* cos(phi_out);
% ky = r_out_k .* sin(phi_out);
% 
% 
% kx = double(kx);
% ky = double(ky);
% s  = double(s);
% 
% % % Interpolate to Cartesian grid
% F = scatteredInterpolant(kx, ky, s, 'nearest'); % Using natural neighbor interpolation
% A_cart1 = F(x_cart, y_cart);
% 
% %OR
% A_cart2 = griddata(kx, ky, s, x_cart, y_cart, 'nearest');
%% IDL Polar_surface attempted in Matalab
Z = double(s);
R = r_out_k;
Theta = phi_out;
grid_flag = 0;
bounds = [kmin, kmin, kmax,kmax];
spacing = [dk,dk];
use_quintic = 0;
missing = 1;



% Handle zero-radius points to avoid duplicate (0,0) locations
R(R == 0) = max(abs(R)) / 1.0e5;


% Convert to Cartesian coordinates
    if grid_flag
        if size(Z,1) ~= numel(R) || size(Z,2) ~= numel(Theta)
            error('Grid size mismatch: Ensure Z is size [numel(R), numel(Theta)]');
        end
        [Theta_mesh, R_mesh] = meshgrid(Theta, R);
        X = R_mesh .* cos(Theta_mesh);
        Y = R_mesh .* sin(Theta_mesh);
    else
        if numel(R) ~= numel(Z) || numel(Theta) ~= numel(Z)
            error('R and Theta must have the same number of elements as Z');
        end
        X = double(R .* cos(Theta));
        Y = double(R .* sin(Theta));
    end

% Triangulate. 
% The IDL version uses Delaunay triangulation

tr = delaunayTriangulation(X, Y);

% Define output grid bounds
if numel(bounds) < 4
    bounds = [min(X(:)),min(Y(:)),max(X(:)),max(Y(:))];
end

% Define output grid spacing
if numel(spacing) < 2
    spacing = [bounds(3) - bounds(1), bounds(4)-bounds(2)] / 50;
end
    % if isempty(spacing)
    %     grid_size = 51;
    %     dx = (x_max - x_min) / (grid_size - 1);
    %     dy = (y_max - y_min) / (grid_size - 1);
    % else
    %     dx = spacing(1);
    %     dy = spacing(2);
    % end

if missing
    missing_value = 0;
end

% Generate Cartesian grid
[Xq, Yq] = meshgrid(linspace(bounds(1), bounds(3), n_cart), linspace(bounds(2), bounds(4), n_cart));

% Interpolation Attempt 1
F = scatteredInterpolant(tr.Points(:,1), tr.Points(:,2), Z, 'natural', 'none');
A_cart1 = F(Xq, Yq);
A_cart1(isnan(A_cart1)) = missing; % Replace NaNs with missing value

% Interpolation Attempt 2
[Xq, Yq] = meshgrid(bounds(1):spacing(1):bounds(3), bounds(2):spacing(2):bounds(4));
A_cart2 = griddata(X, Y, Z, Xq, Yq, 'cubic');  % 'cubic' is smoother; use 'linear' if needed
A_cart2(isnan(A_cart2)) = missing; % Replace NaNs with missing value

% Interpolation Attempt 3
F = scatteredInterpolant(X, Y, Z, 'natural', 'none'); % 'natural' mimics smooth interpolation
A_cart3 = F(Xq, Yq);
A_cart3(isnan(A_cart3)) = missing;

% Interpolation Attempt 4
method = 'linear';
if use_quintic
    method = 'natural';  % Closest MATLAB equivalent
end

F = scatteredInterpolant(double(X(:)), double(Y(:)), double(Z(:)), method, 'none');
A_cart4 = F(Xq, Yq);
A_cart4(isnan(A_cart4)) = missing_value;

end