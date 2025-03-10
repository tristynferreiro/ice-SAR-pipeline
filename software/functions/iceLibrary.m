function iceLibrary = iceLibrary
% ICE: All functions used by F's 
%{
    
    References: 
    [1] 2024. Estimation of Antarctic sea ice thickness through observation of
    wave attenuation. Francesca De Santi a,b,∗, Marcello Vichi b,c, Alberto
    Alberello d.
    [2] 2018. SAR image wave spectra to retrieve the thickness of grease‐pancake sea ice using viscous wave
    propagation models. Giacomo De Carolis 1*, Piero Olla2,3 & Francesca De Santi 1
    [3] Pedersen and Coon 2004 A sea ice model for the marginal ice zone
    with an application to the greenland sea.

    * ALL ANGLES ARE IN DEGREES, all trig functions are in degrees

    Variable Definitions:
    
        k = wave number
        hasselmann_wave_number_spectrum = S_w(k) 
               = directional open ocean wave number spectrum at the boundary of
                 the ice-field from the available SAR image derived using 
                 Hasselmann and Hasselmann (1991)
        waves_in_ice_spectrum = S_i(x,k; h,v)
        distance_traveled = delta(e_k,x)
               = the distance traveled from the ice edge to the position x by
                 the wave of wave number k heading in the direction e_k
        wave_damping = q(h,v) 
               = the wave damping predicted by the viscous wave propagation
                 model considered

%}



iceLibrary.dispersionRelationCPModel = @dispersionRelationCPModel; %[Eq.1 [2]]
iceLibrary.dampingRateCalibratedCPModel = @dampingRateCalibratedCPModel; %[Eq.3 [2]]


iceLibrary.wavesInGreasePancakeIceSpectrum = @wavesInGreasePancakeIceSpectrum;
iceLibrary.distanceTraveled = @distanceTraveled;
iceLibrary.waveDamping = @waveDamping;
end
%%
function [k,damping_rate_ki,wave_dispersion_kr] = dispersionRelationCPModel(open_water_wave_number_k0, gravity,sea_ice_thicknesses,sea_ice_viscosities)
    %DISPERSIONRELATIONCPMODEL Approximation of the Close Packing (CP)
    %model dispersion relation
    %   k0 = open water wave number = omega^2/gravity
    %   rho_i = sea ice density constant
    %   rho_w = water density constant
    %   h = sea ice layer thickness
    %   v = sea ice effective viscosity
    %   g = gravity acceleration constant
    %   k = k_r + ik_i
    %       - real(k) = k_r = defines the wave dispersion 
    %       - imag(k) = k_i = defines the wave attenuation/damping rate

    rho_i = 917; %[kg m^-3], constant
    rho_w = 1025; %[kg m^-3], constant
    k0 = open_water_wave_number_k0;

    k = k0 + rho_i/rho_w .* ...
        (k0.^2 .* sea_ice_thicknesses + (1i * sqrt(gravity)/3) .* (sea_ice_thicknesses.^3 ./ sea_ice_viscosities) .*  k0.^(5/2) ...
        ); % [Eq.1 [1]]

    damping_rate_ki = imag(k);
    wave_dispersion_kr = real(k);

end

function [damping_rate_ki, damping_rate_ki_approximate] = dampingRateCalibratedCPModel(sea_ice_thicknesses, frequencies)
    %DAMPINGRATEFORCALIBRATEDCPMODEL Approximation of the Close Packing (CP)
    %model dispersion relation
    %   h = sea ice layer thickness
    %   rho_i = sea ice density constant
    %   rho_w = water density constant
    %   g = gravity acceleration constant

    rho_i = 917; %[kg m^-3], constant
    rho_w = 1025; %[kg m^-3], constant
    eta = 1; % approximation, it is actually 0.963+-0.093 [3]
    A = rho_i/rho_w * 1/(3 * eta) * (2 * pi / sqrt(gravity))^5; % [Eq.3 [1]]
    damping_rate_ki = A .* sea_ice_thicknesses.^(3/2) .* frequencies.^5; % [Eq.3 [1]]


    A_approximate = 10; % [s^5 m^(-5/2)]
    damping_rate_ki_approximate = A_approximate .* sea_ice_thicknesses.^(3/2) .* frequencies.^5; % [Eq.3 [1]]

end

%%
function [waves_in_ice_spectrum] = wavesInGreasePancakeIceSpectrum(hasselmann_wave_number_spectrum, wave_direction)
    %WAVESINICESPECTRUM simulation of the input ocean wave spectrum in sea 
    %   ice in line with [1]

    distance_traveled = distanceTraveled(wave_direction, new_position);
    wave_damping = waveDamping();

    waves_in_ice_spectrum = hasselmann_wave_number_spectrum .* ...
        exp(-2 .* distance_traveled .* wave_damping); % [Eq.1 [2]]
end

function [distance_traveled] = distanceTraveled(wave_direction, new_position)
    %DISTANCETRAVELED calculate the distance traveled from the ice edge to the position x by 
    %   the wave of wave number k heading in the direction e_k
    
    distance_traveled = 1;
 
end

function [wave_damping] = waveDamping()
    %WAVEDAMPING calculate the wave damping predicted by the viscous wave propagation 
    %       model considered

    wave_damping = 1;

end

function [J] = costFunction(observed_sar_spectrum, simulated_sar_spectrum)
    %COSTFUNCTION Comparison of the modeled and the observed SAR spectrum

    J = trapx(k,trapz(distance.^2, abs(observed_sar_spectrum - simulated_sar_spectrum).^2)); % [Eq.2 [2]]

end



