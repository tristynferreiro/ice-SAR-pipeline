function handh = hAndH1991Functions
% H&H,1991: All functions used in H&H, 1991
%{

ALL ANGLES ARE IN DEGREES, all trig functions are in degrees
ALL k_neg functions are the function when -k is used.

Variable Definitions:
    zeta_k = wave amplitude of the wave spectrum (the fourier coefficient)

    TR_k = RAR MTF
    TR_k_neg = TR(-k) = additive inverse of TR_k
    Tt_k = tilt MTF
    Th_k = hydrodynamic MTF
    k_l = component of the incident wave number vector in the radar look
            direction 
    mu = damping factor, introduced to describe the response of the short
            waves to the long wave modulation, μ
    Y_r+ 1i*Y_i = complex feedback factor, representing the long wave
            modulation of the wind input to the short waves
    k = wave number vector
    omega = gravity wave frequency, ω=√gk
    F_k = F(k) OR E(k) = ocean wave variance spectrum
    F_k_neg = F(-k) = additive negerse of F_k

    xi = azimuthal displacement of the backscattering element, ξ
    beta = β, velocity bunching parameter
    v = satellite_orbital_velocity OR upsilon, υ OR 
                               = range component of the long wave orbital 
                                 velocity 
                               = time average over the period during which 
                                 the scattering element is viewed by the SAR.
                                 BUT in H&H, first order is equal to the 
                                 INSTANTANEOUS orbital velocity in the 
                                 center of the viewing window 
    sar_platform_velocity = U
    sar_slant_range = rho, ρ OR R OR slant range from satellite
    Tv_k = Tυ_k = range velocity (υ) transfer function
    sar_incidence_angle_degrees = theta, θ in degrees
    TS_k = net SAR imaging MTF
    
    fv_r = f^v(r) = orbital velocity covariance function
    r = the x and y coords vector

%}

handh.klUsingSARlook = @klUsingSARlook;
% The Frozen Surface Contribution / RAR
handh.gravityWaveFrequency = @defineGravityWaveFrequency;
handh.rarMTF = @rarMTF;
handh.tiltMTF = @tiltMTF;
handh.hydrodynamicMTF = @hydrodynamicMTF;
handh.hydrodynamicMTFwithComplexFeedback = @hydrodynamicMTFwithComplexFeedback;
handh.rarImageVarianceSpectrum = @rarImageVarianceSpectrum;

% Motion Effects of SAR
handh.azimuthalDisplacement = @azimuthalDisplacement;
handh.beta = @defineBeta;
handh.rangeVelocityMTF = @rangeVelocityMTF;
handh.sarImageAmplitudeSpectrum = @sarImageAmplitudeSpectrum;
handh.velocityBunchingMTF = @velocityBunchingMTF;
handh.sarImagingMTF = @sarImagingMTF;
handh.sarImageVarianceSpectrumLinearMappingTransform = @sarImageVarianceSpectrumLinearMappingTransform; %SAR image variance spectrum

% General Nonlinear Mapping Relation
handh.orbitalVelocityCovarianceFunction = @orbitalVelocityCovarianceFunction;
handh.rarImageIntensityAutocovariance = @rarImageIntensityAutocovariance;
handh.rarImageIntensityCovariance = @rarImageIntensityCovariance;

handh.meanSquareAzimuthalDisplacement = @meanSquareAzimuthalDisplacement;
handh.spectralExpansion2nTerm = @spectralExpansion2nTerm;
handh.spectralExpansion2nMinus1Term = @spectralExpansion2nMinus1Term;
handh.spectralExpansion2nMinus2Term = @spectralExpansion2nMinus2Term;

end

% From GPT



% From paper

function omega = defineGravityWaveFrequency(g,k)
    omega = sqrt(g*k);
end 

function [look,k_l] = klUsingSARlook(sar_look_metadata,k_y)
    % This is defined on the bottom of page 10,715 after [Eq.6, H&H 1991] 
    if(isequal(sar_look_metadata,"right"))
        look = 0;
        k_l = -1*k_y;
        % k_y = -1*k_y;
    elseif(isequal(sar_look_metadata,"left"))
        look =1;
        k_l = k_y;
        % k_y = k_y;
    end
end

% -------------------------------------------------------------------
%% SAR Imaging of Ocean Waves: Frozen Surface Contribution (pg.10,715 - 10,716)

function TR_k = rarMTF(Tt_k, Th_k)
    %[Eq.4, H&H 1991]
    TR_k = Tt_k + Th_k;
end

function Tt_k = tiltMTF(sar_polarisation, sar_incidence_angle_degrees, k_l)
    %[Eq.5, H&H 1991]
    if(strcmp(sar_polarisation, 'VV') && sar_incidence_angle_degrees<=60)
        Tt_k = 4i .* k_l .* cotd(sar_incidence_angle_degrees) .* (1 + sind(sar_incidence_angle_degrees).^2).^(-1); 
    elseif (strcmp(sar_polarisation, 'HH') && sar_incidence_angle_degrees<=60)
        Tt_k = 8i.* k_l .* (sind(2 .* sar_incidence_angle_degrees)).^(-1); 
    end
end

function Th_k = hydrodynamicMTF(omega, mu, k, k_y)
    %[Eq.6, H&H 1991] 
    % from git 4.5 .* omega .* k_y.^2 .* (omega - mu * 1i) ./ abs(k) .* (omega.^2 + mu.^2)
    Th_k = (omega - (1i * mu)) ./ (omega.^2 + mu.^2) .* 4.5 .* k .* omega .* ((k_y.^2./k.^2));
    % We do not use the complex feedback factor: Y_r+1i*Y_i term.
end

function PR_k = rarImageVarianceSpectrum(TR_k, F_k, TR_k_neg, F_k_neg)
    %[Eq.13, H&H 1991]
    PR_k = 0.5 .* (abs(TR_k).^2.*F_k + abs(TR_k_neg).^2.*F_k_neg);
end

% function Th_k = hydrodynamicMTFwithComplexFeedback(omega, mu, k, k_y, Y_r, Y_i)
%     %[Eq.6, H&H 1991] 
%     Th_k = (omega - (1i * mu)) / (omega^2 + mu^2) * 4.5 * k * omega * ((k_y^2/k^2)+Y_r+1i*Y_i);
% end

% function IR_k = imageModulationIntensity(TR_k,zeta_k,Tr_k_neg,zeta_k_neg)
%     %[Eq.10, H&H 1991]
%     IR_k = TR_k*zeta_k+conj(Tr_k_neg,zeta_k_neg);
% end

% -------------------------------------------------------------------
%% SAR Imaging of Ocean Waves: Motion Effects (pg.10,716 - 10,717)
function Tv_k = rangeVelocityMTF(omega, sar_incidence_angle_degrees, kl, k)
    %[Eq.17, H&H 1991]
    Tv_k = -1 * omega .* ((sind(sar_incidence_angle_degrees) .* (kl ./ abs(k))) + (1i .* cosd(sar_incidence_angle_degrees))); 
end

function beta = defineBeta(sar_slant_range, sar_platform_velocity)
    %[Eq.15, H&H 1991], velocity bunching parameter
    beta = sar_slant_range / sar_platform_velocity;
end

function Tvb_k = velocityBunchingMTF(beta, k_x, Tv_k)
    %[Eq.24, H&H 1991]
    Tvb_k = -1i * beta .* k_x .* Tv_k;
end

function TS_k = sarImagingMTF(TR_k, Tvb_k)
    %[Eq.27, H&H 1991]
    TS_k = TR_k + Tvb_k;
end

function PS_k = sarImageVarianceSpectrumLinearMappingTransform(TS_k, F_k, TS_k_neg, F_k_neg)
    %[Eq.26, H&H 1991]
    PS_k = (abs(TS_k).^2 .* F_k./2) + (abs(TS_k_neg).^2 .* F_k_neg./2);
end

% -------------------------------------------------------------------
%% General nonlinear mapping relation pg.5-6

function fv_r = orbitalVelocityCovarianceFunction(F_k, Tv_k)
    %[Eq.43, H&H 1991]
    fv_r = ifftshift(ifft2(F_k.*abs(Tv_k).^2));
end

function fR_r = rarImageIntensityAutocovariance(F_k,TR_k,F_k_neg,TR_k_neg)
    %[Eq.47, H&H 1991]
    % Should I do the rotation of F_k and TR_k in here rather than passing
    % in the rotated variables?
    fR_r = 0.5 * ifftshift(ifft2(F_k .* abs(TR_k).^2 + F_k_neg .* abs(TR_k_neg).^2));
end

function  fRv_r = rarImageIntensityCovariance(F_k, TR_k, Tv_k, F_k_neg, TR_k_neg, Tv_k_neg)
    %[Eq.48, H&H 1991]
    % Covariance between the RAR image intensity, I(x), and the orbital
    % velocity, v(x)
    fRv_r = 0.5 .* ifftshift(ifft2(F_k .* TR_k .* conj(Tv_k) + F_k_neg .* conj(TR_k_neg) .* Tv_k_neg));
end

function PS_2n = spectralExpansion2nTerm(n, fv_r)
    %[Eq.51, H&H 1991]
    PS_2n = (2*pi)^(-2) .* fft2((fv_r.^n)./factorial(n));
end

function PS_2n_minus_1 = spectralExpansion2nMinus1Term(n,fv_r, fRv_r, fRv_neg_r)
    %[Eq.52, H&H 1991]
    PS_2n_minus_1 = (2*pi)^(-2) .* fft2( ...
        (1i .* (fRv_r - fRv_neg_r) .* fv_r.^(n-1)) ./ factorial(n-1) ...
        );
end

function PS_2n_minus_2 = spectralExpansion2nMinus2Term(n, fv_r, fRv_r, fRv_neg_r, fRv_0, fR_r)
    %[Eq.53, H&H 1991]
    if (n == 1)
        % second_factorial_term = 0;
        PS_2n_minus_2 = (2*pi)^(-2) .* fft2( ...
        (1/factorial(n-1)) .* fR_r .* fv_r.^(n-1) ...
        );
    else
        second_factorial_term = 1/factorial(n-2);
        PS_2n_minus_2 = (2*pi)^(-2) .* fft2( ...
        (1/factorial(n-1)) .* fR_r .* fv_r.^(n-1) + ...
        (second_factorial_term) .* (fRv_r - fRv_0) .* (fRv_neg_r - fRv_0) .* fR_r.^(n-2) ...
        );
    end
    
end




function xi_sqr = meanSquareAzimuthalDisplacement(beta, fv_r)
    %[Eq.44, H&H 1991]
    % The integral term in Eq.44 is equivalent to fv(0) in Eq.43.
    xi_sqr = beta.^2 .* fv_r(1,1); 
end

% -------------------------------------------------------------------

% OLD VERSION - INCORRECT?
% function fv_r = orbitalVelocityCovarianceFunction(sar_kx, sar_ky, kx, ky, k, F_k, Tv_k)
%     %[Eq.43, H&H 1991]
%     fv_r = trapz(kx,trapz(ky, F_k.*abs(Tv_k).^2.*exp(1i.*k.*r),1),2);
%     % % In this equation there is a variable 'r' which, as I understand it is
%     % % the the x and y coord vectors of the SAR image.
%     % r = sqrt(sar_kx.^2+sar_ky.^2); % Calculate the coord matrix
%     % 
%     % % The integration is over k, so I do a double integral over kx and ky
%     % fv_r = trapz(kx,trapz(ky, F_k.*abs(Tv_k).^2.*exp(1i.*k.*r),1),2);
% end
%
% function fR_r = rarImageIntensityAutocovariance(r, k, F_k, TR_k, F_k_neg, TR_k_neg)
%     %[Eq.47, H&H 1991]
%     dk = minus(k(1:end-1,:),k(2:end,:)); % subtract col-wise OR diff(k)
%     fR_r = 0.5 * cumtrapz(dk, (F_k  * abs(TR_k)^2 + F_k_neg * abs(TR_k_neg)^2)*exp(1i*k*r));
% end
%
% function  fRv_r = rarImageIntensityCovariance(r, k, F_k, TR_k, Tv_k, F_k_neg, TR_k_neg, Tv_k_neg)
%     %[Eq.48, H&H 1991]
%     % Covariance between the RAR image intensity, I(x), and the orbital
%     % velocity, v(x)
%     dk = minus(k(1:end-1,:),k(2:end,:)); % subtract col-wise OR diff(k)
%     fRv_r = 0.5 * cumtrapz(dk,(F_k*TR_k*conj(Tv_k)+F_k_neg*conj(TR_k_neg)*Tv_k_neg)*exp(1i*k*r));
% end
%
% function xi_sqr = meanSquareAzimuthalDisplacement(k, beta, Tv_k, F_k)
%     %[Eq.44, H&H 1991]
%     dk = minus(k(2:end),k(1:end-1)); % subtract col-wise OR diff(k)
%     xi_sqr = beta^2 * cumtrapz(dk, abs(Tv_k)^2 * F_k);
% end




function xi = azimuthalDisplacement(beta,v)
    % Azimuthal displacement of the apparent position backscattering
    % element in the image plane
    %[Eq.14, H&H 1991]
    xi = beta * v;
end



function v = orbitalVelocity(orbital_velocity_metadata)
    % This is defined below [Eq.15, H&H 1991] on page 10,716
    % The time average over the period during which the scattering element
    % is viewed by the SAR

    % BUT First order v may be set to the instantaneous orbital velocity in the
    % center of the viewing window.
    v = 0; % NEED TO FIGURE OUT HOW TO DO THIS.
end

function xi_mean = azimuthalDisplacementMean(sar_slant_range, v_all)
    % This is defined midway through the first column on page 10,717 
    % v_all = orbital velocity of the cell facet ensemble

    xi_mean = sar_slant_range * mean(v_all);
end

function [IS_k,IS_k_linear] = sarImageAmplitudeSpectrum(IR_k,Tvb_k,zeta_k,Tvb_k_neg,zeta_k_neg,TS_k,TS_k_neg)
    %[Eq.23, H&H 1991]
    IS_k = IR_k + (Tvb_k*zeta_k +conj(Tvb_k_neg * zeta_k_neg));
    %[Eq.25, H&H 1991]: linear approximation
    IS_k_linear = TS_k*zeta_k +conj(TS_k_neg * zeta_k_neg);
end














