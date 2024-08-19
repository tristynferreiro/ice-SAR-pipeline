% Implementation of H&H 1991.
%{ 
ALL ANGLES ARE IN DEGREES, all trig functions are in degrees

Variable Definitions:
    T_R_k = RAR MTF
    T_t_k = tilt MTF
    T_h_k = hydrodynamic MTF
    k_l = component of the incident wave number vector in the radar look
    direction 
    
%}

% RAR Modulation Transfer Function (MTF)

% Tilt MTF
SAR_polarisation = 'VV';
SAR_incidence_angle_degrees = 0;

