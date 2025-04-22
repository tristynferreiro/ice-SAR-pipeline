%% ==================== INITIALISATION ====================

% Import libraries
waveLibrary = waveLibrary;

% Input parameters
numberOfFrequencyBins = 30;
numberOfDirectionBins = 24;
g = 9.81;        % Gravitational acceleration [m/s²]
n = 128;        % Cartesian grid size
showPlots = 1;
ROTATION_ANGLE = 15; % [DEGREES]

%% ==================== SIMULATED SPECTRUM ====================
Hs = 3.3;                             % Significant wave height
T0 = 7.2;                            % Significant wave period, this dictates the peak frequency, f_peak and w0
mean_wave_direction_degrees = 100;
gammaWaveValueInput = 1.3715;       % Some factor between 1 and 5

% Setup the frequency array like in ERA5 data
f = frequencyBins(numberOfFrequencyBins, "ERA5");
[theta_degrees, theta_rad] = directionBins(numberOfDirectionBins, mean_wave_direction_degrees, "ERA5");


[S_1D_f, S_2D_cosm_f_theta_rad] = waveLibrary.simulateWaveSpectrum(Hs,gammaWaveValueInput,f,T0, mean_wave_direction_degrees, theta_degrees, 0, 'JONSWAP');

if showPlots
    figure('Position', [0, 0, 800, 300]);
    subplot(1,2,1); plotLibrary().waveSpectrum1D(S_1D_f,f,"1D Simulated JONSWAP Wave Spectrum, E(f)");
        
    subplot(1,2,2); 
    plotLibrary().waveSpectrum2D(1, abs(S_2D_cosm_f_theta_rad), f, theta_degrees, 0.5, "Simulated 2D wave spectrum using JONSWAP and  cos^2(theta) model");
end


%% ==================== RUN TESTS ====================
disp(".... STARTING TESTS ....")
test1DWaveSpectrum(f,S_1D_f,[Hs,T0]);

test2DWaveSpectrum(theta_rad,f,S_2D_cosm_f_theta_rad,[Hs,T0,mean_wave_direction_degrees]);

testWaveNumberSpectrum(f,theta_rad,S_2D_cosm_f_theta_rad,[Hs,T0,mean_wave_direction_degrees],showPlots)

rotation_angle_rad = ROTATION_ANGLE .* (pi/180);
testRotation(f,theta_rad,S_2D_cosm_f_theta_rad,rotation_angle_rad, [Hs,T0,mean_wave_direction_degrees], showPlots)

testRotationAndWaveNumber(f,theta_rad,S_2D_cosm_f_theta_rad,rotation_angle_rad, [Hs,T0,mean_wave_direction_degrees], showPlots);

testInterpolationEftheta(f,theta_rad,S_2D_cosm_f_theta_rad,128,[Hs,T0,mean_wave_direction_degrees,gammaWaveValueInput],showPlots)

%% ==================== FUNCTIONS ====================

function f = frequencyBins(numberOfFrequencyBins, Type)
    if(Type == "ERA5")
        if numberOfFrequencyBins<=30
            % Setup the frequency array as per the ERA5 documentation.
            f = zeros(numberOfFrequencyBins,1);
            f(1) = 0.03453;
            for i=2:length(f)
                f(i) = f(i-1)*1.1;
            end
            f = f';
        else 
            % Need to preserve the statistical property regardless of the
            % size of the array. i.e. it still needs to mirror the case as
            % if there were only 24 bins (as we expect in ERA5 data).
            f_max = 0.5478; % Calculated from the 24 bin case
            % Setup the frequency array as per the ERA5 documentation.
            f = zeros(numberOfFrequencyBins,1);
            f_min = 0.03453;
            f(1) = f_min;
            a = (f_max / f_min)^(1 / (numberOfFrequencyBins - 1));

            for i=2:length(f)
                f(i) = f(i-1)*a;
            end
            f = f';
        end
    elseif (Type == "linear")
        % Generic array
         f = linspace(0.01,1,numberOfFrequencyBins);
         f = f';
    end
end

function [theta_degrees, theta_rad] = directionBins(numberOfDirectionBins, mean_wave_direction_degrees, Type)

    if(Type == "ERA5")
        max_degree_bin = 352.5; % ERA5 maximum
        min_degree_bin = 7.5;
        
        theta_degrees = linspace(min_degree_bin,max_degree_bin,numberOfDirectionBins);
        theta_degrees =theta_degrees';

        theta_rad = theta_degrees.* (pi/180);

    elseif (Type == "full circle")
        max_degree_bin = 360; % ERA5 maximum
        min_degree_bin = 0;
        
        bin_spacing = (max_degree_bin-min_degree_bin)/numberOfDirectionBins; % Calculate the grid spacing
        
        
        theta_degrees = linspace(mean_wave_direction_degrees-(bin_spacing*numberOfDirectionBins/2), ...
                                            mean_wave_direction_degrees+(bin_spacing*numberOfDirectionBins/2) ...
                                            ,numberOfDirectionBins);
        theta_rad = theta_degrees.* (pi/180);
    end

end

function [E_kx_ky, Xq, Yq] = polar_surface(R, Theta, Z,sar_azimuth_resolution,sar_range_resolution,interp_size)
    % Convert polar coordinates to Cartesian
    [X, Y, Z_cart] = pol2cart(Theta, R, Z); % TH must be in radians. (according to function docs)
    
    % Swap X and Y because they need to be the correct orientation (Y = cos
    % and X = sin so that the Y axis is azimuth and the X axis is range).
    temp = Y;
    Y=X;
    X=temp;

    % Define Cartesian grid for interpolation
    xq = linspace(min(X(:)), max(X(:)), interp_size); % Need to add in range and azimuth resolution
    yq = linspace(min(Y(:)), max(Y(:)), interp_size);
    [Xq, Yq] = meshgrid(xq, yq);
    
    % Interpolate Z values onto Cartesian grid
    E_kx_ky = griddata(X, Y, Z_cart, Xq, Yq, 'cubic');
    E_kx_ky(isnan(E_kx_ky)) = 0;

    %     % Plot the surface
    % figure;
    % contour(Xq, Yq, E_kx_ky);
    % colorbar;
    % xlabel('X'); ylabel('Y');
    % title('Polar Surface E(k, theta) in Cartesian Coordinates E(kx,ky)');
    % axis equal;
end

function test1DWaveSpectrum(f,S_1D_f,SimInputParams)
    %TEST1DWAVESPECTRUM Compare the parameters of the simulated 1D spectrum
    %to the input parameters
        % S_1D          = E(f)
        % f             = frequency array
        % SimInputParams   = [Hs,T0]

    % Calc parameters of 1D to check that it matches the input
    m0_1D = trapz(f,S_1D_f,2);
    Hs_derived = 4 * sqrt(m0_1D);
    m1_1D = trapz(f,S_1D_f .* f.^1,2);
    Tm_derived = (m1_1D/m0_1D)^(-1); % [Holthuijsen Eq.5] 

    % Display Test results
    disp("------------------------------")
    disp("      Test 1D Spectrum        ")
    disp("------------------------------")
    disp("Param | Input | 1D Spectrum")
    disp("------------------------------")
    disp([' Hs   |   ', num2str(SimInputParams(1)), '   |     ', num2str(Hs_derived)])
    disp([' T0   |   ', num2str(SimInputParams(2)), '  |     -']);
    disp([' Tm   |   -   |    ', num2str(Tm_derived)]);
    
end

function test2DWaveSpectrum(theta_rad,f,S_2D,SimInputParams)
    %TEST2DWAVESPECTRUM Compare the parameters of the simulated 2D spectrum
    %to the input parameters
        % S_2D          = E(f, theta) where theta is in radians
        % f             = frequency array
        % theta_rad     = direction array in radians
        % SimInputParams   = [Hs,T0,meanWaveDirectionDegrees]
   
    % Calc parameters of 2D to check that it matches the input
    [Hs_derived,Tm_derived,direction_derived,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(theta_rad,f,S_2D);
    
    % Display Test results
    disp("------------------------------")
    disp("      Test 2D Spectrum        ")
    disp("------------------------------")
    disp("Param | Input | 2D Spectrum")
    disp("------------------------------")
    disp([' Hs   |   ', num2str(SimInputParams(1)), '   |     ', num2str(Hs_derived)])
    disp([' T0   |   ', num2str(SimInputParams(2)), '  |     -']);
    disp([' Tm   |   -   |    ', num2str(Tm_derived)]);
    disp([' Dir_m |   ', num2str(SimInputParams(3)), '  |     ', num2str(direction_derived)])
    

end

function testWaveNumberSpectrum(f,theta_rad,S_2D,SimInputParams,showPlots)
    %TESTWAVENUMBERSPECTRUM Compare the parameters of the 2D spectrum
    %E(f,theta) after wave number conversion and back to the input 
    %parameters
        % S_2D          = E(f, theta) where theta is in radians
        % f             = frequency array
        % theta_rad     = direction array in radians
        % SimInputParams   = [Hs,T0,meanWaveDirectionDegrees]

    w = f.*(2*pi);
    k = (w).^2 ./ 9.81;

    [E_kx_ky, Jacobian_kx_ky, kx_matrix,ky_matrix] = waveLibrary().waveNumberSpectrum(S_2D,w,k,theta_rad);
    
    [wave_spectrum,Jacobian_f_theta] = waveLibrary().waveSpectrum(E_kx_ky, w, k);
    
    if showPlots
        figure('Position', [0, 0, 1200, 300]);
        subplot(1,3,1); plotLibrary().waveSpectrum2D(1, abs(S_2D), f, ...
            theta_rad.*(180/pi), 0.5, "Simulated 2D wave spectrum E(f,theta)");
        subplot(1,3,2); plotLibrary().waveNumberSpectrum(E_kx_ky, ...
            kx_matrix, ky_matrix, "Wave number spectrum, E(kx,ky)");
        subplot(1,3,3); 
            plotLibrary().waveSpectrum2D(1, abs(wave_spectrum), f, ...
                theta_rad.*(180/pi), 0.5, "Wave Spectrum after conversion " + ...
                "E(Kx,ky) --> E(f,theta)");
   
    
        % Show that the wave spectrum, calculated from the wave number spectrum, looks like the original wave spectrum    
        figure('Position', [0, 0, 800, 300]);
            subplot(1,2,1);
            contour(f,theta_rad.*(180/pi), S_2D);
                xlabel("Frequency"); ylabel('Direction [degrees]');
                title("Final 2-D Wave Spectrum, E(f,theta)", "calculated from E(k)"); 
            hold on;
            contour(f,theta_rad.*(180/pi), wave_spectrum);
                xlabel("Frequency"); ylabel('Direction [degrees]');
                title("Final 2-D Wave Spectrum, E(f,theta)", "calculated from E(k)"); c = colorbar();c.Label.String = '[m^2 s / rad]';
            legend("Original 2D Spectrum","Tested 2D Spectrum (E(k)-->E(f,theta))");
            c = colorbar();c.Label.String = '[m^2 s / rad]';
            subplot(1,2,2);
                plotLibrary().waveSpectrum2D(1, abs(S_2D), f, theta_rad.*(180/pi), 0.5, "Simulated 2D wave spectrum E(f,theta)");
                hold on;
                plotLibrary().waveSpectrum2D(1, abs(wave_spectrum), f, theta_rad.*(180/pi), 0.5, ["Final 2-D Wave Spectrum, E(f,theta)", "calculated from E(k)"]);
                hold off;

    end

    
    % COMPARE the parameter values
    [Hs_derived,Tm_derived,direction_derived,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(theta_rad,f,wave_spectrum);
    
    % Display Test results
    disp("------------------------------")
    disp("  Test E(kx,ky) Calculation   ")
    disp("------------------------------")
    disp("Param | Input | 2D Spectrum   ")
    disp("------------------------------")
    disp([' Hs   |   ', num2str(SimInputParams(1)), '   |     ', num2str(Hs_derived)])
    disp([' T0   |   ', num2str(SimInputParams(2)), '  |     -']);
    disp([' Tm   |   -   |    ', num2str(Tm_derived)]);
    disp([' Dir_m |   ', num2str(SimInputParams(3)), '  |     ', num2str(direction_derived)])
    

end

function testRotation(f,theta_rad,S_2D,rotation_angle_rad, SimInputParams, showPlots)
    %TESTROTATION Apply the rotation angle and then remove it, compare the 
    % spectrum parameters at each point and output.
        % S_2D              = E(f, theta) where theta is in radians
        % f                 = frequency array
        % theta_rad         = direction array in radians
        % rotationAngle     = angle to rotate the spectum by in radians
        % SimInputParams    = [Hs,T0,meanWaveDirectionDegrees]

    adjusted_theta_rad = theta_rad + rotation_angle_rad;

    %COMPARE spectra with angle adjustment
    if showPlots
        figure('Position', [0, 0, 1200, 300]);
        subplot(1,3,1); plotLibrary().waveSpectrumAdjustedToSAR(S_2D, f, theta_rad.*(180/pi), 1, "Original 2-D Wave Spectrum, E(f,theta)");
        subplot(1,3,2); plotLibrary().waveSpectrumAdjustedToSAR(S_2D, f, adjusted_theta_rad.*(180/pi), 1, "Adjusted to SAR scene geometry 2-D Wave Spectrum, E(f,theta)");
        subplot(1,3,3);
            plotLibrary().waveSpectrumAdjustedToSAR(S_2D, f, adjusted_theta_rad.*(180/pi), 1, "");
            hold on;
            plotLibrary().waveSpectrumAdjustedToSAR(S_2D, f, theta_rad.*(180/pi), 1, "Comparison of original vs adjusted spectrum");
    end
       
    
    % COMPARE the parameter values
    [Hs_derived,Tm_derived,direction_derived,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(adjusted_theta_rad,f,S_2D);

    readjusted_theta_rad = adjusted_theta_rad - rotation_angle_rad;
    [Hs_derived_readj,Tm_derived_readj,direction_derived_readj,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(readjusted_theta_rad,f,S_2D);
    
    % Display Test results
    disp("--------------------------------------------")
    disp("         Test 2D Spectrum Rotation          ")
    disp("--------------------------------------------")
    disp("Param | Input | Rotated 2D | Un-rotated 2D ")
    disp("--------------------------------------------")
    disp([' Hs   |   ', num2str(SimInputParams(1)), '   |     ', num2str(Hs_derived), '      |     ',num2str(Hs_derived_readj)])
    disp([' T0   |   ', num2str(SimInputParams(2)), '  |     -' ,'      |     -']);
    disp([' Tm   |   -   |    ', num2str(Tm_derived),'  |   ', num2str(Tm_derived_readj)]);
    disp([' Dir_m |   ', num2str(SimInputParams(3)), '  |     ', num2str(direction_derived), '     |     ', num2str(direction_derived_readj)])

end

function testRotationAndWaveNumber(f,theta_rad,S_2D,rotation_angle_rad, SimInputParams, showPlots)

    w = f.*(2*pi);
    k = (w).^2 ./ 9.81;
    adjusted_theta_rad = theta_rad + rotation_angle_rad;


    [E_kx_ky, ~, kx_matrix,ky_matrix] = waveLibrary().waveNumberSpectrum(S_2D,w,k,theta_rad);
    
    [adjusted_E_kx_ky, ~,adjusted_kx_matrix,adjusted_ky_matrix] = waveLibrary().waveNumberSpectrum(S_2D,w,k, adjusted_theta_rad);
    [adjusted_wave_spectrum, ~] = waveLibrary().waveSpectrum(adjusted_E_kx_ky, w, k);
    
 
    if showPlots
        figure('Position', [0, 0, 1200, 300]);
        subplot(1,3,1); plotLibrary().waveNumberSpectrum(E_kx_ky, kx_matrix, ky_matrix, "Original Wave number spectrum, E(kx,ky)");
        
        subplot(1,3,2); plotLibrary().waveNumberSpectrum(adjusted_E_kx_ky, adjusted_kx_matrix, adjusted_ky_matrix, "Adjusted E(kx,ky)");
            
        subplot(1,3,3); 
            plotLibrary().waveNumberSpectrum(E_kx_ky, kx_matrix, ky_matrix, "");
            hold on;
            plotLibrary().waveNumberSpectrum(adjusted_E_kx_ky, adjusted_kx_matrix, adjusted_ky_matrix, "Comparison Test");
            hold off;
    end
        
    % COMPARE the parameter values
    [Hs_derived,Tm_derived,direction_derived,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(adjusted_theta_rad,f,adjusted_wave_spectrum);

    readjusted_theta_rad = adjusted_theta_rad - rotation_angle_rad;
    [Hs_derived_readj,Tm_derived_readj,direction_derived_readj,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(readjusted_theta_rad,f,adjusted_wave_spectrum);
    
    % Display Test results
    disp("----------------------------------------------")
    disp(" Test 2D Spectrum Rotation & Wave Number Conv ")
    disp("----------------------------------------------")
    disp("Param | Input | Rotated 2D | Un-rotated 2D    ")
    disp("----------------------------------------------")
    disp([' Hs   |   ', num2str(SimInputParams(1)), '   |     ', num2str(Hs_derived), '      |     ',num2str(Hs_derived_readj)])
    disp([' T0   |   ', num2str(SimInputParams(2)), '  |     -' ,'      |     -']);
    disp([' Tm   |   -   |    ', num2str(Tm_derived),'  |   ', num2str(Tm_derived_readj)]);
    disp([' Dir_m |   ', num2str(SimInputParams(3)), '  |     ', num2str(direction_derived), '     |     ', num2str(direction_derived_readj)])
end

function testInterpolationEftheta(f,theta_rad,S_2D,interp_size,SimInputParams,showPlots)
    %TESTINTERPOLATIONEFTHETA want to test interpolation of the ERA5
    %spectrum using simulated data
        % S_2D              = E(f, theta) where theta is in radians
        % f                 = frequency array
        % theta_rad         = direction array in radians
        % interp_size       = size of interpolated spectrum
        % SimInputParams    = [Hs,T0,meanWaveDirectionDegrees,gammaWaveValueInput]


    % Define a new simulated wave spectrum which is the same size as the
    % interpolated wave spectrum this is the IDEAL INTERP SPECTRUM
    ideal_interp_f = frequencyBins(interp_size, "ERA5");
    [ideal_interp_theta_degrees, ideal_interp_theta_rad] = directionBins(interp_size, ...
        SimInputParams(3), "ERA5");

    [~, ideal_interp_S_2D] = waveLibrary().simulateWaveSpectrum( ...
        SimInputParams(1),SimInputParams(4),ideal_interp_f,SimInputParams(2),  SimInputParams(3), ...
        ideal_interp_theta_degrees, 0, 'JONSWAP');

    
    % Now interpolate the original spectrum
    [interp_theta_rad, interp_f] = meshgrid(ideal_interp_theta_rad,ideal_interp_f);
    
    % Interpolate Z values onto Cartesian grid
    interp_S_2D = griddata(f,theta_rad, S_2D, interp_f,interp_theta_rad, 'cubic');
    interp_S_2D(isnan(interp_S_2D)) = 0;
    interp_S_2D(interp_S_2D<0) = 0;
    
    % COMPARE the parameter values
    [Hs_derived,Tm_derived,direction_derived,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(interp_theta_rad(1,:),interp_f(:,1),interp_S_2D);

    [Hs_derived_ideal,Tm_derived_ideal,direction_derived_ideal,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(ideal_interp_theta_rad,ideal_interp_f,ideal_interp_S_2D);
    
    % Display Test results
    disp("----------------------------------------------")
    disp(" Testing Interpolation of Simulated ERA5 data ")
    disp("----------------------------------------------")
    disp("Param | Input | Ideal 'interp' | Interp       ")
    disp("----------------------------------------------")
    disp([' Hs   |   ', num2str(SimInputParams(1)), '   |     ', num2str(Hs_derived_ideal), '          |     ',num2str(Hs_derived)])
    disp([' T0   |   ', num2str(SimInputParams(2)), '  |     -' ,'          |        -']);
    disp([' Tm   |   -   |    ', num2str(Tm_derived_ideal),'      |      ', num2str(Tm_derived)]);
    disp([' Dir_m |   ', num2str(SimInputParams(3)), '  |     ', num2str(direction_derived_ideal), '         |     ', num2str(direction_derived)])

end



