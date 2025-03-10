function plots = plotLibrary
% Plot  Library 
%   This contains all of the functions relating to plots

plots.sarImage = @sarImage; % single SAR image plotted on pixel axis
plots.sarImageAzRange = @sarImageAzRange; % single SAR image on azimuth-range grid
plots.sarImageLatLon = @sarImageLatLon; % single SAR image on lat-lon grid

plots.sarTransect = @sarTransect; % single SAR image on lat-lon grid
plots.sarTransectWithBuoy = @sarTransectWithBuoy; % single SAR image on lat-lon grid with wave buoy location marker
plots.sarTransectOnSARImage = @sarTransectOnSARImage; % image plotted on another image
plots.sarTransectOnSARImageWithBuoyLocation = @sarTransectOnSARImageWithBuoyLocation;
plots.threeSarTransectsOnSARImageWithBuoyLocation = @threeSarTransectsOnSARImageWithBuoyLocation;

plots.waveSpectrum1D = @waveSpectrum1D; % single wave spectrum E(f)
plots.waveSpectrum2D = @waveSpectrum2D; % single wave spectrum E(f, theta)
plots.waveSpectrumAdjustedToSAR = @waveSpectrumAdjustedToSAR; %
plots.waveNumberSpectrum = @waveNumberSpectrum; % 

plots.generalSpectrumPlots = @generalSpectrumPlots; %
plots.interpolationPlots = @interpolationPlots; %


end

function [] = sarImage(sentinel, sar_data, plot_title)
%sarImage plot the SAR image as pixels  

    if sentinel
        imshow(sar_data);
        title(plot_title);
        % title("Sentinel 1A (C-band) Calibrated SAR image - DD/MM/YYYY")
        ylabel("Across Range Pixels");
        xlabel("Along Range Pixels");
    else
        pcolor(sar_data);
        clim([0 2]);
        shading("flat");
        colormap("gray");
        title("Sentinel 1A (C-band) Calibrated SAR image - DD/MM/YYYY")
        ylabel("Across Range Pixels");
        xlabel("Along Range Pixels");
    end
end

function [] = sarImageAzRange(sentinel, sar_data, sar_azimuth_resolution, sar_range_resolution, plot_title)
%sarImageAzRange plot the SAR image on the azimuth-range grids    
    
    if sentinel
        imshow(sar_data);
        title(plot_title);
        % title("Sentinel 1A (C-band) Calibrated SAR image - DD/MM/YYYY")
        ylabel("Across Range Pixels");
        xlabel("Along Range Pixels");
    else
        pcolor(sar_data);
        clim([0 2]);
        shading("flat");
        colormap("gray");
        title("Sentinel 1A (C-band) Calibrated SAR image - DD/MM/YYYY")
        ylabel("Across Range (m)");
        xlabel("Along Range (m)");
    end

end

function [] = sarImageLatLon(sentinel, sar_data, latitude_grid, longitude_grid, plot_title)
%sarImageLatLon plot the SAR image on the latitude-longitude grids

    if sentinel
        imshow(sar_data);
        title(plot_title);
        % title("Sentinel 1A (C-band) Calibrated SAR image - DD/MM/YYYY")
        ylabel("Across Range Pixels");
        xlabel("Along Range Pixels");
    else
        pcolor(sar_data,longitude_grid,latitude_grid);
        clim([0 2]);
        shading("flat");
        colormap("gray");
        title("Sentinel 1A (C-band) Calibrated SAR image - DD/MM/YYYY")
        xlabel("Latitude");
        ylabel("Longitude");
    end

end

function [] = sarTransect(sar_transect_data, transect_latitude_grid, transect_longitude_grid, plot_title)
%sarTransect plot the SAR transect on latitude-longitude grids
    
    sar_transect_data(isnan(sar_transect_data))=0;

    pcolor(transect_longitude_grid, transect_latitude_grid, sar_transect_data./mean(sar_transect_data(:)));
    shading("flat"); colormap("jet"); clim([0 2]);
    title(plot_title);
    xlabel("Longitude"); ylabel("Latitude");

end

function [] = sarTransectWithBuoy(sar_transect_data, transect_latitude_grid, transect_longitude_grid, wave_buoy_lon, wave_buoy_lat, plot_title)
%sarTransect plot the SAR image on lat-lon grid with wave buoy location marker
    
    sar_transect_data(isnan(sar_transect_data)) = 0;

    pcolor(transect_longitude_grid, transect_latitude_grid, sar_transect_data./mean(sar_transect_data(:)));
    shading("flat"); colormap("gray"); clim([0 2]);
    title(plot_title);
    xlabel("Longitude"); ylabel("Latitude");
    hold on;
    plot(round(wave_buoy_lon,2), round(wave_buoy_lat,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    hold off;

end

function [] = sarTransectOnSARImage(sar_data, latitude_grid, longitude_grid, sar_transect_data, transect_latitude_grid, transect_longitude_grid, plot_title)
%sarTransect plot the SAR transect on the original image
%   Note: if the original data is too large, the pcolor will not render.
%   The original image will then need to be transected (around the area of interest)
%   before passing it into this function

    sar_data(isnan(sar_data)) = 0;
    sar_transect_data(isnan(sar_transect_data)) = 0;

    pcolor(longitude_grid, latitude_grid, sar_data./mean(sar_data(:)));
    hold on; 
    pcolor(transect_longitude_grid, transect_latitude_grid, sar_transect_data./mean(sar_transect_data(:))); 
    shading("flat");
    colormap("gray"); clim([0 2]);
    xlabel("Longitude"); ylabel("Latitude");
    title(plot_title);
    hold off;

end

function [] = sarTransectOnSARImageWithBuoyLocation(sar_data, latitude_grid, longitude_grid, sar_transect_data, transect_latitude_grid, transect_longitude_grid,wave_buoy_lat,wave_buoy_lon, plot_title)
%sarTransect plot the SAR transect on the original image
%   Note: if the original data is too large, the pcolor will not render.
%   The original image will then need to be transected (around the area of interest)
%   before passing it into this function

    sar_data(isnan(sar_data)) = 0;
    sar_transect_data(isnan(sar_transect_data)) = 0;

    pcolor(longitude_grid, latitude_grid, sar_data./mean(sar_data(:))); 
    hold on; 
    pcolor(transect_longitude_grid, transect_latitude_grid, sar_transect_data./mean(sar_transect_data(:))); 
    shading("flat"); colormap("gray"); clim([0 2]);
    xlabel("Longitude"); ylabel("Latitude");
    title(plot_title);

    plot(round(wave_buoy_lon,3), round(wave_buoy_lat,3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    hold off;
end

function [] = threeSarTransectsOnSARImageWithBuoyLocation(sar_data, latitude_grid, longitude_grid, sar_transect_data_1, transect_latitude_grid_1, transect_longitude_grid_1, latitude_of_interest_1 ,longitude_of_interest_1, sar_transect_data_2, transect_latitude_grid_2, transect_longitude_grid_2, latitude_of_interest_2 ,longitude_of_interest_2, sar_transect_data_3, transect_latitude_grid_3, transect_longitude_grid_3, latitude_of_interest_3 ,longitude_of_interest_3, plot_title)
%sarTransect plot the SAR transect on the original image
%   Note: if the original data is too large, the pcolor will not render.
%   The original image will then need to be transected (around the area of interest)
%   before passing it into this function
    
    sar_data(isnan(sar_data)) = 0;
    sar_transect_data_1(isnan(sar_transect_data_1)) = 0;
    sar_transect_data_2(isnan(sar_transect_data_2)) = 0;
    sar_transect_data_3(isnan(sar_transect_data_3)) = 0;

    pcolor(longitude_grid, latitude_grid, sar_data./mean(sar_data(:))); 
    shading("flat"); colormap("gray"); clim([0 2]);
    xlabel("Longitude"); ylabel("Latitude");
    title(plot_title);

    hold on; 
    pcolor(transect_longitude_grid_1, transect_latitude_grid_1, sar_transect_data_1./mean(sar_transect_data_1(:))); 
    pcolor(transect_longitude_grid_2, transect_latitude_grid_2, sar_transect_data_2./mean(sar_transect_data_2(:))); 
    pcolor(transect_longitude_grid_3, transect_latitude_grid_3, sar_transect_data_3./mean(sar_transect_data_3(:))); 
    % shading("flat"); colormap("jet"); clim([0 2]);
    

    plot(round(longitude_of_interest_1,3), round(latitude_of_interest_1,3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(round(longitude_of_interest_2,3), round(latitude_of_interest_2,3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot(round(longitude_of_interest_3,3), round(latitude_of_interest_3,3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    
    % Show where ERA5 grid point seperation is
    % [a] = find(sar_latGrid(min_lat_index:max_lat_index, min_lon_index:max_lon_index) == -34)
    x = [17.9 18.5]; y = [-34 -34];
    line(x,y,'Color',"#EDB120",'LineStyle','--');

    x = [18 18]; y = [-34.3 -33.3];
    line(x,y,'Color',"#EDB120",'LineStyle','--');
    
    hold off;
    legend("","","","","Cape Point Buoy","Test location 2","Test location 3","ERA5 data boundary lines");

end



%%
function [] = waveSpectrum1D(wave_spectrum, frequency_bins, plot_title)
%waveSpectrum plot the wave spectrum E(f, theta)

    plot(frequency_bins,wave_spectrum)
    title(plot_title)
    xlabel("Frequency, f (Hz)")
    ylabel("Magnitude, E(f) [m^2/Hz/rad]")
end

function [] = waveSpectrum2D(polar, wave_spectrum, frequency_bins, direction_bins, max_frequency_to_plot, plot_title)
%waveSpectrum plot the wave spectrum E(f, theta)

    if polar
        PolarContour(wave_spectrum, frequency_bins, max_frequency_to_plot, direction_bins, 40);
        title(plot_title);
    else
        contour(frequency_bins,direction_bins,wave_spectrum,20);
        xlabel("Frequency"); ylabel('Direction [degrees]');
        title(plot_title); 
        c = colorbar();c.Label.String = '[m^2/Hz/rad]';
    end
end


function [] = waveSpectrumAdjustedToSAR(wave_spectrum, frequency_bins, direction_bins, max_frequency_to_plot, plot_title)
%waveSpectrum version of wave spectrum plot that is adjusted to SAR image
%orientation

% TO DO: Reformat to  -> 360 and add lines for every 30 degrees

    % NEED TO ADD AZ AND RANGE AXIS LABELS
    PolarContour(wave_spectrum, frequency_bins, max_frequency_to_plot, direction_bins, 40);
    title(plot_title);
end

function [] = waveNumberSpectrum(wave_number_spectrum, kx_range, ky_azimuth, plot_title)
%waveNumberSpectrum plot the wave number spectrum

    contour(kx_range, ky_azimuth, wave_number_spectrum,40);
    xlabel("k_{x = range}"); ylabel("k_{y = azimuth}"); title(plot_title)
    c = colorbar();c.Label.String = '[units?]';

    % xlim([-0.05 0.05]); ylim([-0.05 0.05]); % FIND CLEVER WAY TO DO THIS
end

%%
function [] = generalSpectrumPlots(dB, data, kx_range, ky_azimuth, plot_title)
%generalSpectrumPlots plot the spectra in the H&H 1991 process
    
    if dB
        contour(kx_range, ky_azimuth, 20*log10(data));
        xlabel("k_{x = range}"); ylabel("k_{y = azimuth}"); title(plot_title);
        c = colorbar(); c.Label.String = 'dB Scale'; 
        grid on;
        % xlim([-0.08 0.08]); ylim([-0.08 0.08]);
    else
        contour(kx_range, ky_azimuth, data,40);
        xlabel("k_{x = range}"); ylabel("k_{y = azimuth}"); title(plot_title); 
        c = colorbar(); c.Label.String = 'Linear Scale';
        grid on;
        % xlim([-0.08 0.08]); ylim([-0.08 0.08]);
    end

end

function [] = interpolationPlots(era5_k,era5_direction_bins_degrees,era5_kx_matrix,era5_ky_matrix,era5_wave_number_spectrum,interp_k_matrix,interp_kx_matrix,interp_ky_matrix,interp_E_kx_ky)
    figure('Position', [0, 0, 1600, 900]);
    subplot(3,4,1)
    plot(era5_k);
    title("ERA5 k"); xlabel('bin number'); ylabel('magnitude'); colorbar;
    subplot(3,4,2);
    contour(era5_k, era5_direction_bins_degrees,era5_kx_matrix,40);
    title("ERA5 kx"); xlabel('ERA5 k'); ylabel('Direction [degrees]'); colorbar;
    subplot(3,4,3);
    contour(era5_k, era5_direction_bins_degrees, era5_ky_matrix,40);
    title("ERA5 ky"); xlabel('ERA5 k'); ylabel('Direction [degrees]'); colorbar;
    subplot(3,4,4); % Positive 2-D wave number spectrum E(kx,ky)=E(k)
    contour(era5_kx_matrix, era5_ky_matrix, era5_wave_number_spectrum,40);
    xlim([-0.05 0.05]);ylim([-0.08 0.08]);
    xlabel('k_{x = range} [rad/m]'); ylabel('k_{y = azimuth} [rad/m]'); title("E(k) [m^2 s/rad]"); colorbar;
    % 
    subplot(3,4,5)
    plot(interp_k_matrix(:,64));
    title("Interpolated ERA5 k"); xlabel('bin number'); ylabel('Magnitude'); 
    subplot(3,4,6);
    plot(interp_kx_matrix(64,:));
    title("Interpolated kx"); xlabel('bin number'); ylabel('Magnitude of kx'); 
    subplot(3,4,7);
    plot(interp_ky_matrix(:,64));
    title("Interpolated ky"); xlabel('bin number'); ylabel('Magnitude of ky'); 
    subplot(3,4,8); % Positive 2-D wave number spectrum E(kx,ky)=E(k)
    contour(interp_kx_matrix(1,:), interp_ky_matrix(:,1), interp_E_kx_ky,40);
    xlim([-0.05 0.05]);ylim([-0.08 0.08]);
    xlabel('k_{x = range} [rad/m]'); ylabel('k_{y = azimuth} [rad/m]'); title("Interpolated E(k) [m^2 s/rad]"); colorbar;

    subplot(3,4,9)
    contour(interp_k_matrix,40);
    title("Interpolated ERA5 k matrix"); xlabel('bin number'); ylabel('bin number'); colorbar;
    subplot(3,4,10);
    contour(interp_kx_matrix,40)
    title("Interpolated kx matrix"); xlabel('bin number'); ylabel('bin number'); colorbar;
    subplot(3,4,11);
    contour(interp_ky_matrix,40)
    title("Interpolated ky matrix"); xlabel('bin number'); ylabel('bin number'); colorbar;
    subplot(3,4,12); % Positive 2-D wave number spectrum E(kx,ky)=E(k)
    contour(interp_kx_matrix, interp_ky_matrix, interp_E_kx_ky,40);
    xlim([-0.05 0.05]);ylim([-0.08 0.08]);
    xlabel('k_{x = range} [rad/m]'); ylabel('k_{y = azimuth} [rad/m]'); title("Interpolated E(k) [m^2 s/rad]"); colorbar;

end
