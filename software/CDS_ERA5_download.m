% Used to download data from ERA5 CDS
%% ==================== INITIALISATION ====================
clear; close all; clc;
year = '2025';
month = '03';
day = '07';
coords = [-30, 14, -35, 19]; % [-55, -10, -70, 10]

% Select which datatypes to download
downloadSingles = 1;
downloadFull = 1;
plotSingles = 1;
plotFull = 1;

% Setup parameters for the file naming
mission = "CapePoint"; % or 'SCALE'
date_str = year+""+month+""+day;

% For plotting
% latitude_of_interest = -58;
% longitude_of_interest = -5;
latitude_of_interest = -34.204;
longitude_of_interest = 18.28666944;
month_array = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
time_of_interest = datetime(day+"-"+month_array(str2num(month))+"-"+year+" 17:00:00"); %'13-Dec-2024 17:00:00'

disp('*==== Set Parameters ====*');
disp(time_of_interest);
disp("Chosen lat: " + num2str(latitude_of_interest));
disp("Chosen lon: " + num2str(longitude_of_interest));


%% ==================== Python Environment ====================
% 1. Follow the instructions in cds_setup.txt
% 2. Then restart matlab, uncomment the code below and run it. Make sure to 
% set the base path to where your virtual environment is located on your 
% computer.
% Note: If you have opted to not use virtual environment then use pyenv to 
% ensure that your MATLAB is running the correction version of python which 
% has the cds library installed.

% Path to your virtual environment's base directory
venv_base_path = '/Users/tris/Documents/MSc/';

% Command to activate the virtual environment and run a Python command
% Note: this should be the path to the env with csd installed
command = ['cd ' venv_base_path ' && source .venv/bin/activate && ' ...
    'python -c "import sys; print(sys.executable)"'];

% Execute the command
[status, cmdout] = system(command);

% Display the output to check the correct python path has been used 
disp(cmdout);
% Set the Python environment in MATLAB to the correct version
pyenv('Version', '/Users/tris/Documents/MSc/.venv/bin/python');
pyenv
%% ==================== CDS API ====================
% Now that the Python environment in MATLAB is correctly set, specify 
% the parameters that you want to download and run the code below.

% Check previous requests: https://cds.climate.copernicus.eu/cdsapp#!/yourrequests

% Import the cdsapi module
cdsapi = py.importlib.import_module('cdsapi');

% Create a Client instance
c = cdsapi.Client();
%% ==================== ERA5 Singles ====================
if downloadSingles
    params = py.dict(pyargs( ...
        'product_type', 'reanalysis', ...
        'variable', py.list({'mean_direction_of_total_swell', 'mean_direction_of_wind_waves', 'mean_period_of_total_swell', ...
                              'mean_period_of_wind_waves', 'mean_wave_direction', 'mean_wave_period', ...
                              'significant_height_of_combined_wind_waves_and_swell', 'significant_height_of_total_swell', ...
                              'significant_height_of_wind_waves', 'wave_spectral_kurtosis'}), ...
        'year', year, ...
        'month', month, ...
        'day', day, ...
        'time', py.list({'00:00', '01:00', '02:00', '03:00', '04:00', '05:00', ...
                         '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', ...
                         '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', ...
                         '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'}), ...
        'area', py.list(coords), ...
        'format', 'netcdf' ...
    ));
    
    % Define the output file name
    output_file_path_singles = mission+'_ERA5_single_2Dws_'+date_str+".nc";
    
    % Retrieve the values
    c.retrieve('reanalysis-era5-single-levels', params, output_file_path_singles);
    
    disp(['Single data downloaded:', output_file_path_singles]);
end

%% ==================== ERA5 Full 2D Spectrum ====================

if downloadFull
    full_date = year+"-"+month+"-"+day;

    % Define the parameters for the retrieve function
    params = py.dict(pyargs( ...
        'class', 'ea', ...
        'date', full_date, ...
        'direction', '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24', ...
        'domain', 'g', ...
        'expver', '1', ...
        'frequency', '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30', ...
        'param', '251.140', ...
        'stream', 'wave', ...
        'time', '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00', ...
        'type', 'an', ...
        'area', py.list(coords), ...  % North, West, South, East
        'grid', '1.0/1.0', ... % Latitude/longitude
        'format', 'netcdf' ... % Output format
        ));
    
    
    % Define the output file name
    output_file_path_full = mission+'_ERA5_2Dws_'+date_str+".nc";
    
    % Call the retrieve method
    c.retrieve('reanalysis-era5-complete', params, output_file_path_full);

    disp(['Full data downloaded:', output_file_path_full]);
end

%% ==================== Display outputs of downloads ====================
if plotFull
    disp('*==== Full Spectrum Outputs ====*');

    % Get Parameters of the complete data
    if downloadFull
        filepath = "/Users/tris/Documents/MSc/ice-SAR-pipeline/software/"+output_file_path_full;
    else
        filepath = "/Users/tris/Documents/MSc/data/ERA5-2Dws_20220720.nc";
    end

    % Store file info in MATLAB struct
    info = ncinfo(filepath);

    % Create an ERA5 object
    era5 = ERA5(filepath);
    
   % Get attributes
    longitudes = era5.getLongitude;
    latitudes = era5.getLatitude;
    time = era5.getTime;
    
    freq_bins = era5.FrequencyBins;
    direction_bins_degrees = era5.DirectionBins;
    
    [era5_d2fd, era5_lat, era5_lon, era5_datetime] = era5.getSlicedWaveSpectrumD2Fd(latitude_of_interest,longitude_of_interest,time_of_interest);

    % Plot the 2D spectra
    figure('Position', [100, 100, 800, 300]);
    % subplot(1,2,1);
    %     plotLibrary().waveSpectrum2D(0,era5_d2fd, freq_bins, ...
    %         direction_bins_degrees, 0.3, "2-D Wave Spectrum, E(f,theta)");
    % subplot(1,2,2);
        plotLibrary().waveSpectrum2D(1,era5_d2fd, freq_bins, ...
            direction_bins_degrees, 0.3, "2-D Wave Spectrum, E(f,theta)");

    
    % Calculate spectrum characteristics
    [Hs_init,Tm_init,direction_init,~] = ...
        waveLibrary().calculateSpectrumCharacteristics(direction_bins_degrees.*(pi/180),freq_bins,era5_d2fd);
    
    disp(era5_datetime);
    disp("Lat: " + num2str(era5_lat));
    disp("Lon: " + num2str(era5_lon));
    disp("  ");
    
    disp(['Hs [m]: ', num2str(Hs_init)]);
    disp(['Tm [m/s]: ', num2str(Tm_init)]);
    disp(['Dir [degrees]: ', num2str(direction_init)]);

end

if plotSingles
    disp('*==== Singles Outputs ====*');

    % Get Parameters of the complete data
    if downloadSingles
        filepath_single = "/Users/tris/Documents/MSc/ice-SAR-pipeline/software/"+output_file_path_singles;
    else
        filepath_single = "/Users/tris/Documents/MSc/data/ERA5-2Dws_20220720.nc";
    end

    % Store file info in MATLAB struct
    info_single = ncinfo(filepath_single);

    % Calculate spectrum characteristics
    if plotFull
        % want to ensure that the same latitude and longitude are chosen
        % (full spectrum does every 1 degree and singles does every 0.5
        % degree)
        [era5_datetime_single, era5_lon_single, era5_lat_single, ...
        era5_direction,era5_Tm,era5_Hs, era5_mdts,era5_mpts,era5_shts, ...
        era5_mdww, era5_mpww, era5_shww] = ...
            getSingleValues(filepath_single, ...
                era5_lat, ...
                era5_lon, ...
                time_of_interest);
    else
        [era5_datetime_single, era5_lon_single, era5_lat_single, ...
        era5_direction,era5_Tm,era5_Hs, era5_mdts,era5_mpts,era5_shts, ...
        era5_mdww, era5_mpww, era5_shww] = ...
            getSingleValues(filepath_single, ...
                latitude_of_interest, ...
                longitude_of_interest, ...
                time_of_interest);
    end
   
    disp(era5_datetime_single);
    disp("Lat: " + num2str(era5_lat_single));
    disp("Lon: " + num2str(era5_lon_single));
    disp("  ");
    disp(['Singles ERA5 Hs: ', num2str(era5_Hs)]);
    disp(['Singles ERA5 Tm: ', num2str(era5_Tm)]);
    disp(['Singles ERA5 Dir: ', num2str(era5_direction)]);
    disp("  ");
    disp(['Singles ERA5 mdts: ', num2str(era5_mdts)]);
    disp(['Singles ERA5 mpts: ', num2str(era5_mpts)]);
    disp(['Singles ERA5 shts: ', num2str(era5_shts)]);
    disp("  ");
    disp(['Singles ERA5 mdww: ', num2str(era5_mdww)]);
    disp(['Singles ERA5 mpww: ', num2str(era5_mpww)]);
    disp(['Singles ERA5 shww: ', num2str(era5_shww)]);


end

%% ==================== FUNCTIONS ====================

function [date, longitude, latitude, mwd,mwp,swh, mdts,mpts,shts mdww,mpww,shww] = getSingleValues(filepath_single,sar_center_latitude,sar_center_longitude,sar_time)

    % This method returns the wave spectrum
    % characteristics/integral values calculated by the EMCWF and
    % downloaded from CDS. Specifically, matching the lat, long and
    % time values of the SAR image.

    %    filepath_singles must be the datafile containing the
    %           single values

    info_single = ncinfo(filepath_single);

    %data storing 2D wave spectra with 24 directions (starting from 7.5N
    %direction towards) and 30 frequencie starting from 0.0345 Hz;
    era5_lon_single = ncread(filepath_single,'longitude');
    era5_lat_single = ncread(filepath_single,'latitude');
    era5_time_epoch_single = ncread(filepath_single,'valid_time');
    
    % find matches to the SAR transect
    [~,lat_match_index] = min(abs(era5_lat_single - sar_center_latitude));
    [~,lon_match_index] = min(abs(era5_lon_single - sar_center_longitude));
    % lon_match_index = lon_match_index-1;
    longitude = era5_lon_single(lon_match_index);
    latitude = era5_lat_single(lat_match_index);

    % Find the slice that matches the SAR image
    reference_date = '1970-01-01';
    era5_time_single = datetime(int64(era5_time_epoch_single),'ConvertFrom','epochtime','Epoch',reference_date);
    era5_date_only = dateshift(era5_time_single, 'start', 'day');
    sar_date_only = dateshift(sar_time, 'start', 'day');
    time_match_index = era5_date_only == sar_date_only;
    time_match_index = find(hour(era5_time_single(time_match_index)) == hour(sar_time));
    date = era5_time_single(time_match_index);

    % Mean Wave Direction (θ): The average direction from which the swell waves are coming.
    % Mean Wave Period (T): The average period of the swell waves.
    mean_wave_direction = ncread(filepath_single,'mwd');
    mwd = mean_wave_direction(lat_match_index, lon_match_index, time_match_index); % [lon lat time]
    
    mean_wave_period = ncread(filepath_single,'mwp');
    mwp = mean_wave_period(lat_match_index, lon_match_index, time_match_index); % [lon lat time]
    
    significant_wave_height = ncread(filepath_single,'swh');
    swh = significant_wave_height(lat_match_index, lon_match_index, time_match_index); % [lon lat time]

    mean_direction_total_swell = ncread(filepath_single,'mdts'); %121 31 48 = lon lat time
    mdts = mean_direction_total_swell(lat_match_index, lon_match_index, time_match_index); % [lon lat time]

    mean_period_total_swell = ncread(filepath_single,'mpts');
    mpts = mean_period_total_swell(lat_match_index, lon_match_index, time_match_index); % [lon lat time]

    significant_wave_height_total_swell = ncread(filepath_single,'shts');
    shts = significant_wave_height_total_swell(lat_match_index, lon_match_index, time_match_index); % [lon lat time]

    mean_direction_wind_waves = ncread(filepath_single,'mdww');
    mdww = mean_direction_wind_waves(lat_match_index, lon_match_index, time_match_index); % [lon lat time]

    mean_period_wind_waves = ncread(filepath_single,'mpww');
    mpww = mean_period_wind_waves(lat_match_index, lon_match_index, time_match_index); % [lon lat time]

    significant_wave_height_wind_waves = ncread(filepath_single,'shww');
    shww = significant_wave_height_wind_waves(lat_match_index, lon_match_index, time_match_index); % [lon lat time]
    
    
    % Plot the Lat-Long x-y grid 
    [lonGrid_single, latGrid_single] = meshgrid(era5_lon_single,era5_lat_single);
    % Create a 2x2 grid of subplots
    figure;
    % WAVE DIRECTION
    minVal = min([
        min(mean_direction_total_swell(:));
        min(mean_direction_wind_waves(:));
        min(mean_wave_direction(:));
    ]);

    % maxVal = max([
    %     max(mean_direction_total_swell(:));
    %     max(mean_direction_wind_waves(:));
    %     max(mean_wave_direction(:));
    % ]);
    maxVal = 360;

    subplot(3, 3, 1); % 2 rows, 2 columns, 1st subplot
    pcolor(lonGrid_single, latGrid_single, mean_direction_total_swell(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Mean Direction - Total Swell');

    subplot(3, 3, 2); % 2 rows, 2 columns, 1st subplot
    pcolor(lonGrid_single, latGrid_single, mean_direction_wind_waves(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Mean Direction - Wind Waves');

    subplot(3, 3, 3); % 2 rows, 2 columns, 1st subplot
    pcolor(lonGrid_single, latGrid_single, mean_wave_direction(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Mean Direction - Wave');

    % WAVE PERIOD
    minVal = min([
        min(mean_period_total_swell(:));
        min(mean_period_wind_waves(:));
        min(mean_wave_period(:));
    ]);
    maxVal = max([
        max(mean_period_total_swell(:));
        max(mean_period_wind_waves(:));
        max(mean_wave_period(:));
    ]); 

    subplot(3, 3, 4); % 2 rows, 2 columns, 1st subplot
    pcolor(lonGrid_single, latGrid_single, mean_period_total_swell(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Mean Period - Total Swell');

    subplot(3, 3, 5); % 2 rows, 2 columns, 1st subplot
    pcolor(lonGrid_single, latGrid_single, mean_period_wind_waves(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Mean Period - Wind Waves');

    subplot(3, 3, 6); % 2 rows, 2 columns, 1st subplot
    pcolor(lonGrid_single, latGrid_single, mean_wave_period(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Mean Period - Wave');

    % SIGNIFICANT WAVE HEIGHT
    minVal = min([ ...
        min(significant_wave_height_total_swell(:));
        min(significant_wave_height_wind_waves(:));
        min(significant_wave_height(:)) ...
    ]);
    maxVal = max([
        max(significant_wave_height_total_swell(:));
        max(significant_wave_height_wind_waves(:));
        max(significant_wave_height(:))
    ]);
    subplot(3, 3, 7); 
    pcolor(lonGrid_single, latGrid_single, significant_wave_height_total_swell(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Significant Wave Height - Total Swell');

    subplot(3, 3, 8); 
    pcolor(lonGrid_single, latGrid_single, significant_wave_height_wind_waves(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Significant Wave Height - Wind Waves');

    subplot(3, 3, 9);
    pcolor(lonGrid_single, latGrid_single, significant_wave_height(:,:,time_match_index)'); % Note: transpose 'test' to match the grid orientation
    clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
    xlabel('Longitude'); ylabel('Latitude'); title('Significant Wave Height - Wave');

end


% % Create a Client instance
% c = cdsapi.Client();
% 
% params = py.dict(pyargs( ...
%     'product_type', 'reanalysis', ...
%     'variable', py.list({'mean_direction_of_total_swell', 'mean_direction_of_wind_waves', 'mean_period_of_total_swell', ...
%                          'mean_period_of_wind_waves', 'mean_wave_direction', 'mean_wave_period', ...
%                          'significant_height_of_combined_wind_waves_and_swell', 'significant_height_of_total_swell', ...
%                          'significant_height_of_wind_waves', 'wave_spectral_kurtosis'}), ...
%     'year', '2022', ...
%     'month', '07', ...
%     'day', '22', ...
%     'time', py.list({'00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', ...
%                      '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', ...
%                      '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'}), ...
%     'area', py.list([-55, -10, -70, 10]), ...
%     'format', 'netcdf' ...
%     ));
% 
% % Define the output file name
% date = "20220722";
% output_file = 'ERA5-single-2Dws_'+date+".nc";
% 
% % c.retrieve('SCALE_reanalysis-era5-single-levels', params, output_file);
% 

% % Import the cdsapi module
% cdsapi = py.importlib.import_module('cdsapi');
% 
% % Create a Client instance
% c = cdsapi.Client();
% 
% params = py.dict(pyargs( ...
%     'product_type', 'reanalysis', ...
%     'variable', py.list({'mean_direction_of_total_swell', 'mean_direction_of_wind_waves', 'mean_period_of_total_swell', ...
%                            'mean_period_of_wind_waves', 'mean_wave_direction', 'mean_wave_period', ...
%                            'significant_height_of_combined_wind_waves_and_swell', 'significant_height_of_total_swell', ...
%                            'significant_height_of_wind_waves', 'wave_spectral_kurtosis'}), ...
%     'year', '2022', ...
%     'month', '07', ...
%     'day', '22', ...
%     'time', py.list({'00:00', '01:00', '02:00', '03:00', '04:00', '05:00', ...
%                      '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', ...
%                      '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', ...
%                      '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'}), ...
%     'area', py.list([-55, -10, -70, 10]), ...
%     'format', 'netcdf' ...
% ));
% 
% 
% % Define the output file name
% date = "20220722_2";
% output_file = 'SCALE_ERA5-single-2Dws_'+date+".nc";
% 
% c.retrieve('reanalysis-era5-single-levels', params, output_file);

% % Import the cdsapi module
% cdsapi = py.importlib.import_module('cdsapi');
% 
% % Create a Client instance
% c = cdsapi.Client();
% 
% % Define the parameters for the retrieve function
% params = py.dict(pyargs( ...
%     'class', 'ea', ...
%     'date', '2022-07-22', ...
%     'direction', '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24', ...
%     'domain', 'g', ...
%     'expver', '1', ...
%     'frequency', '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30', ...
%     'param', '251.140', ...
%     'stream', 'wave', ...
%     'time', '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00', ...
%     'type', 'an', ...
%     'area', py.list([-55, -10, -70, 10]), ...  % North, West, South, East
%     'grid', '1.0/1.0', ... % Latitude/longitude
%     'format', 'netcdf' ... % Output format
%     ));
% 
% 
% % Define the output file name
% date = "20220722";
% output_file = 'SCALE_ERA5-2Dws_'+date+".nc";
% 
% % Call the retrieve method
% c.retrieve('reanalysis-era5-complete', params, output_file);