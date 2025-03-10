classdef ERA5
%ERA5 Summary of this class goes here
    %   Detailed explanation goes here
    %
    % All the properties are set using the official ERA5 model
    % documentation: https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
    %   Specifically:
        % "
        % The NetCDF wave spectra file will have the dimensions longitude, latitude, direction, frequency and time.
        % However, the direction and frequency bins are simply given as 1 to 24 and 1 to 30, respectively.
        % The direction bins start at 7.5 degree and increase by 15 degrees until 352.5, with 90 degree being towards the east (Oceanographic convention).
        % The frequency bins are non-linearly spaced. The first bin is 0.03453 Hz and the following bins are: f(n) = f(n-1)*1.1; n=2,30. 
        % The data provided is the log10 of spectra density. To obtain the spectral density one has to take to the power 10 (10 ** data). This will give the units 2D wave spectra as m**2 s radian**-1 . Very small values are discarded and set as missing values. These are essentially 0 m**2 s radian**-1.
        % This recoding can be done with the Python xarray package, for example:
        % import xarray as xr
        % import numpy as np
        % da = xr.open_dataarray('2d_spectra_201601.nc')
        % da = da.assign_coords(direction=np.arange(7.5, 352.5 + 15, 15))
        % da = da.assign_coords(frequency=np.full(30, 0.03453) * (1.1 ** np.arange(0, 30)))
        % da = 10 ** da
        % da = da.fillna(0)
        % da.to_netcdf(path='2d_spectra_201601_recoded.nc')
        % Units of 2D wave spectra
        % Once decoded, the units of 2D wave spectra are m2 s radian-1
        % "

    properties
        FilePath
        DirectionBins
        FrequencyBins
    end

    methods
        function obj = ERA5(filepath)
            %ERA5 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.FilePath = filepath;

            % Decode the directions according to the ERA5 documentation
            obj.DirectionBins = ((7.5: 15: 352.5))'; % From ERA5 documentation
            
            % Decode the frequencies according to the ERA5 documentation
            era5_frequency = ncread(obj.FilePath,'frequencyNumber');
            obj.FrequencyBins = zeros(1,size(era5_frequency,1));
            obj.FrequencyBins(1) = 0.03453;
            for i=2:length(obj.FrequencyBins)
                obj.FrequencyBins(i) = obj.FrequencyBins(i-1)*1.1;
            end

        end

        function era5_d2fd_all_5dims = getAllWaveSpectrumD2FD(obj)
            % This method returns the ERA5 wave spectrum, E(f, theta) for
            % all dimentions [lon lat frq dir time]

            era5_d2fd_all_5dims = ncread(obj.FilePath,'d2fd'); % 'frequency direction spectrum' Wave energy spectrum OR variance density spectrum [m^2 s / rad] : how the wave energy is distributed across different spatial and temporal scales
            
            % Decode the spectra according to the ERA5 documentation
            era5_d2fd_all_5dims = 10.^era5_d2fd_all_5dims; % [m^2 s /rad]
            era5_d2fd_all_5dims(isnan(era5_d2fd_all_5dims)) = 0; % Necessary step according to the ERA5 documentation

        end

        function era5_lon = getLongitude(obj)
            % This method returns the longitude grid
            era5_lon = ncread(obj.FilePath,'longitude');

        end
        
        function era5_lat = getLatitude(obj)
            % This method returns the latitude grid
            era5_lat = ncread(obj.FilePath,'latitude');

        end

        function era5_time = getTime(obj)
            % This method returns the time values. The values are converted
            % from EPOCH time to real datetime using the ERA5 metadata.
            
            attribute_name = 'valid_time';
            % Read in the epoch time from the data
            era5_time_epoch = ncread(obj.FilePath,attribute_name);

            % Read in the EPOCH reference date - this can change
            % and should not be hardcoded
            variables_struct = ncinfo(obj.FilePath).Variables;
            for i = 1:length(variables_struct)
                if isequal(variables_struct(i).Name,attribute_name)
                    attributes = variables_struct(i).Attributes;
                    for j = 1: length(attributes)
                        if isequal(attributes(j).Name, 'units')
                            reference_date = ncinfo(obj.FilePath).Variables(i).Attributes(j).Value;

                            % Need to clean the value so that it is only
                            % the date-time string in YYYY-MM-DD format.
                            reference_date = reference_date(end-9:end);
                        end
                      
                    end
                   
                    
                end
            end
           
            % Convert from epoch time to real time values
            era5_time = datetime(int64(era5_time_epoch),'ConvertFrom','epochtime','Epoch',reference_date);

        end

        
        function [era5_d2fd, era5_lat,era5_lon, era5_time] = getSlicedWaveSpectrumD2Fd(obj,sar_center_latitude,sar_center_longitude,sar_time)
            % This method returns the ERA5 wave spectrum, E(f, theta)
            % sliced at the best match of the specified lat,lon and time
            % NOTE: lat and lon should be center of the sar transect.

            era5_d2fd_all_5dims = obj.getAllWaveSpectrumD2FD;
            era5_lon_all = obj.getLongitude;
            era5_lat_all = obj.getLatitude;
            
            [~,lat_match_index] = min(abs(era5_lat_all - sar_center_latitude));
            [~,lon_match_index] = min(abs(era5_lon_all - sar_center_longitude));

            era5_time_all = obj.getTime;

            % Find the slice that matches the SAR image
            era5_date_only = dateshift(era5_time_all, 'start', 'day');
            sar_date_only = dateshift(sar_time, 'start', 'day');
            time_match_index = era5_date_only == sar_date_only;
            time_match_index = find(hour(era5_time_all(time_match_index)) == hour(sar_time));

            % Slice the data
            era5_d2fd(:,:) = era5_d2fd_all_5dims(lon_match_index,lat_match_index,:,:,time_match_index);
            
            if(all(era5_d2fd == 0))
                era5_d2fd(:,:) = era5_d2fd_all_5dims(lat_match_index,lon_match_index,:,:,time_match_index);
            end

            % Transpose the data so that it is [dir frq] not [frq dir]
            era5_d2fd = era5_d2fd';

            % Return the latitudes and longitudes
            era5_lon = era5_lon_all(lon_match_index);
            era5_lat = era5_lat_all(lat_match_index);
            era5_time = era5_time_all(time_match_index);

        end

        function [date, longitude, latitude, mwd,mwp,swh, mdts,mpts,shts mdww,mpww,shww] = getSingleValues(obj,filepath_single,sar_center_latitude,sar_center_longitude,sar_time)
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
            
            
            % % Plot the Lat-Long x-y grid 
            % [lonGrid_single, latGrid_single] = meshgrid(lon_single,lat_single);
            % % Create a 2x2 grid of subplots
            % figure;
            % % WAVE DIRECTION
            % minVal = min([
            %     min(mean_direction_total_swell(:));
            %     min(mean_direction_wind_waves(:));
            %     min(mean_wave_direction(:));
            % ]);
            % 
            % % maxVal = max([
            % %     max(mean_direction_total_swell(:));
            % %     max(mean_direction_wind_waves(:));
            % %     max(mean_wave_direction(:));
            % % ]);
            % maxVal = 360;
            % 
            % subplot(3, 3, 1); % 2 rows, 2 columns, 1st subplot
            % pcolor(lonGrid_single, latGrid_single, mean_direction_total_swell(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Mean Direction - Total Swell');
            % 
            % subplot(3, 3, 2); % 2 rows, 2 columns, 1st subplot
            % pcolor(lonGrid_single, latGrid_single, mean_direction_wind_waves(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Mean Direction - Wind Waves');
            % 
            % subplot(3, 3, 3); % 2 rows, 2 columns, 1st subplot
            % pcolor(lonGrid_single, latGrid_single, mean_wave_direction(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Mean Direction - Wave');
            % 
            % % WAVE PERIOD
            % minVal = min([
            %     min(mean_period_total_swell(:));
            %     min(mean_period_wind_waves(:));
            %     min(mean_wave_period(:));
            % ]);
            % maxVal = max([
            %     max(mean_period_total_swell(:));
            %     max(mean_period_wind_waves(:));
            %     max(mean_wave_period(:));
            % ]); 
            % 
            % subplot(3, 3, 4); % 2 rows, 2 columns, 1st subplot
            % pcolor(lonGrid_single, latGrid_single, mean_period_total_swell(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Mean Period - Total Swell');
            % 
            % subplot(3, 3, 5); % 2 rows, 2 columns, 1st subplot
            % pcolor(lonGrid_single, latGrid_single, mean_period_wind_waves(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Mean Period - Wind Waves');
            % 
            % subplot(3, 3, 6); % 2 rows, 2 columns, 1st subplot
            % pcolor(lonGrid_single, latGrid_single, mean_wave_period(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Mean Period - Wave');
            % 
            % % SIGNIFICANT WAVE HEIGHT
            % minVal = min([ ...
            %     min(significant_wave_height_total_swell(:));
            %     min(significant_wave_height_wind_waves(:));
            %     min(significant_wave_height(:)) ...
            % ]);
            % maxVal = max([
            %     max(significant_wave_height_total_swell(:));
            %     max(significant_wave_height_wind_waves(:));
            %     max(significant_wave_height(:))
            % ]);
            % subplot(3, 3, 7); 
            % pcolor(lonGrid_single, latGrid_single, significant_wave_height_total_swell(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Significant Wave Height - Total Swell');
            % 
            % subplot(3, 3, 8); 
            % pcolor(lonGrid_single, latGrid_single, significant_wave_height_wind_waves(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Significant Wave Height - Wind Waves');
            % 
            % subplot(3, 3, 9);
            % pcolor(lonGrid_single, latGrid_single, significant_wave_height(:,:,time_index_single)'); % Note: transpose 'test' to match the grid orientation
            % clim([minVal maxVal]); colorbar; % Add a colorbar to show the data scale
            % xlabel('Longitude'); ylabel('Latitude'); title('Significant Wave Height - Wave');


        end
        
    end
end

%% Other stuff:
% era5_direction = ncread(filepath,'directionNumber');
% era5_frequency = ncread(filepath,'frequencyNumber');