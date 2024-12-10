classdef ERA5
    %ERA5 Summary of this class goes here
    %   Detailed explanation goes here
    %
    % All the properties are set using the official ERA5 model
    % documentation: https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation

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

        
        function era5_direction_bins_adjusted = getDirectionsInSARGeometry(obj, sar_azimuth_to_north_angle)
            % We need to rotate the wave spectrum so that it is aligned with the angle
            % at which the SAR data has been taken.
            %   NOTE: the angle needs to be the azimuth to north angle not
            %   the raw scene orientation value from the Metadata (i.e. 180
            %   - scene orientation angle OR 180+scene orientation angle)

            era5_direction_bins_adjusted = obj.DirectionBins + sar_azimuth_to_north_angle;
        end
        
        function [era5_d2fd, lat_match_index,lon_match_index, time_match_index] = getSlicedWaveSpectrumD2Fd(obj,sar_center_latitude,sar_center_longitude,sar_time)
            % This method returns the ERA5 wave spectrum, E(f, theta)
            % sliced at the best match of the specified lat,lon and time
            % NOTE: lat and lon should be center of the sar transect.

            era5_d2fd_all_5dims = obj.getAllWaveSpectrumD2FD;
            era5_lon = obj.getLongitude;
            era5_lat = obj.getLatitude;
            
            lat_match_index = find(era5_lat == ceil(sar_center_latitude));
            lon_match_index = find(era5_lon == ceil(sar_center_longitude));

            
            era5_time = obj.getTime;

            % Find the slice that matches the SAR image
            era5_date_only = dateshift(era5_time, 'start', 'day');
            sar_date_only = dateshift(sar_time, 'start', 'day');
            time_match_index = era5_date_only == sar_date_only;
            time_match_index = find(hour(era5_time(time_match_index)) == hour(sar_time));

            % Slice the data
            era5_d2fd(:,:) = era5_d2fd_all_5dims(lon_match_index,lat_match_index,:,:,time_match_index);
            
            if(all(era5_d2fd == 0))
                era5_d2fd(:,:) = era5_d2fd_all_5dims(lat_match_index,lon_match_index,:,:,time_match_index);
            end
            % Transpose the data so that it is [dir frq] not [frq dir]
            era5_d2fd = era5_d2fd';

        end

        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

%% Other stuff:
% era5_direction = ncread(filepath,'directionNumber');
% era5_frequency = ncread(filepath,'frequencyNumber');