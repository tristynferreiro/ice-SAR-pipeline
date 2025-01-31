classdef Sentinel1A
    %SENTINEL1A Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Filepath
        ProductType
        AcquisitionMode
        OrbitalPass             %Ascending/Descending
        PlatformHeading         %Azimuthal heading relative to North from S1A original data

        %%%%% Properties for inversion
        AcquisitionStartDatetime
        AcquisitionStopDatetime
        Polarisation
        LookDirection
        AzimuthResolution
        RangeResolution
        % % SlantRangeToFirstPixel
        % SlantRangeToTransect
        % IncidenceAngleTransect
        SatelliteVelocity           % [m/s] 
        SceneOrientationAngle       % [degrees to North] Azimuthal Angle / Platform Heading / Scene Orientation from SNAP export
        TransectSlantRanges
        TransectIncidenceAngles
        
        

    end

    methods
        function obj = Sentinel1A(filepath)
            %SENTINEL1A Construct an instance of this class
            %   INPUT PARAMS: sentinel1A_filepath = path to the .nc file
            
            obj.Filepath = filepath;
            
            % Extract the attributes table from the metadata MATLAB
            % structure
            metadata = ncinfo(filepath,'metadata');
            metadata_attributes = metadata.Attributes;
            
            % Set the object properties using the metadata
            % List of required attributes
            required_attributes = ["PROC_TIME", ...                         % Capture Date
                "first_line_time", "last_line_time", ...                    % Capture Time
                "mds1_tx_rx_polar","mds2_tx_rx_polar",...                   % Polarisation
                "antenna_pointing",...                                      % Look
                "slant_range_to_first_pixel", ...                           % Slant Range
                "Orbit_State_Vectors:orbit_vector1:time",...                % Sat Velocity start
                "Orbit_State_Vectors:orbit_vector17:z_vel",...              % Sat Velocity end
                "PRODUCT_TYPE",...
                "ACQUISITION_MODE",...
                "PASS",...
                "range_spacing","azimuth_spacing",...
                "centre_heading","centre_heading2"                          % Azimuthal Angle / Platform Heading / Scene Orientation; center heading is angle to North (0 or 360) and center heading 2 is angle to South (180)
                ];  

            % Update attribute list values to include prefix found in the
            % metadat_attributes structure
            prefix = 'Abstracted_Metadata:';
            required_attributes_prefix = prefix + required_attributes;

            % Find where in the structure the required attributes are 
            [~,matching_row_number] = ismember(required_attributes_prefix, {metadata_attributes.Name});

            % ERROR check: missing attributes in structure
            if any(matching_row_number == 0)
                % Update attribute list values to include prefix found in the
                % metadat_attributes structure
                prefix = 'Metadata_Group:Abstracted_Metadata:'; 
                required_attributes_prefix = prefix + required_attributes;

                % Find where in the structure the required attributes are 
                [~,matching_row_number] = ismember(required_attributes_prefix, {metadata_attributes.Name});

                if any(matching_row_number == 0)
                    error('Some required attributes are not found in metadata_attributes.');
                end
            end

            % Capture Date
            obj.AcquisitionStartDatetime = datetime(metadata_attributes(matching_row_number(2)).Value, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSSSSS');
            obj.AcquisitionStopDatetime = datetime(metadata_attributes(matching_row_number(3)).Value, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSSSSS');
            
            % Capture Time
            % obj.CaptureTime = mean([ obj.AcquisitionStartDatetime, obj.AcquisitionStopDatetime]);
            
            % Polarisation
            % obj.Polarisation = [...
            %     metadata_attributes(matching_row_number(4)).Value,...
            %     metadata_attributes(matching_row_number(5)).Value];
            obj.Polarisation = [...
                metadata_attributes(matching_row_number(5)).Value];

            % Look
            obj.LookDirection = metadata_attributes(matching_row_number(6)).Value;

            % Radar Slant Range to first pixel
            % obj.SlantRangeToFirstPixel = metadata_attributes(matching_row_number(7)).Value;
            
            % Product Type
            obj.ProductType  = metadata_attributes(matching_row_number(10)).Value;
            
            % Acquisition Mode
            obj.AcquisitionMode  = metadata_attributes(matching_row_number(11)).Value;

            % Pass direction 
            obj.OrbitalPass  = metadata_attributes(matching_row_number(12)).Value;
            
            % Range Resolution
            obj.RangeResolution  = metadata_attributes(matching_row_number(13)).Value;

            % Azimuth Resolution
            obj.AzimuthResolution  = metadata_attributes(matching_row_number(14)).Value;

            % Scene Orientation Angle / Azimuthal Angle / Platform Heading 
            % obj.SceneOrientationAngle  = metadata_attributes(matching_row_number(15)).Value;
            obj.SceneOrientationAngle  = metadata_attributes(matching_row_number(15)).Value;
            
            name = 'platformHeading'; % 'Original_Product_Metadata:annotation:s1a-iw-grd-vh-20241014t173427-20241014t173452-056101-06dd48-002_xml:product:generalAnnotation:productInformation:platformHeading'
            % Find where in the structure the required attributes are 
            [matching_index] = find(endsWith({metadata_attributes.Name},name) ==1, 1); % Find first match
            obj.PlatformHeading = str2num(metadata_attributes(matching_index(1)).Value);

            % Satellite Velocity in m/s
            % Compute offsets for all indices
            orbit_vector_indices = (matching_row_number(8):7:matching_row_number(9)-5);  
            
            % Compute indices for time, x_vel, y_vel, and z_vel
            time_indices = orbit_vector_indices;
            x_vel_indices = 4 + orbit_vector_indices;
            y_vel_indices = 5 + orbit_vector_indices;
            z_vel_indices = 6 + orbit_vector_indices;
    
            
            meta_orb_x = [metadata_attributes(x_vel_indices).Value];
            meta_orb_y = [metadata_attributes(y_vel_indices).Value];
            meta_orb_z = [metadata_attributes(z_vel_indices).Value];

            obj.SatelliteVelocity = mean(sqrt(meta_orb_x.^2 + meta_orb_y.^2 + meta_orb_z.^2));

            % ****************************************************
            % NOTE: accuracy can be improved by only averaging the state
            % vector velocities for the time periods that align with the
            % transects used.

            % format = 'dd-MMM-yyyy HH:mm:ss.SSSSSS';
            % meta_orb_time = datetime([metadata_attributes(time_indices).Value], 'InputFormat', format);
            % 
            % % Calculate time differences
            % mean_capture_time = mean([ obj.AcquisitionStartDatetime, obj.AcquisitionStopDatetime]);
            % time_diff = abs(meta_orb_time - mean_capture_time);
            % % Find the index of the closest datetime
            % [~, time_index] = min(time_diff);

            % meta_orb_x = [metadata_attributes(x_vel_indices(time_index)).Value];
            % meta_orb_y = [metadata_attributes(y_vel_indices(time_index)).Value];
            % meta_orb_z =
            % [metadata_attributes(z_vel_indices(time_index)).Value];
            % ****************************************************
            obj.TransectSlantRanges = table([],[], 'VariableNames', {'Field', 'Value'});
            obj.TransectIncidenceAngles = table([],[], 'VariableNames', {'Field', 'Value'});

        end

        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % function obj = getSatelliteVelocity(obj)
        %     %GETSATELLITEVELOCITY Summary of this method goes here
        %     %   Detailed explanation goes here
        % 
        % end
 
        function sar_data = getCalibratedSARImage(obj,polarisation)
            %GETCALIBRATEDSARImage returns calibrated Sentinel 1A SAR image
            %   Imports the SAR Image that has been calibrated using
            %   SNAP IDE.
            
            [found, ~] = ismember(obj.Polarisation,polarisation);
            if(~found)
                error('Specified polarisation value "%s" is not valid OR does not exist in the object.',polarisation);
            else
                sar_data = ncread(obj.Filepath,'Sigma0_'+polarisation)';
            end
        end
        
 

        function lat_grid = getLatitudeGrid(obj)
            %GETLATITUDEGRID 
            %   This returns the latitude of the full image
            lat_grid = ncread(obj.Filepath,'lat')';
        end
       
        function lon_grid = getLongitudeGrid(obj)
            %GETLONGITUDEGRID 
            %   This returns the longitude of the full image
            lon_grid = ncread(obj.Filepath,'lon')';
        end

        function [slant_range_to_transect_center] = getSlantRange(obj,sar_transect_lat_indices,sar_transect_lon_indices)
            c = physconst('LightSpeed');
            RTT = c /2;
            slant_range_time = ncread(obj.Filepath,"slant_range_t") .* 10e-9;
            slant_range_time = slant_range_time(sar_transect_lat_indices,sar_transect_lon_indices); % center of the transect
            % sar_slant_range = sar_slant_range_time(sar_transect_size/2,sar_transect_size/2) .* RTT;

            slant_range_to_transect_center = mean(slant_range_time(:)) .* RTT;
        end

        function [incidence_angle_at_transect] = getIncidenceAngle(obj,sar_transect_lat_indices,sar_transect_lon_indices)
            %
            %   sar_transect_lat_indices = all the latitude indices (should
            %   match what is used to slice the sar data)
            %   sar_transect_lon_indices = all the longitude indices (should
            %   match what is used to slice the sar data)
            
            incident_angle_grid = ncread(obj.Filepath,"i_angle");

            incidence_angle_degrees_transect = incident_angle_grid(sar_transect_lat_indices,sar_transect_lon_indices); % center of the transect
            % sar_center_incidence_angle_degrees =
            % sar_center_incidence_angle_degrees_transect(sar_transect_size/2,sar_transect_size/2);
            % % get center pixel
            % sar_center_incidence_angle_degrees =
            % sar_center_incidence_angle_degrees *
            % ones(sar_sub_transect_size); % greate matrix

            incidence_angle_at_transect = mean(incidence_angle_degrees_transect(:));
        end
        
        %------------------------------------------------------------------
        % function thermal_calibrated_image = thermalNoiseCalibration(sar_image)
        %     %THERMALNOISECALIBRATION Preprocessing procedure
        %     %   Used to remove thermal noise (speckle) from image.
        % 
        % end

        function azimuth_to_north_angle = azimuthToNorthAngleConversion(obj)
            % this is the angle between the azimuth plane of the SAR image
            % and the true north direction.
            if obj.LookDirection == "right" && obj.OrbitalPass == 'ASCENDING'
                azimuth_to_north_angle = -1 * (360 - obj.SceneOrientationAngle) ; % degrees
                % The -1 is applied because the SAR geometry means that the
                % angle is a counterclockwise rotation.
            end
        end

    end
end