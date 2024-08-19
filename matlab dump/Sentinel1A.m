classdef Sentinel1A
    %SENTINEL1A Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Filepath
        CaptureDate
        CaptureTime
        Polarisation
        Look
        SlantRange
        SatelliteVelocity
        Beta

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
            obj = obj.setAttributes(metadata_attributes);
        end

        function obj = setAttributes(obj,metadata_attributes)
            %SETATTRIBUTES Set object properties from metadata in Sentinel 1A data file
            %   metadata_attributes = the Attribute structure of the metadata MATLAB structure
                        
            % List of required attributes
            required_attributes = ["PROC_TIME", ...                     % Capture Date
                "first_line_time", "last_line_time", ...                % Capture Time
                "mds1_tx_rx_polar","mds2_tx_rx_polar",...               % Polarisation
                "antenna_pointing",...                                  % Look
                "slant_range_to_first_pixel", ...                       % Slant Range
                "Orbit_State_Vectors:orbit_vector1:time",...            % Sat Velocity start
                "SRGR_Coefficients:srgr_coef_list_1:zero_doppler_time"  % Sat Velocity end
                ];  

            % Update attribute list values to include prefix found in the
            % metadat_attributes structure
            prefix = 'Abstracted_Metadata:';
            required_attributes_prefix = prefix + required_attributes;

            % Find where in the structure the required attributes are 
            [~,matching_row_number] = ismember(required_attributes_prefix, {metadata_attributes.Name});

            % ERROR check: missing attributes in structure
            if any(matching_row_number == 0)
                error('Some required attributes are not found in metadata_attributes.');
            end

            % Capture Date
            format = 'dd-MMM-yyyy HH:mm:ss.SSSSSS';
            obj.CaptureDate = datetime(metadata_attributes(matching_row_number(1)).Value, 'InputFormat', format);

            % Capture Time
            format = 'dd-MMM-yyyy HH:mm:ss.SSSSSS';
            captureTime(1) = datetime(metadata_attributes(matching_row_number(2)).Value, 'InputFormat', format);
            captureTime(2) = datetime(metadata_attributes(matching_row_number(3)).Value, 'InputFormat', format);
            obj.CaptureTime = mean(captureTime);
            
            % Polarisation
            obj.Polarisation = [...
                metadata_attributes(matching_row_number(4)).Value,...
                metadata_attributes(matching_row_number(5)).Value];

            % Look
            obj.Look = metadata_attributes(matching_row_number(6)).Value;

            % Radar Slant Range
            obj.SlantRange = metadata_attributes(matching_row_number(7)).Value;
            
            % Satellite Velocity in m/s
            obj = obj.getSatelliteVelocity(metadata_attributes,matching_row_number);

            % Beta
            obj.Beta = obj.SlantRange/obj.SatelliteVelocity;
            
        end

        function obj = getSatelliteVelocity(obj,metadata_attributes,matching_row_number)
            %GETSATELLITEVELOCITY Summary of this method goes here
            %   Detailed explanation goes here
            
            % Calculate how many orbit vectors there are: each has 7
            % associated fields
            % num_orbit_state_vectors = (matching_row_number(9) - matching_row_number(8))/7;

            % Compute offsets for all indices
            orbit_vector_indices = (matching_row_number(8):7:matching_row_number(9)-6);  
            
            % Compute indices for time, x_vel, y_vel, and z_vel
            time_indices = orbit_vector_indices;
            x_vel_indices = 4 + orbit_vector_indices;
            y_vel_indices = 5 + orbit_vector_indices;
            z_vel_indices = 6 + orbit_vector_indices;
    
            format = 'dd-MMM-yyyy HH:mm:ss.SSSSSS';
            meta_orb_time = datetime([metadata_attributes(time_indices).Value], 'InputFormat', format);

            % Calculate time differences
            time_diff = abs(meta_orb_time - obj.CaptureTime);
            % Find the index of the closest datetime
            [~, time_index] = min(time_diff);

            meta_orb_x = [metadata_attributes(x_vel_indices(time_index)).Value];
            meta_orb_y = [metadata_attributes(y_vel_indices(time_index)).Value];
            meta_orb_z = [metadata_attributes(z_vel_indices(time_index)).Value];

            obj.SatelliteVelocity = sqrt(meta_orb_x.^2 + meta_orb_y.^2 + meta_orb_z.^2);
 
        end

        function sar_data = getSARdata(obj,polarisation)
            %GETSARDATA Summary of this method goes here
            %   Detailed explanation goes here
            [found, ~] = ismember(obj.Polarisation,polarisation);
            if(~found)
                error('Specified polarisation value "%s" is not valid OR does not exist in the object.',polarisation);
            else
                sar_data = ncread(obj.Filepath,'Sigma0_'+polarisation)';
            end
        end
        
        function thermal_calibrated_image = thermalNoiseCalibration(sar_image)
            %THERMALNOISECALIBRATION Preprocessing procedure
            %   Used to remove thermal noise (speckle) from image.
            
        end

        function [radiometric_calibrated_image] = radiometricCalibration(sar_image)
            %RADIOMETRICCALIBRATION Preprocessing procedure
            %   more detail....
            
        end

        function incidence_angle_grid_full = getIncidenceAngleGrid(obj)
            %GETINCIDENCEANGLEGRID Summary of this method goes here
            %   This returns the incidence angle of the full image
            incidence_angle_grid_full = ncread(obj.Filepath,'Incidence_Angle')';
        end

        function lat_grid = getLatitudeGrid(obj)
            %GETLATITUDEGRID Summary of this method goes here
            %   This returns the incidence angle of the full image
            lat_grid = ncread(obj.Filepath,'Lat')';
        end
       
        function lon_grid = getLongitudeGrid(obj)
            %GETLONGITUDEGRID Summary of this method goes here
            %   This returns the incidence angle of the full image
            lon_grid = ncread(obj.Filepath,'Lon')';
        end


    end
end