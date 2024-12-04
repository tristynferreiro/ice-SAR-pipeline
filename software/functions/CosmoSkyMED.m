classdef CosmoSkyMED
    %COSMOSKYMED Summary of this class goes here
    %   Detailed explanation goes here
    %  
    %   ----------------------------------------------------------
    %   In accordance with the e-geos documentation:
    %       acquisition_mode = A_mode: SCS_B, DGM_B, GEC_B and GTC_B images
    %       range_spreading_loss_compensation_flag = Rsl_flag
    %       ref_slant_range = R_ref in meters
    %       ref_slant_range_exponent = R_exp = used to do the range
    %               spreading loss correction
    %       incidence_angle_compensation_flag = Inc_flag
    %       ref_incidence_angle = α_ref in degrees
    %       rescaling_factor = F
    %       calibration_constant_compensation_flag = K_flag = 
    %               Calibration Constant Compensation Flag
    %       calibration_constant = K
    %       calibrated_image = σ^0

    properties
        Filepath
        ProductType
        AcquisitionMode
        %OrbitalPass              %Ascending/Descending 

        %%%%% Properties for inversion
        AcquisitionStartDatetime
        AcquisitionStopDatetime
        Polarisation
        LookDirection
        RangeResolution
        AzimuthResolution
        SatelliteVelocity = 7541.89; %m/s from Giacomo's output file.
        SceneOrientationAngle
        
        
        %%%%% Properties for calibration
        RefSlantRange
        RefSlantRangeExponent
        RefIncidenceAngle
        CalibrationConstant
        RescalingFactor
        RangeSpreadingLossCompensationFlag
        IncidenceAngleCompensationFlag
        CalibrationConstantCompensationFlag
    end

    methods
        function obj = CosmoSkyMED(filepath,product_type)
            %COSMOSKYMED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Filepath = filepath;

            obj.ProductType = product_type;
            
            required_attributes = ["ACQUISITION_MODE",...   % Acquisition mode
                      "ref_slant_range", ...
                      "ref_slant_range_exp",...
                      "ref_inc_angle",...
                      "calibration_factor",...
                      "rescaling_factor",...
                      "range_spread_comp_flag",...          % Range Spreading Loss Compensation Flag
                      "inc_angle_comp_flag",...
                      "abs_calibration_flag",...
                      "mds1_tx_rx_polar",...
                      "antenna_pointing",...
                      "range_spacing","azimuth_spacing",...
                      "slant_range_to_first_pixel"...
                ];  

            if product_type=="CSG"
                % List of required attributes
                required_attributes = "Abstracted_Metadata:"+required_attributes;
                metadata_attributes = ncinfo(filepath,'metadata').Attributes;
            elseif product_type=="CSK"
                metadata_attributes = ncinfo(filepath,'Metadata_Group').Groups(1).Attributes; 
            end
            
            % Find where in the structure the required attributes are 
            [~,matching_row_number] = ismember(required_attributes, {metadata_attributes.Name});

            % ERROR check: missing attributes in structure
            if any(matching_row_number == 0)
                error('Some required attributes are not found in metadata_attributes.');
            end

            % Acquisition mode A_mode
            obj.AcquisitionMode = metadata_attributes(matching_row_number(1)).Value;

            % Reference Slant Range
            obj.RefSlantRange = metadata_attributes(matching_row_number(2)).Value;

            % Reference Slant Range exponent
            obj.RefSlantRangeExponent = metadata_attributes(matching_row_number(3)).Value;

            % Reference Incidence Angle
            obj.RefIncidenceAngle = metadata_attributes(matching_row_number(4)).Value;

            % Calibration Constant
            obj.CalibrationConstant = metadata_attributes(matching_row_number(5)).Value;

            % Rescaling Factor
            obj.RescalingFactor = metadata_attributes(matching_row_number(6)).Value;

            % Range Spreading Loss Compensation Geometry Flag - Rsl_flag
            obj.RangeSpreadingLossCompensationFlag = metadata_attributes(matching_row_number(7)).Value;
            
            % Incidence Angle Compensation Geometry Flag - Inc_flag
            obj.IncidenceAngleCompensationFlag = metadata_attributes(matching_row_number(8)).Value;

            % Calibration Constant Compensation Flag - K_flag
            obj.CalibrationConstantCompensationFlag = metadata_attributes(matching_row_number(9)).Value;
            
            % Polarisation
            obj.Polarisation = metadata_attributes(matching_row_number(10)).Value;

            % Look Direction (antenna pointing direction)
            obj.LookDirection = metadata_attributes(matching_row_number(11)).Value;

            % Resolution
            obj.RangeResolution = metadata_attributes(matching_row_number(12)).Value;
            obj.AzimuthResolution = metadata_attributes(matching_row_number(13)).Value;
            
            % Slant Range to First Pixel
            % obj.SlantRangeToFirstPixel = metadata_attributes(matching_row_number(14)).Value;

            % Acquisition times
            [obj.AcquisitionStartDatetime, obj.AcquisitionStopDatetime] = obj.getAcquisitiontime;

            % Scene Orientation Angle / Azimuthal Angle to True north
            obj.SceneOrientationAngle = obj.getSceneOrientationAngle;
   
        end
        
        %------------------------------------------------------------------
        

        
        %------------------------------------------------------------------
        function [sar_data] = getSARImageAmplitude(obj)
            %GETImage Get the sar data from file
            sar_data = ncread(obj.Filepath,'Amplitude'); 
        end

        function [radiometric_calibrated_image] = getRadiometricallyCalibratedImage(obj)
            %RADIOMETRICCALIBRATION Preprocessing procedure 
            %   Uses the single pixel backscattering coefficient, σ^0 
            %   In accordance with the e-geos documentation
            
            % Step 1: Evaluate the power image (intensity)
            power_image = abs(obj.getSARImageAmplitude).^2;
        
            % Step 2:
            scaling_factor = 0;
            if(~isequal(obj.RangeSpreadingLossCompensationFlag, "NONE"))
                % Remove the reference slant angle
                scaling_factor = obj.RefSlantRange.^(2.*obj.RefSlantRangeExponent);
            end
        
            % Step 3:
            if(~isequal(obj.IncidenceAngleCompensationFlag,"NONE"))
                % Remove the reference incidence angle
                scaling_factor = scaling_factor.*sind(obj.RefIncidenceAngle);
            end
        
            % Step 4: Remove the rescaling factor
            scaling_factor = scaling_factor.* (1/obj.RescalingFactor.^2);
        
            % Step 5: 
            if(obj.CalibrationConstantCompensationFlag==0)
                % Apply the calibration factor 
                scaling_factor = scaling_factor.*(1/obj.CalibrationConstant);
            end
        
            % Step 6: Apply the total scaling factor
            radiometric_calibrated_image = power_image.*scaling_factor;
               
        end

        function [acquisition_start_datetime, acquisition_stop_datetime] = getAcquisitiontime(obj)
            % Extract the acquisition times as a datetime object
            
            % Names of the required attributes
            required_attributes = ["start_date","stop_date"];

            if obj.ProductType=="CSG"
                % List of required attributes
                % required_attributes = "Abstracted_Metadata:"+required_attributes;
                % metadata_attributes = ncinfo(filepath,'metadata').Attributes;
            elseif obj.ProductType=="CSK"
                metadata_attributes = ncinfo(obj.Filepath).Attributes; 
            end
            
            % Find where in the structure the required attributes are 
            [~,matching_row_number] = ismember(required_attributes, {metadata_attributes.Name});

            % Return the times
            acquisition_information = metadata_attributes(matching_row_number(1)).Value;
            acquisition_start_datetime = datetime(acquisition_information, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSSSSS');
            acquisition_information = metadata_attributes(matching_row_number(2)).Value;
            acquisition_stop_datetime = datetime(acquisition_information, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSSSSS');
            
        end

        function [scene_orientation_angle] = getSceneOrientationAngle(obj)
            % Extract the acquisition times as a datetime object
            
            % Names of the required attributes
            required_attributes = "Scene_Orientation";

            if obj.ProductType=="CSG"
                % List of required attributes
                % required_attributes = "Abstracted_Metadata:"+required_attributes;
                % metadata_attributes = ncinfo(filepath,'metadata').Attributes;
            elseif obj.ProductType=="CSK"
                metadata_attributes = ncinfo(obj.Filepath,'Metadata_Group').Groups(2).Groups(1).Attributes; 
            end
            
            [~,matching_row_number] = ismember(required_attributes, {metadata_attributes.Name});
            scene_orientation_angle = metadata_attributes(matching_row_number(1)).Value; % [degrees]
            % IN Sentinel this is called the "centre_heading" and "centre_heading2"  
        end
        
        function [lat_grid] = getLatitudeGrid(obj)
            %GETLATITUDEGRID 
            %   This returns the latitude of the full image
            lat_grid = ncread(obj.Filepath,'lat')';
        end
       
        function [lon_grid] = getLongitudeGrid(obj)
            %GETLONGITUDEGRID 
            %   This returns the longitude of the full image
            lon_grid = ncread(obj.Filepath,'lon')';
        end
        
        function [slant_range_to_transect_center] = getSlantRange(obj,sar_transect_lat_end_index,sar_transect_lat_start_index,sar_transect_lon_start_index,sar_transect_lon_end_index)
            c = physconst('LightSpeed');
            RTT = c /2;
            slant_range_time = ncread(obj.Filepath,"slant_range_time") .* 10e-9;
            slant_range_time = slant_range_time(sar_transect_lat_end_index:sar_transect_lat_start_index,sar_transect_lon_start_index:sar_transect_lon_end_index); % center of the transect
            % sar_slant_range = sar_slant_range_time(sar_transect_size/2,sar_transect_size/2) .* RTT;

            slant_range_to_transect_center = mean(slant_range_time(:)) .* RTT;
        end

        function [incidence_angle_at_transect] = getIncidenceAngle(obj,sar_transect_lat_end_index,sar_transect_lat_start_index,sar_transect_lon_start_index,sar_transect_lon_end_index)
            incident_angle_grid = ncread(obj.Filepath,"incident_angle");

            incidence_angle_degrees_transect = incident_angle_grid(sar_transect_lat_end_index:sar_transect_lat_start_index,sar_transect_lon_start_index:sar_transect_lon_end_index); % center of the transect
            % sar_center_incidence_angle_degrees =
            % sar_center_incidence_angle_degrees_transect(sar_transect_size/2,sar_transect_size/2);
            % % get center pixel
            % sar_center_incidence_angle_degrees =
            % sar_center_incidence_angle_degrees *
            % ones(sar_sub_transect_size); % greate matrix

            incidence_angle_at_transect = mean(incidence_angle_degrees_transect(:));
        end

        %------------------------------------------------------------------
        function dB_sar_image = convertToDbs(linear_sar_image)
            %CONVERTTODBS convert the linear image to decibels
            dB_sar_image = 10*log10(linear_sar_image);
        end

        function azimuth_to_north_angle = azimuthToNorthAngleConversion(obj)
            % this is the angle between the azimuth plane of the SAR image
            % and the true north direction.
            if obj.LookDirection == "right"
                azimuth_to_north_angle = 180 - obj.SceneOrientationAngle ; % degrees
            end
        end
        
    end
end