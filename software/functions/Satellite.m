classdef Satellite
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Filepath
        ProductType                    % 

        AcquisitionStartDatetime
        AcquisitionStopDatetime
        Polarisation
        LookDirection
        RangeResolution
        AzimuthResolution
        SlantRangeToFirstPixel
        SlantRange
        SatelliteVelocity               %m/s from Giacomo's output file.
        SceneOrientationAngle
        OrbitType                       %Ascending/Descending 
    end

    methods
        function obj = Satellite(filepath)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Filepath = filepath;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end