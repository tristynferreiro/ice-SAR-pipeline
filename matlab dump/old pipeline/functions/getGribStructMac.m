function struc = getGribStructMac(GribPath)
 % Read in the .grib2 wave buoy data file on MacOS

OutputColumns = ["latitude", "longitude", "time", "significantWaveHeight", "significantWavePeriod", "direction", "windSpeed", "windDirection"];
full_path = dir(strcat(GribPath,'/*.grib2'));
S = struct();
for i = 1:length(OutputColumns)
    S = setfield(S,OutputColumns(i),[]);
end

for i = 1:length(full_path)
    % Define the command to run
    file_path = strcat(full_path(i).folder, '/', full_path(i).name);
    output_file = strcat(full_path(i).folder, '/',full_path(i).name(1:end-6), '.nc') ;
    if ~exist(output_file,"file")
        command = sprintf('wgrib2 %s -netcdf %s >nul 2>&1',file_path, output_file);
        setenv('PATH', getenv('PATH')+":/opt/homebrew/bin/"); % MATLAB needs to know the wgrib path
        system(command);  %convert grib2 to netcdf
    end

    % Get the nc file info 
    fileinfo = ncinfo(output_file);
    % get lon/lat data
    lon = ncread(output_file,'longitude'); % X
    lat = ncread(output_file,'latitude'); % Y
    lat = sort(lat,'descend')';

    % get time data
    time = fileinfo.Variables(3).Attributes(5).Value;
    datevec = datetime(time(1:end-3),"InputFormat","yyyy.MM.dd HH:mm:ss","TimeZone","UTC");
    datevec.Format = "dd-MMM-uuuu HH:mm:ss";
    S(i).latitude = lat;
    S(i).longitude = lon;
    S(i).time = datevec;
    S(i).significantWaveHeight = ncread(output_file,'HTSGW_surface');
    S(i).significantWavePeriod = ncread(output_file,'PERPW_surface');
    S(i).direction = ncread(output_file,'DIRPW_surface');
    S(i).windSpeed = ncread(output_file,'WIND_surface');
    S(i).windDirection = ncread(output_file,'WDIR_surface');
    
    % To save space delete the newly created .nc file
    delete(output_file)
end

struc = S;
end