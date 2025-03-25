function transects = transectLibrary
% Wave Function Library 
%   This contains all of the functions relating to waves

transects.transectFromSARImage = @transectFromSARImage; % 

end

function [sar_transect,sar_transect_lonGrid,sar_transect_latGrid, sar_transect_center_longitude,lon_indices,sar_transect_center_latitude,lat_indices] = transectFromSARImage(latitude_of_interest,longitude_of_interest,sar_latGrid,sar_lonGrid,sar_Data, sar_large_transect_size)
%TransectFromSARImage take subset of image centered around specified coords
%   sar_large_transect_size = the transect size before filtering (usually
%           512) down to the smaller transect size

% Find the matches for the Longitude grid
    lat_difference = abs(sar_latGrid(:,1) - latitude_of_interest);
    min_difference = min(lat_difference(:));
    lat_start_i_all = find(lat_difference == min_difference); % Alternatively: find(sar_latGrid == min_difference, 1, 'first');
    sar_transect_lat_center_index = lat_start_i_all(1);

    lon_difference = abs(sar_lonGrid(sar_transect_lat_center_index,:) - longitude_of_interest);
    min_difference = min(lon_difference(:));
    lon_start_i_all = find(lon_difference == min_difference);
    sar_transect_lon_center_index = lon_start_i_all(1);
    sar_transect_center_longitude = sar_lonGrid(sar_transect_lat_center_index,sar_transect_lon_center_index);

    lat_indices_1 = sar_transect_lat_center_index-sar_large_transect_size/2: sar_transect_lat_center_index+sar_large_transect_size/2-1;
    lon_indices = sar_transect_lon_center_index-sar_large_transect_size/2: sar_transect_lon_center_index+sar_large_transect_size/2-1;
    if(max(lon_indices)>size(sar_lonGrid,2))
        shift = (max(lon_indices)-size(sar_lonGrid,2));
        lon_indices = lon_indices - (shift);
        sar_transect_center_longitude = sar_lonGrid(sar_transect_lat_center_index, sar_transect_lon_center_index-shift);
    end
    sar_transect_lonGrid = sar_lonGrid(lat_indices_1,lon_indices);


% Now find the matches for lat grid
    lat_difference = abs(sar_latGrid(:,sar_transect_lon_center_index) - latitude_of_interest);
    min_difference = min(lat_difference(:));
    lat_start_i_all = find(lat_difference == min_difference); % Alternatively: find(sar_latGrid == min_difference, 1, 'first');
    sar_transect_lat_center_index = lat_start_i_all(1);
    sar_transect_center_latitude = sar_latGrid(sar_transect_lat_center_index,sar_transect_lon_center_index);

    lat_indices = sar_transect_lat_center_index-sar_large_transect_size/2: sar_transect_lat_center_index+sar_large_transect_size/2-1;
    sar_transect_latGrid = sar_latGrid(lat_indices,lon_indices);
    
    sar_transect = sar_Data(lat_indices,lon_indices);

%% Original version: new version does selection based on smallest difference. This solved a bug where no match was found because of rounding.
    % Find the matches for the Longitude grid
    %     [lat_start_i_all] = find(round(sar_latGrid(:,1),4)== round(latitude_of_interest,4));
    %     sar_transect_lat_center_index = lat_start_i_all(1);
    %     [lon_start_i_all] = find(round(sar_lonGrid(sar_transect_lat_center_index,:),4)== round(longitude_of_interest,4));
    %     sar_transect_lon_center_index = lon_start_i_all(1);
    %     sar_transect_center_longitude = sar_lonGrid(sar_transect_lat_center_index,sar_transect_lon_center_index);
    % 
    %     lat_indices_1 = sar_transect_lat_center_index-sar_large_transect_size/2: sar_transect_lat_center_index+sar_large_transect_size/2-1;
    %     lon_indices = sar_transect_lon_center_index-sar_large_transect_size/2: sar_transect_lon_center_index+sar_large_transect_size/2-1;
    % 
    %     sar_transect_lonGrid = sar_lonGrid(lat_indices_1,lon_indices);
    % 
    % 
    % % Now find the matches for lat grid
    %     [lat_start_i_all] = find(round(sar_latGrid(:,sar_transect_lon_center_index),4)== round(latitude_of_interest,4));
    %     sar_transect_lat_center_index = lat_start_i_all(1);
    %     sar_transect_center_latitude = sar_latGrid(sar_transect_lat_center_index,sar_transect_lon_center_index);
    % 
    %     lat_indices =sar_transect_lat_center_index-sar_large_transect_size/2: sar_transect_lat_center_index+sar_large_transect_size/2-1;
    %     sar_transect_latGrid = sar_latGrid(lat_indices,lon_indices);
    % 
    %     % lat_difference = abs(sar_latGrid - latitude_of_interest);
    %     % min_difference = min(lat_difference(:));
    %     % [min_row, min_column] = find(lat_difference == min_difference); % Alternatively: find(sar_latGrid == min_difference, 1, 'first');
    %     % sar_transect_latGrid = sar_latGrid(min_row(1),min_column(1));
    %     % 
    %     sar_transect = sar_Data(lat_indices,lon_indices);

end