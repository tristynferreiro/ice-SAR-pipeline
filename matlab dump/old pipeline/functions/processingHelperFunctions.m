function func = processingHelperFunctions
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    func.latToPixel = @latToPixel;
    func.compareWithinTolerance = @compareWithinTolerance;
end

function [topLeftPixelx, topLeftPixely] = latToPixel(lat_grid, lon_grid, top_left_lat, top_left_lon)
    lat_index_col = []; lon_index_col = [];

    for i = 10:-1:0 % number of decimal places / tolerance
        
        % Exit the loop if matches have been found
        if(~isempty(lat_index_col) && ~isempty(lon_index_col))
            break
        end
        if(isempty(lat_index_col))
            [lat_index_row,lat_index_col] = compareWithinTolerance(lat_grid, top_left_lat, i);
        end
        if (isempty(lon_index_col)) 
            [lon_index_row,lon_index_col] = compareWithinTolerance(lon_grid, top_left_lon, i);
        end
    end
    
    % Find the indices that are common between latitude and longitude
    commonRowIndices = intersect(lat_index_row, lon_index_row);
    commonColIndices = intersect(lat_index_col, lon_index_col);
    
    topLeftPixelx = max(commonRowIndices)-256;
    topLeftPixely = max(commonColIndices)-256;

end

function [index_row, index_col] = compareWithinTolerance(grid, top_left_coord, decimal_places)
    % Do a comparison of values to find matching indices in the grid matrix
    % within certain decimal place tolerance.

    % Truncate the extra decimal places:
    grid_truncated = floor(grid*10^decimal_places)/10^decimal_places;
    top_lef_coord_truncated = floor(top_left_coord*10^decimal_places)/10^decimal_places;
    [index_row, index_col] = find( grid_truncated== top_lef_coord_truncated);
end