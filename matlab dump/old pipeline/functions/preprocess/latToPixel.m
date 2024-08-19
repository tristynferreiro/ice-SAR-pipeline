function [topLeftPixelx, topLeftPixely] = latToPixel(lat_grid, lon_grid, top_left_lat, top_left_lon)
    lat_index_col = []; lon_index_col = [];

    for i = 5:-1:0 % number of decimal places / tolerance
        % NOTE: number of decimals also affects/tells us about the satellite resolution https://en.wikipedia.org/wiki/Decimal_degrees

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
    %[index_row, index_col] = find(round(grid,decimal_places) == top_left_coord);
    % Truncate the extra decimal places:
    grid_truncated = floor(grid*10^decimal_places)/10^decimal_places;
    top_lef_coord_truncated = floor(top_left_coord*10^decimal_places)/10^decimal_places;
    [index_row, index_col] = find( grid_truncated == top_lef_coord_truncated);
end