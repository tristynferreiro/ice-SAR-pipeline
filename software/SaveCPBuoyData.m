%%This script is used to write CSIR data text files into .mat files
%% ==================== INITIALISATION ====================
clear; close all; clc;

%% ==================== Raw file data ====================
% The data must manually be copied from the file. Only the data for a
% single date and time should be copied in each time this code is run.
    % StID : Unique 4-character station id
    % Latitude: Latitude (south positive) (decimal degrees)
    % Longitude: Longitude (east positive) (decimal degrees)
    % Date/Time: Date and time (UTC)
    % Hs: Average height of highest 1/3 of all waves (meters)
    % Tp: Spectral peak wave period (seconds)
    % Directn: Peak direction (degree TN)
    % SpreadF: Peak spreading factor (degrees)
    % Instr Code: Instrument code
    % Frequency: Frequency of spectral values (hertz)
    % Energy: Spectral energy density (meters**2 / hertz)
    % Direction: Direction spectrum (degrees TN)
    % Spreading: Spreading factor spectrum (degrees)

% 1. Copy the first row of the data record into the string below
% StID  Date   Time      Hs     Tp  Directn SpreadF Instr Code
first_row = " ";

% 2. Copy the data matrix into the array below
% Frequency  Energy  Direction  Spreading
data_matrix = [
 
];

%% ==================== Process data ====================
% This slicing is done according to the convention given in the file:
%       StID  Date         Time      Hs     Tp  Directn SpreadF Instr Code
%       Frequency      Energy  Direction  Spreading
% 
% If this changes the code will need to be updated.

% Handle the record information
data_first_row = split(first_row);

year = data_first_row(2);
month = str2num(data_first_row(3));
day = data_first_row(4);
time = data_first_row(5); 
reformat_time = extractBetween(time,1,2)+":"+extractBetween(time,3,4)+":00";
month_array = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
buoyData.date = datetime(day+"-"+month_array(month)+"-"+year+" "+reformat_time); %'13-Dec-2024 17:00:00'

buoyData.Hs = str2num(data_first_row(6));
buoyData.Tp = str2num(data_first_row(7));
buoyData.Dir = str2num(data_first_row(8));
buoyData.SpreadF = str2num(data_first_row(9));

% Slice the data matrix to the correct variables
buoyData.frequencies = data_matrix(:,1); 
buoyData.energy = data_matrix(:,2); 
buoyData.directions = data_matrix(:,3); 
buoyData.spreading = data_matrix(:,4); 


%% ==================== Save to file ====================
filename = "/Users/tris/Documents/MSc/data/CapePointBuoy_"+year+ month+day+time+".mat";

% Save the variables into a .mat file
save(filename,'-struct','buoyData'); 

disp("File saved: "+ filename);
% Clear them out of the workspace
clear buoyData;

%% ==================== Test file ====================
% load(filename)