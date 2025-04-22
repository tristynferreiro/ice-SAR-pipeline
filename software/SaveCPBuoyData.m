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
first_row = "CP01  2025  02  11 1830   2.040   10.50  224.42   40.18     32";

% 2. Copy the data matrix into the array below
% Frequency  Energy  Direction  Spreading
data_matrix = [
 0.025      0.012    111.900     74.410
     0.030      0.030    201.900     62.220
     0.035      0.029    164.000     59.310
     0.040      0.008    175.200     72.850
     0.045      0.012    154.100     59.760
     0.050      0.025    189.300     51.250
     0.055      0.024    185.100     53.820
     0.060      0.048    182.300     46.440
     0.065      0.069    237.100     68.150
     0.070      0.054    228.700     64.450
     0.075      0.265    221.600     61.210
     0.080      0.982    217.400     33.350
     0.085      1.268    249.800     46.330
     0.090      0.977    252.600     45.540
     0.095      3.118    223.000     24.060
     0.100      2.165    234.300     33.350
     0.110      1.949    211.800     28.310
     0.120      1.119    209.000     33.570
     0.130      0.742    197.700     35.360
     0.140      0.649    189.300     35.140
     0.150      1.182    176.600     23.390
     0.160      2.345    166.800     16.670
     0.170      2.631    164.000     15.440
     0.180      1.394    164.000     17.680
     0.190      1.206    165.400     16.450
     0.200      1.064    166.800     14.550
     0.210      0.973    168.200     17.120
     0.220      0.682    166.800     18.800
     0.230      0.420    168.200     26.740
     0.240      0.630    164.000     18.580
     0.250      0.469    168.200     18.350
     0.260      0.420    166.800     19.140
     0.270      0.289    166.800     22.830
     0.280      0.289    168.200     27.640
     0.290      0.264    159.800     30.100
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
month =data_first_row(3);
day = data_first_row(4);
time = data_first_row(5); 
reformat_time = extractBetween(time,1,2)+":"+extractBetween(time,3,4)+":00";
month_array = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
buoyData.date = datetime(day+"-"+month_array(str2num(data_first_row(3)))+"-"+year+" "+reformat_time); %'13-Dec-2024 17:00:00'

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