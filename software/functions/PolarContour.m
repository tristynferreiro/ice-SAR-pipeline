function [px, py] = PolarContour(dirSpec, freq, maxf, dir, contourLines)
% PolarContour creates a polar/contour plot for the presentation of directional data.
% This plot was originally designed to plot directional frequency spectra from ocean waves.
%   dirSpec - A matrix of directional frequency data, with frequency data
%             propagating along the x-axis (horizontal) and direction across
%             the y-axis (vertical).
%   freq    - An index of frequency intervals in a vertical list.
%   maxf    - The maximum frequency to be displayed on the plot.

% Define direction bins based on dirSpec dimensions (in radians)
% temp = 360 - (360 / size(dirSpec, 1));
% a = linspace(0, temp, size(dirSpec, 1)); 
a = dir;
dirs = degtorad(a);

% Find the index of the frequencies less than or equal to maxf
freqIndex = find(freq <= maxf);

% Select only the directional data for frequencies less than or equal to maxf
dirSpecAtFreq = dirSpec(:, freqIndex);

% Create an artificial dataset to plot in polar coordinates
[df, ddir] = meshgrid(freq(freqIndex), a); 
ddir = degtorad(ddir);                                                        
[px, py] = pol2cart(ddir, df);                                                  
h = polar(px, py);                                                           
% Remove the artificial dataset but keep the axis
delete(h);                                                                  
hold on;                                                                     
% Plot the contour graph onto the polar axis
contour(px, py, dirSpecAtFreq, contourLines);  
% Plot aesthetics
% Rotate and flip the plot to enable 0 degrees = North
view([90 -90]);
% Set axis ticks for defined values (frequency values in this case)
set(gca, 'YTick', []);                                      
colorbar('vert');                                                           
% ylabel('Direction [degrees] / Frequency [Hz]');                             
xlabel('m^2/Hz/rad');
end