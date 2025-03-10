[dir,freq] = ndgrid(era5_direction_bins_original_degrees,era5_freq_bins);
 G = griddedInterpolant(dir,freq,era5_d2fd);
 dir_new = linspace(min(era5_direction_bins_original_degrees), max(era5_direction_bins_original_degrees), 128);
 freq_new = logspace(min(era5_freq_bins), max(era5_freq_bins), 128);