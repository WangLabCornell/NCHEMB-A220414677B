Constant-force stretching
1. showtracestretch_hold.m
Script to load the converted data including force and extension from optical tweezers instrument for the constant force experiments and analyze trace-by-trace to obtain the hold time at the set force of 60 pN

2. data_analysis_2.m
Script to fit the remained fraction of tether under constant force over time with a double-explonential function to obtain the tether's half-life (t_1/2) 

3. x_f_WLC.mat
Look-up table of normalized extension and stretching force of dsDNA following the modified Marko-Siggia model for DNA stretching.

Constant-velocity stretching
1. showtracestretch_old.m
Script to load the converted data including force and extension from optical tweezers instrument for the constant velocity experiments and analyze trace-by-trace to obtain tether entension at 0.5pN and tether brakage force

2. print_stretching_traces.m
Script to analyze the force-extension curve trace-by-trace to detect disruption peaks of DNA loop rupture and print the single-molecule traces with information of 0.5-pN extension, loop rupture force, tether breakage force, loop numbers, etc.  

3. data_statistic_stretching_all.m
Script to summarize the data and give statistics for DNA loop distuption force, tether breakage force, 0.5-pN extension

Twisting
1. analyze_topo_dna_with_pause_analysis.m
Script to load the magnetic tweezers data and analyze trace-by-trace to obtain the winding curve of DNA tether and the overall relaxation activity of topo II  

2. pause_analysis_new.m
Script to characterize pausing behavior of topo II during supercoil relaxation based on dwell-time analysis and print out single-molecule traces of magnetic tweezers experiment

3. plot_pause_with_level.m
Script to plot the pausing duration of topo II on time-elapsed traces of tether extension 

4. data_plot_summary_210519.m
Script to summarize the data and give statistics for pausing frequency, pause-free speed, pause duration of topo II.

Unzipping
1. alignment_stalldetection5.m
Script to load the converted data from optical tweezers instrument and perform sequence alignment with unzipping theory 

2. print_saved_trace.m
Script to print out unzipping trace and calculate the sliding distance of bound topo II

3. data_plot_summary_191122.m
Script to summarize the data and give statistics for tether breakage probability, sliding distance, max unzipping force, number of sliding regions.   

Source for other MATLAB codes needed to run our scripts:
movingmean.m: https://www.mathworks.com/matlabcentral/fileexchange/41859-moving-average-function?s_tid=srchtitle
movingslope.m: https://www.mathworks.com/matlabcentral/fileexchange/16997-movingslope
peakfinder.m: https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate
export_fig.m: https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig
