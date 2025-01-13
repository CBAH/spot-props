118---IC5070 contains files related to HOYS IX publication, though the methodology was established in HOYS VII  --- https://ui.adsabs.harvard.edu/abs/2023MNRAS.520.5433H/abstract, https://ui.adsabs.harvard.edu/abs/2024MNRAS.529.4856H/abstract 

mk_lc_fmt_5 is the file which makes the lightcurves the right format (and is the 5th iteration). this slices years long lightcurves into 6 month slices, half overlapping with the previous slice, to identify the amplitude of the variable light curve over time

118_make_error_files is run to generate text files which vary observed amplitudes within their photometric errors 10000 times and fits the best spot model each time, these are put to text

118_plot_error_files takes those 10000 iterations, plots and saves the output. Medians and median absolute errors are taken and saved as the spot properties and uncertainties. 

phaseplots_13 ensures slices meet the criteria for phase and snr (and creates seperate tables for the data that doesn't meet these criteria), then plot all suitable data points in a 4-panel plot for each object which shows the amplitude, phase, spot temperature and spot size over time. 
