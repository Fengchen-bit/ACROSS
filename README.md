# ACROSS
A program collection for ACROSS data analysis and visualization

This project provides a toolkit to help researchers with **ACROSS data analysis and visualization**. Due to the large size of the datasets, we have stored the data separately. This repository primarily contains the programs and scripts we use for data processing and visualization.

Due to the manual quality control required for raw data processing and transfer function generation (including verification of source oscillation frequency, weight phase, and seismometer operational status) , we start with the processed transfer functions as our baseline. This data can be downloaded from [Mendeley Data](https://data.mendeley.com/preview/p7nw36tbbf?a=67361f27-f595-4482-a322-f00a0f994cab), which includes one month of data. Given the substantial volume of data, additional data can be obtained by contacting the *corresponding author*.

##  Code Analysis Programs
_The following are brief descriptions of the data analysis programs in **code_analysis** folder. Some results are saved in the **data** folder, while others require separate downloads due to file size constraints._\
**Extract_transfer_function.m**: Compiles and saves hourly transfer functions for the next processing step.\
**Select_transfer_function_Sum_power_spetral_time_ver_XX.m** :Filters for higher signal-to-noise ratio data and stacks them to obtain daily and reference data for subsequent velocity change calculations in the **XX** component of transfer function.(Output: download from [here](https://data.mendeley.com/preview/p7nw36tbbf?a=67361f27-f595-4482-a322-f00a0f994cab))\
**Cal_dv_v_stretch_XX.m**: Calculates velocity changes (dv/v) for **XX** component. (Output: **dv_v_XX**)\
**Cal_dv_v_PPC_Relation_z_2019_2021_XX.m**: Simulates horizontal or vertical components. (Output: download from [here](https://data.mendeley.com/preview/sjzhv63hbc?a=d9b90b93-8b96-4e30-a23f-232af9632684))
##  Data Files
_The following are brief comments for the **data** files, which are generated from **code_analysis** or prepared in advance for plotting figure_\
**Coda_XX.mat**: Measured end time of coda wave in **XX** component.\
**dv_v_XX.mat**: Measured seismic velocity changes (dv/v) in **XX** component.\
**Mycolormaps.mat**: Colorbar data.\
**nan_day.mat**: Excluded day data.\
**new_seasonal_data.mat**: Rainfall and temperature data\
**rain_2019.mat**: Rainfall data for 2019.\
**rain_2020.mat**: Rainfall data for 2020 (before observation period).\
**rain_data.mat**: rainfall data during the observation period. 
##  Plotting Scripts
_The following are brief descriptions of plotting figure used in the GRL publication_\
**Plot_dv_v_all_TF_cmp_colorbar.m**: Visualizes seismic velocity changes across different components (Ur, Rr, Tr, Ut, Rt, Tt) over time.\
**Plot_dv_v_short_long_term_variation.m**: Compares Seismic Velocity Change (dv/v) with Environmental Factors Correlation Analysis.\
**Plot_dv_v_PPC_2019_2021_4subfigure.m**: Compares measured seismic velocity changes with synthetic velocity changes (Simulated results).\
**Plot_PPC_relationship_2019_2021_contour.m**: Generates RMS(δ²) contour plots for different zeta (ζ) and starting dates. The analysis helps identify optimal parameters for poroelastic modeling of seismic velocity changes in response to precipitation loading.
