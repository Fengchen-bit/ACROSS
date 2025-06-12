# ACROSS
A program collection for ACROSS data analysis and visualization

This project provides a toolkit to help researchers with **ACROSS data analysis and visualization**. Due to the large size of the datasets, we have stored the data separately. This repository primarily contains the programs and scripts we use for data processing and visualization.

Due to the manual quality control required for raw data processing and transfer function generation (including verification of source oscillation frequency, weight phase, and seismometer operational status) , we start with the processed transfer functions as our baseline. This data can be downloaded from [Mendeley Data](data.mendeley.com/preview/p7nw36tbbf?a=373f3eb2-7963-49ac-9005-5fd84556772b) for 2020.10~2021.5 and [Mendeley Data](data.mendeley.com/preview/vkxs7nh7wd?a=2085f575-0a33-46ea-afb1-308d65ba0bc5) for 2021.6~2021.8.

##  Code Analysis Programs
_The following are brief descriptions of the data analysis programs in **code_analysis** folder. Some results are saved in the **data** folder, while others require separate downloads due to file size constraints._\
**Extract_transfer_function.m**: Compiles and saves hourly transfer functions for the next processing step.\
**Select_transfer_function_Sum_power_spetral_time_ver_XX.m** :Filters for higher signal-to-noise ratio data and stacks them to obtain daily and reference data for subsequent velocity change calculations in the **XX** component of transfer function.(**Output: download from [here](https://data.mendeley.com/preview/sjzhv63hbc?a=cbb09aeb-b900-4629-9506-23311be0317f)**)\
**Cal_dv_v_stretch_XX.m**: Calculates velocity changes (dv/v) for **XX** component. (**Output:dv_v_XX.mat**)\
**Cal_dv_v_PPC_Relation_z_2019_2021_XX.m**: Simulates horizontal or vertical components. (**Output: download from [here](https://data.mendeley.com/preview/6bxptzrkph?a=ff4f98cf-dba1-4de2-9c39-0190cbeeccec)**)
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
