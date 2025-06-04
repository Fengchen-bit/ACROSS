# ACROSS
A program collection for ACROSS data analysis and visualization

This project provides a toolkit to help researchers with **ACROSS data analysis and visualization**. Due to the large size of the datasets, we have stored the data separately. This repository primarily contains the programs and scripts we use for data processing and visualization.

Considering the manual time correction required for raw data processing and transfer function generation, we start with the processed transfer functions as our baseline. This data can be downloaded from [Mendeley Data](https://data.mendeley.com/preview/p7nw36tbbf?a=67361f27-f595-4482-a322-f00a0f994cab), which includes one month of data. Given the substantial volume of data, additional data can be obtained by contacting the <ins>corresponding author<ins>.

**Extract_transfer_function**: The hourly transfer function can be compiled and saved for next step.\
**Select_transfer_function_Sum_power_spetral_time_ver_XX** :To filter for higher signal-to-noise ratio data and stack them to obtain daily and reference data for subsequent velocity change calculations in **XX** component of transfer function.\
**Cal_dv_v_stretch_XX**: To calculate velocity changes (dv/v) for **XX** componnet .\
**Cal_dv_v_PPC_Relation_z_2019_2021_XX**: To simulatite horizontal or vertical components.\

**Coda_XX**: Measured end time of coda wave in **XX** component.\
**dv_v_XX**: Measured seismic velocity changes in **XX** component.\
**Mycolormaps**: colorbar\
**nan_day**: no used day\
**new_seasonal_data**: rainfall and temperature data\
**rain_2019**: rainfall data in 2019\
**rain_2020**: rainfall data in 2020 before observation \
**rain_data**: rainfall data in observation period\

**Plot_dv_v_all_TF_cmp_colorbar**: Visualizes seismic velocity changes across different components (Ur, Rr, Tr, Ut, Rt, Tt) over time\
**Plot_dv_v_PPC_2019_2021_4subfigure**: Compares measured seismic velocity changes with synthetic velocity changes (Simulated results)\
**Plot_dv_v_short_long_term_variation**: Compares Seismic Velocity Change (dv/v) with Environmental Factors Correlation Analysis\
**Plot_PPC_relationship_2019_2021_contour**: RMS(δ²) for different zeta (ζ) and starting dates. The analysis helps identify optimal parameters for poroelastic modeling of seismic velocity changes in response to precipitation loading.\
