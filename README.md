# ACROSS
A program collection for ACROSS data analysis and visualization

This project provides a toolkit to help researchers with **ACROSS data analysis and visualization**. Due to the large size of the datasets, we have stored the data separately. This repository primarily contains the programs and scripts we use for data processing and visualization.

Considering the manual time correction required for raw data processing and transfer function generation, we start with the processed transfer functions as our baseline. This data can be downloaded from [Mendeley Data](https://data.mendeley.com/preview/p7nw36tbbf?a=67361f27-f595-4482-a322-f00a0f994cab), which includes one month of data. Given the substantial volume of data, additional data can be obtained by contacting the <ins>corresponding author<ins>.

**Extract_transfer_function**: The hourly transfer function can be compiled and saved for next step.
**Select_transfer_function_Sum_power_spetral_time_ver_**,To filter for higher signal-to-noise ratio data and stack them to obtain daily and reference data for subsequent velocity change calculations, we can use  to process different transfer function components.
Subsequently, we can calculate velocity changes (dv/v) using Cal_dv_v_stretch_**.
For the simulation part, we can simulate horizontal or vertical components using Cal_dv_v_PPC_Relation_z_2019_2021_**.
