%% Pore Pressure Change Calculation from Precipitation and dv/v Modeling
% This script calculates pore pressure changes in the subsurface caused by
% precipitation loading and establishes the relationship with measured seismic
% velocity changes (dv/v) from vertical component. 
%
% Physical basis:
% 1. Precipitation adds surface loading, creating pore pressure changes
% 2. Pressure diffuses into subsurface following 1D diffusion equation
% 3. Pore pressure changes affect effective stress and seismic velocities
% 4. Linear relationship
%
% Main objectives:
% 1. Calculate pore pressure changes from daily precipitation data
% 2. Model pressure diffusion using complementary error function solution
% 3. Establish empirical relationship between pressure and velocity changes
% 4. Optimize penetration depth parameter through systematic parameter sweep
% 5. Generate synthetic dv/v predictions for different starting dates
%
% Input files required:
% - dv_v_*.mat: Velocity change data for each seismic component
% - rain_data.mat: Daily precipitation measurements
% - rain_2019.mat, rain_2020.mat: Annual precipitation data
% - nan_day.mat: Days with missing data
%
% Output: Relation_PPC_z_Vel.mat containing optimization results
%
% Author: FENG CHEN

clc
clear
close all
load('dv_v_Ur.mat');
load('nan_day.mat');
load('rain_data.mat','Rain_data'); % Observation period precipitation dataset
load('rain_2019.mat','year_array1'); % 2019 annual data
rain_data_2019=year_array1;
load('rain_2020.mat','year_array'); % 2020 annual data 
rain_data_2020=year_array;

components = {'Ur', 'Rr', 'Tr', 'Ut', 'Rt', 'Tt'};
dv_mean=nan(322,6);
for i = 1:length(components)
    % Load component data
    load(['dv_v_' components{i} '.mat']);
        % Get the variable name dynamically
    var_name = ['dv_v_' components{i}];
        % Calculate mean for non-NaN values for each day
    for day_idx = 1:length(t)
        st_idx = find(~isnan(eval([var_name '(day_idx,:)'])));
        if ~isempty(st_idx)
            dv_mean(day_idx, i) = mean(eval([var_name '(day_idx,st_idx)']));
        end
    end
end

new_fitted_dv_v_vel=mean(dv_mean(:,[1,4]),2); % Average of vertical components
nan_data=isnan(new_fitted_dv_v_vel);
new_dv_v_clean_vel=new_fitted_dv_v_vel(~nan_data);
new_dv_v_mean_vel=mean(new_dv_v_clean_vel);

%% calculated fitted z for vel
Rain_data(1:2)=[];
Rain_data(323:end)=[];
% Combine multi-year precipitation data
new_Rain_data=[rain_data_2019 rain_data_2020 Rain_data];
% This represents precipitation change from previous day
tmp_new_Rain_data=[0 new_Rain_data(1:end-1)];

Tdata=t(1:end);
rho = 1000; % water density (kg/m^3)
g = 9.81; % g (m/s^2)
% Calculate pressure changes from precipitation loading
% ΔP = ρgh, where h is precipitation height change
% Factor 0.001 converts mm to m
delta_p=rho*g*(new_Rain_data-tmp_new_Rain_data)*0.001;

delta_tt = 86400; % Time increment: 1 day in seconds
n = length(delta_p);  % Total number of days in analysis

z_min = 8;  % Minimum value for z
z_max = 7900;  % Maximum value for z 
num_points = 1000;  % Keep the same number of points

% Create logarithmically spaced z values
z_values = logspace(log10(z_min), log10(z_max), num_points);

for day_idx=1:n-length(Tdata)+1
    for jj=1:1000
        z=z_values(jj);
        P = zeros(n, 1); 
        for tt = day_idx+1:n
            tmp_P=0;
            for  ii=day_idx:tt
                tmp_P=tmp_P+delta_p(ii) * erfc(z*1/sqrt(delta_tt*(tt-ii)));
            end
            P(tt)=tmp_P;
        end
        P(1:end-322)=[];
        P_clean=P(~nan_data);
        P_mean=mean(P);
        cov_matrix = cov(new_dv_v_clean_vel, P_clean);
        cov_dv_v_P = cov_matrix(1, 2);
        var_P = var(P_clean);
        % Formula: dv/v_syn = mean(dv/v) + (cov(dv/v,P)/var(P)) × (P - mean(P))
        tmp_dv_v_syn(:, jj, day_idx) =new_dv_v_mean_vel+ cov_dv_v_P / var_P * (P_clean - P_mean);
        delta_sigma_vel(jj, day_idx) = sum((tmp_dv_v_syn(:, jj, day_idx) - new_dv_v_clean_vel).^2);
    end
        % Display progress for current starting date
    disp(['delta_sigma(:, ', num2str(day_idx), ') = ', num2str(delta_sigma_vel(:, day_idx)')]); 
end
save Relation_PPC_z_Vel.mat
