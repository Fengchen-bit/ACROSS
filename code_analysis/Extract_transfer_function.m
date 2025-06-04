%% ACROSS Transfer Function Data Extraction and Integration
% This script extracts and integrates transfer function data from multiple individual files
% into a unified dataset suitable for seismic velocity change analysis.
%
% - Data collected every 2 hours during operational periods
%
% Main objectives:
% 1. Scan directory for all ACROSS transfer function data files
% 2. Extract temporal information from filenames
% 3. Load and organize transfer function data systematically
% 4. Handle missing data periods with appropriate NaN filling
% 5. Integrate multi-component data (URT components in r and t directions)
% 6. Create unified dataset for velocity change analysis
%
% Input requirements:
% - URT_Tser_*.mat files containing transfer function time series
% - Files named with timestamp format: URT_Tser_yyyyMMddHHmmss.mat
% - Each file contains traces_URTr and traces_URTt structures
%
% Output: tf_nama_data.mat containing integrated transfer function dataset
%
% Data structure:
% - tfdata{time_index}(samples, stations, components)
% - Components 1-3: URT in radial direction URT r
% - Components 4-6: URT in transverse direction URT T
%
% Author: FENG CHEN

clc
clear
basedir='K:\Seismic_Interfermetory_ACROSS\Ikuta_transfer_function\';
matfiles = dir(fullfile(basedir,'URT_Tser_*.mat')); % Expected filename format: URT_Tser_yyyyMMddHHmmss.mat
fnames=char(matfiles.name);
% Timestamp format: yyyyMMddHHmmss (year-month-day-hour-minute-second)
fdates=datetime(fnames(:,10:23),'InputFormat','yyyyMMddHHmmss'); 
t1=datetime([2020 10 14]); % start time of observation period
% t2=datetime([2020 11 13]); 
t2=datetime([2021 8 31]); % end time of observation period

% Create expected measurement times: 1, 3, 5, 7, ..., 23 hours daily
expected_dates = t1+hours(1):hours(2):t2+hours(23); 

% Define daily processing boundaries
current_day = dateshift(t1, 'start', 'day');
end_day = dateshift(t2, 'start', 'day');

while current_day <= end_day
    daily_data_present = false;  % Flag to track data availability for current day
    for hour = 1:12  % ACROSS system measures every 2 hours: 01, 03, 05, ..., 23 hours
        current_time = current_day + hours(2 * (hour - 1) + 1);  
        idx = find(fdates == current_time); % Find corresponding data file for current time
        expected_idx = find(expected_dates == current_time);
        if ~isempty(idx)  
           % Data file exists for current time - load and process
            load(fullfile(basedir,matfiles(idx).name));
            tfdata{expected_idx,1}(1:5000,1:14,1:3)=traces_URTr.data; % 1-3 URT r
            tfdata{expected_idx,1}(1:5000,1:14,4:6)=traces_URTt.data; % 4-6 URT t
            daily_data_present = true;  
       else
            % No data file exists for current time - fill with NaN
            tfdata{expected_idx,1} = nan(5000,14,6);  
        end
    end
    if ~daily_data_present
        % Entire day missing - important for identifying system outages
        fprintf('No data: %s\n', datestr(current_day, 'yyyy-mm-dd'));
    end
    % Advance to next day
    current_day = current_day + days(1);
end

save('tf_nama_data.mat'); % Output file contains complete transfer function dataset (raw) ready for analysis



