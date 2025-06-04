%% Transfer Function Selection and Processing Using ACROSS Method
% This script processes ACROSS data to extract high-quality transfer functions for velocity change analysis.
% The method applies statistical outlier detection and temporal segmentation to ensure
% data quality and construct reliable reference traces for interferometric analysis.
%
% ACROSS Method Background:
% - Uses controlled seismic sources with precisely known timing and amplitude
% - Provides high signal-to-noise ratio for transfer function estimation
% - Enables detection of small velocity changes through repeat measurements
% - Transfer functions represent impulse response of the medium
%
% Input files required:
% - URT_Tser_*.mat: Time series transfer function data files
% - tf_nama_data.mat: Collected transfer function dataset
%
% Output: Ur_data.mat containing reference and daily transfer functions
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
t2=datetime([2021 8 31]); % end time of observation period
t=t1:t2;
Seisnum={'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14'};% Station numbering (14 seismic stations)
expected_dates = t1+hours(1):hours(2):t2+hours(23); % Create time axis for expected measurements (every 2 hours)
% Define daily boundaries for processing
current_day = dateshift(t1, 'start', 'day');
end_day = dateshift(t2, 'start', 'day');

%% Extract Rr component transfer functions
load('tf_nama_data.mat','tfdata','traces_URTt'); % Collected transfer function dataset
for ii=1:length(tfdata)
    Ur_data(:,ii,:)=tfdata{ii}(:,:,1); % 5000*3864*14 Ur  % - Components Ur-1 Rr2- Tr-3 Ut-4 Ur-5 Tt-6
    % Data structure: (samples × time × stations)
end

tt=traces_URTt.t;
Power_data=nan(length(expected_dates),14);

while current_day <= end_day
    daily_data_present = false;
    for hour = 1:12  % Check 12 measurement times per day (every 2 hours: 1, 3, 5, ..., 23)
        current_time = current_day + hours(2 * (hour - 1) + 1);
        idx = find(fdates == current_time);  % Find corresponding data file
        expected_idx = find(expected_dates == current_time);
        if ~isempty(idx)
            % Calculate FFT for current time
            FFt_data=fft(squeeze(Ur_data(:,expected_idx,:))); % 5000*14 fft
            Power_data(expected_idx,:)=sum(FFt_data(3.01*fsamp:7.15*fsamp,:).*conj(FFt_data(3.01*fsamp:7.15*fsamp,:)),1);
            daily_data_present = true;
        end
    end
    if ~daily_data_present
        % Log days with missing data
        fprintf('No data: %s\n', datestr(current_day, 'yyyy-mm-dd'));
    end
    current_day = current_day + days(1);
end

%% Define operational periods based on power outages and system maintenance
% These periods reflect times when ACROSS system was operational
% Gaps correspond to power outages, maintenance, or system failures

% 2020/10/14 15:00 2020/11/12 11:00
Period{1}=[datetime([2020 10 14 15 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2020 11 12 11 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2020/11/12 17:00 2020/12/17 07:00
Period{2}=[datetime([2020 11 12 17 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2020 12 17 07 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2020/12/17 15:00 2021/01/07 11:00
Period{3}=[datetime([2020 12 17 15 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 01 08 11 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2021/02/03 17:00 2021/03/18 07:00
Period{4}=[datetime([2021 02 03 17 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 03 18 07 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2021/03/18 23:00 2021/04/01 09:00
Period{5}=[datetime([2021 03 18 23 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 04 01 09 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2021/04/01 17:00 2021/04/22 19:00
Period{6}=[datetime([2021 04 01 17 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 04 22 19 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2021/04/23 03:00 2021/05/06 09:00
Period{7}=[datetime([2021 04 23 03 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 05 06 09 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2021/05/06 17:00 2021/06/08 15:00
Period{8}=[datetime([2021 05 06 17 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 06 08 15 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2021/06/08 19:00 2021/06/18 07:00
Period{9}=[datetime([2021 06 08 19 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 06 18 07 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2021/06/18 13:00 2021/07/05 09:00
Period{10}=[datetime([2021 06 18 13 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 07 05 09 00 00],'InputFormat','yyyyMMddHHmmss')]; %
%2021/07/16 13:00 2021/08/31 23:00
Period{11}=[datetime([2021 07 16 13 00 00],'InputFormat','yyyyMMddHHmmss')  datetime([2021 08 31 23 00 00],'InputFormat','yyyyMMddHHmmss')]; %

%% Apply statistical outlier detection for each operational period
% Calculate period-wise statistics to identify problematic stations/times
for ii=1:length(Period)
    % Find data indices within current period
    Period_idx{ii}=find(expected_dates>=Period{ii}(1) & expected_dates<=Period{ii}(2));
    nanidx=find(isnan(Power_data(Period_idx{ii},1)));
    Period_idx{ii}(nanidx)=[];
    for jj=1:length(Seisnum) % 14
        tmp_power=Power_data(Period_idx{ii},jj);
        q1 = quantile(tmp_power,0.25);
        q2=  quantile(tmp_power,0.5);
        q3 = quantile(tmp_power,0.75);
        iqr = q3 - q1;
        upper_whisker = min(max(tmp_power),q3+1.5*iqr);
        lower_whisker = max(min(tmp_power),q1-1.5*iqr);
        % Calculate mean excluding outliers
        Power_data_period_mean(jj,ii)=mean(tmp_power(tmp_power>=lower_whisker & tmp_power<=upper_whisker));
    end
end

%% Global outlier detection across all periods and stations
% Identify periods/stations with anomalous power levels
q1 = quantile(Power_data_period_mean(:),0.25);
q2=  quantile(Power_data_period_mean(:),0.5);
q3 = quantile(Power_data_period_mean(:),0.75);
iqr = q3 - q1;
upper_whisker = min(max(Power_data_period_mean(:)),q3+1.5*iqr);
lower_whisker = max(min(Power_data_period_mean(:)),q1-1.5*iqr);

for ii=1:length(Period) %11
    no_use_st{ii}=find(Power_data_period_mean(:,ii)<=lower_whisker | Power_data_period_mean(:,ii)>=upper_whisker);
end

nan_base_idx=find(isnan(Power_data(:,1)));
for seisidx=1:length(Seisnum) % 14
    seis_nan_idx{seisidx}=nan_base_idx;
end

for ii=1:length(Period) %11
    for seisidx=1:length(no_use_st{ii})
        seis_nan_idx{no_use_st{ii}(seisidx)}=union(seis_nan_idx{no_use_st{ii}(seisidx)},Period_idx{ii});
    end
end

%% Calculate reference wavefrom for each station
% Reference wavefrom serve as baseline for velocity change detection
for seisidx=1:14
    tmp_valid_idx = setdiff(1:size(Ur_data, 2), seis_nan_idx{seisidx});
    tmp_power=Power_data(tmp_valid_idx,seisidx);
    q1 = quantile(tmp_power,0.25);
    q3 = quantile(tmp_power,0.75);
    iqr = q3 - q1;
    upper_whisker = min(max(tmp_power),q3+1.5*iqr);
    lower_whisker = max(min(tmp_power),q1-1.5*iqr);
    % Find remaining outliers
    outlier_positions=find(Power_data(:,seisidx) <= lower_whisker | (Power_data(:,seisidx) >= upper_whisker)& ~isnan(Power_data(:,seisidx)));
    valid_idx{seisidx}=setdiff(tmp_valid_idx, outlier_positions); % Final valid indices after all quality control
    hour_idx{seisidx}=expected_dates(valid_idx{seisidx});
    ref_data(:,seisidx)=mean(squeeze(Ur_data(:,valid_idx{seisidx},seisidx)),2);
    daily_data(:,:,seisidx)=nan(5000,length(t));
    % plot(tt,ref_data(:,seisidx));
    % hold on
end

%% Generate daily averaged transfer functions
% Create daily averages for velocity change analysis
current_day = dateshift(t1, 'start', 'day');
day_idx=1;
while current_day <= end_day
    daily_data_present = false;  % Process each station for current day
    for seisidx=1:14
        exeidx=find(hour_idx{seisidx}>=current_day&hour_idx{seisidx}<=current_day+caldays(1));
        % Require minimum 6 measurements per day for reliable average
        if ~isempty(exeidx) & length(exeidx)>=6
            daily_data(:,day_idx,seisidx)=mean(squeeze(Ur_data(:,valid_idx{seisidx}(exeidx),seisidx)),2);
        end
    end
    day_idx=day_idx+1;
    current_day = current_day + days(1);
end

save('Ur_data','t','tt','ref_data','daily_data');