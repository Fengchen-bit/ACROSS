%% 6 Seismic Velocity Change (dv/v) of each transfer function
% This script processes and visualizes seismic velocity changes across 
% different components (Ur, Rr, Tr, Ut, Rt, Tt) over time.
% 
% Input files required:
% - Mycolormaps.mat: Custom colormap for visualization
% - nan_day.mat: Days with NaN values to be handled
% - dv_v_*.mat: Velocity change data for each component
%
% Output: Multi-panel figure showing dv/v variations with error bars
%
% Author: FENG CHEN

clc         
clear      
close all  
load('Mycolormaps.mat');    % Load custom colormap for plotting
load('nan_day.mat');        % Load indices of days with NaN values
load('dv_v_Ur.mat');        % Load initial Ur component data 

% Initialize arrays for statistical analysis
dv_mean = nan(length(t), 6);    % Mean dv/v values for each component and day
dv_se = nan(length(t), 6);      % Standard error for each component and day

%% Define seismic components
% Seismic wave components: Radial (R), Transverse (T), Vertical (Z/U)
% Subscripts: r = Rayleigh wave, t = transverse/Love wave
components = {'Ur', 'Rr', 'Tr', 'Ut', 'Rt', 'Tt'};

%% Process each seismic component
for i = 1:length(components)
    % Load component-specific velocity change data
    load(['dv_v_' components{i} '.mat']);
    % Construct variable name for dynamic evaluation
    var_name = ['dv_v_' components{i}];
    % Calculate daily statistics for each component
    for day_idx = 1:length(t)
        % Find valid (non-NaN) station indices for current day
        st_idx = find(~isnan(eval([var_name '(day_idx,:)'])));
        if ~isempty(st_idx)
            % Calculate mean dv/v across all valid stations
            dv_mean(day_idx, i) = mean(eval([var_name '(day_idx,st_idx)']));
            % Calculate standard error: SE = std/sqrt(n)
            dv_se(day_idx, i) = std(eval([var_name '(day_idx,st_idx)'])) / sqrt(length(st_idx));
        end
    end
    % Store raw data in 3D array for plotting individual stations
    dv_stack(:, :, i) = eval(var_name);
end

%% Create multi-panel figure
fg = figure('Position', [50, 100, 1300, 900]);  % Set figure size and position
ax = gobjects(1, 6);                            % Pre-allocate axes objects

% Define subplot layout parameters
st = [0.1 0.1];         % Starting position [left, bottom]
bd = [0.75 0.75/7];     % Panel dimensions [width, height]

%% Set colormap parameters for error visualization
color_min = 0.01;       % Minimum error value for colormap (1%)
color_max = 0.05;       % Maximum error value for colormap (5%)

% Set error values to zero for NaN days to avoid plotting issues
dv_se(nan_day, :) = 0;

ch_num = 14;            % Number of seismic stations/channels

%% Create subplots for each seismic component
for fig_index = 1:6
    % Create subplot with specified position
    ax(fig_index) = axes(fg, 'Position', [st(1), st(2) + bd(2) * (fig_index - 1), bd]);
    
    if fig_index == 1  % Bottom panel (Ur component) with time axis
        hold on
        
        % Plot individual station data as gray scatter points
        for kk = 1:ch_num
           scatter(ax(fig_index), datenum(t), dv_stack(:, kk, fig_index), 15, ...
               'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
            hold on
        end
        
        % Plot mean values with color-coded standard errors
        scatter(ax(fig_index), datenum(t), dv_mean(:, fig_index), 15, dv_se(:, fig_index), 'filled');
        
        % Apply colormap and set color limits
        colormap(ax(fig_index), mycmap);
        clim(ax(fig_index), [color_min color_max]);
        
        % Set time axis limits
        xlim([datenum(t(1)) datenum(t(end))]);
        box on
        
        % Configure time axis ticks (monthly intervals)
        set(ax(1), 'Xtick', [datenum(t(19)) datenum(t(49)) datenum(t(80)) datenum(t(111)) datenum(t(139))...
                           datenum(t(170)) datenum(t(200)) datenum(t(231)) datenum(t(261)) datenum(t(292))]...
                           , 'Xticklabel', {}, 'FontWeight', 'normal', 'FontSize', 13);
        
        % Configure axis properties        
        gx = gca; 
        gx.XAxisLocation = 'bottom'; 
        gx.XAxis.TickDirection = 'out'; 
        
        % Set axis labels
        xlabelObj = xlabel('Calendar time (month:year)');  
        ylabel('dv/v (%)');
        ylim([-1 1])        % Set y-axis limits to ±1%
        
        % Adjust xlabel position
        pos = get(xlabelObj, 'Position');  
        set(xlabelObj, 'Position', [pos(1), pos(2) - 0.65, pos(3)]);  
        
        %% Create custom time axis with month/year labels
        Tdata = t(1:end);
        years = year(Tdata);
        months = month(Tdata);
        % Find first occurrence of each month in each year
        firstOccurrence_month = nan(max(years) - min(years) + 1, 12);  
        for i = min(years):max(years)
            for j = 1:12
                idx = find(years == i & months == j, 1, 'first');
                if ~isempty(idx)
                    firstOccurrence_month(i - min(years) + 1, j) = idx;
                end
            end
        end
        
        % Reshape and clean month occurrence data
        firstOccurrence_month = reshape(firstOccurrence_month', 1, []);
        firstOccurrence_month(find(isnan(firstOccurrence_month))) = [];
        firstOccurrence_month = firstOccurrence_month + 14;  % Offset for positioning
        
        % Create auxiliary axis for custom time labels
        axPos = get(ax(1), 'Position');
        ax1 = axes('Position', [axPos(1), axPos(2) - 0.03, axPos(3), 0.03], 'Visible', 'off');
        
        yeardata = datenum(Tdata(firstOccurrence_month));
        date = nan(length(Tdata), 1);
        
        % Create month labels
        bar(ax1, datenum(Tdata), date, 1, 'white');
        for i = 1:length(firstOccurrence_month)
            text(datenum(Tdata(firstOccurrence_month(i))), 1, ...
                datestr(Tdata(firstOccurrence_month(i)), 'mm'), 'Parent', ax1, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 10);
        end
        
        % Create year labels and markers
        years = unique(year(Tdata));
        firstOccurrence_year2 = zeros(size(years)); 
        for i = 1:length(years)
            firstOccurrence_year2(i) = find(year(Tdata) == years(i), 1, 'first');
        end
        
        if firstOccurrence_year2(1) == 1
            firstOccurrence_year2(1) = [];
        end
        
        % Calculate year midpoints for label positioning
        years = year(Tdata);
        midpoint = nan(max(years) - min(years) + 1, 1);  
        for i = min(years):max(years)
            idxs = find(years == i);  
            if ~isempty(idxs)
                midpoint(i - min(years) + 1) = idxs(round(length(idxs)/2));  
            end
        end
        
        % Add year separators and labels
        date(firstOccurrence_year2) = 1.5;
        hold on
        bar(ax1, datenum(Tdata), date, 1, 'k');
       years = unique(year(Tdata));
        for i = 1:length(midpoint)
            text(datenum(Tdata(midpoint(i))), 0.1, num2str(years(i)), 'Parent', ax1, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 12);
        end
        % Configure auxiliary axis
        xlim([datenum(Tdata(1)) datenum(Tdata(end))])
        set(gca, 'Visible', 'off')
        box off
        % Set figure background to white
        fig = fg;
        fig.Color = 'white';  
        ax = gca;  
        ax.XColor = 'white';  
        ax.YColor = 'white';  

    else  % Upper panels (other components) without time labels
        % Plot individual station data
        for kk = 1:ch_num
           scatter(ax(fig_index), datenum(t), dv_stack(:, kk, fig_index), 15, ...
               'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
            hold on
        end
        % Plot mean values with error coloring
        scatter(ax(fig_index), datenum(t), dv_mean(:, fig_index), 15, dv_se(:, fig_index), 'filled');
        % Configure axis properties
        xlim([datenum(t(1)) datenum(t(end))]);
        set(ax(fig_index), 'Xtick', [], 'Xticklabel', {}, 'FontWeight', 'normal', 'FontSize', 13);
        ylabel('dv/v (%)');
        ylim([-1 1]);       % Consistent y-axis limits across all panels
        box on
    end
    
    % Apply consistent colormap to all panels
    colormap(ax(fig_index), mycmap);
    clim(ax(fig_index), [color_min color_max]);
end

%% Add colorbar for error visualization
c = colorbar;
c.Position = [.88 .095 .02 .325];      % Position colorbar on right side
c.Ticks = [0.01, 0.02, 0.03, 0.04, 0.05];  % Set tick marks at 1-5%

% Calculate position for colorbar label
xPos = c.Position(1) + c.Position(3) + 0.037;  % X position 
yPos = c.Position(2) + c.Position(4)/2.5;      % Y position 

% Add rotated colorbar label
annotation('textbox', [xPos, yPos - 0.05, 0.12, 0.2], ... 
    'String', 'error (%)', 'FontSize', 12, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Rotation', 90);
