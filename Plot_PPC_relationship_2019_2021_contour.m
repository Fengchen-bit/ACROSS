%% Poroelastic Model Parameter Optimization Visualization
% This script creates time-z heatmaps showing the evolution of RMS(δ²) for different zeta (ζ) and starting
% dates. The analysis helps identify optimal parameters for poroelastic
% modeling of seismic velocity changes in response to precipitation loading.
%
% Main objectives:
% 1. Visualize temporal evolution of model fitting quality
% 2. Compare horizontal vs vertical component sensitivities
% 3. Identify optimal and alternative zeta over time
% 4. Create heatmaps with custom colormaps
%
% Input files required:
% - Relation_PPC_z_Vel.mat: Vertical component simulation results
% - Relation_PPC_z_Hor.mat: Horizontal component simulation results  
% - extract_rain_index.mat: Time indices for analysis period
%
% Author: FENG CHEN

clc
clear
close all

% Load vertical data
load('Relation_PPC_z_Vel.mat','delta_sigma_vel','z_values') % r*c
delta_sigma_vel(isnan(delta_sigma_vel)) = Inf;
% Load horizontal data
load('Relation_PPC_z_Hor.mat','delta_sigma_hor')
delta_sigma_hor(isnan(delta_sigma_hor)) = Inf;
load('extract_rain_index.mat');
t1=datetime([2019 1 1]);
% t2=datetime([2020 11 17]);
t2=datetime([2020 10 14]);
tmp_TT=t1:t2;
TT=tmp_TT(extracted_index);
new_delta_sigma_vel=delta_sigma_vel(:,extracted_index);
new_delta_sigma_hor=delta_sigma_hor(:,extracted_index);

vel_finite = new_delta_sigma_vel(isfinite(new_delta_sigma_vel));
hor_finite = new_delta_sigma_hor(isfinite(new_delta_sigma_hor));

%% Create figure with doubled height for two subfigure
fg = figure('Position', [50, 100, 500, 600]); % Doubled height from 300 to 600
% Create custom colormaps with improved transitions
% Create yellow-red colormap for vertical with SMOOTHER transitions and less darkening
red_steps = 100;
red_map = zeros(red_steps, 3);
for i = 1:red_steps
    % Smoothly transition from bright yellow to red (but not dark red)
    position = (i-1)/(red_steps-1);
    % Red channel stays high throughout (never goes below 0.8)
    red_map(i,1) = 1 - 0.2*position; % Limit darkening
    % Green channel (controls yellow) gradually decreases
    red_map(i,2) = max(0, 1 - 1.5*position);
    % Blue channel stays at 0
    red_map(i,3) = 0;
end

blue_steps = 100;
blue_map = zeros(blue_steps, 3);
for i = 1:blue_steps
    % Smoothly transition from bright blue to green
    position = (i-1)/(blue_steps-1);
    adjusted_position = position * 0.7; 
    % Blue channel gradually decreases (but not below 0.8 at top)
    blue_map(i,3) = 1 - 0.2*adjusted_position; % Limit darkening
    % Green channel gradually increases
    blue_map(i,2) = max(0, 1 - 1.5*adjusted_position);
    % Red channel stays at 0
    blue_map(i,1) = 0;
end
%% Plot vertical change (top) with yellow-red colormap
ax1 = axes(fg, 'Position', [0.12, 0.57, 0.75, 0.38]); % Top plot
TT2=1:length(TT);
h1 = imagesc(ax1, TT2, z_values, new_delta_sigma_vel);
set(ax1, 'YDir', 'normal');
xlim([TT2(1) TT2(end)]);
set(ax1, 'YScale', 'log');
set(ax1, 'XTick', [TT2(1) TT2(191) TT2(415)], 'XTickLabel', {datestr(TT(1), 'yyyy/mm/dd'), datestr(TT(191), 'yyyy/mm/dd'), datestr(TT(415), 'yyyy/mm/dd')});
yticks([1e1 1e2 320 1e3 4200]);
yticklabels({'10^1','10^2', '320','10^3','4200'});
ylim([8 7900]);
ylabel('$\zeta \ (\mathrm{s \cdot m^{\frac{1}{2}}})$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$\mathrm{Starting \ date} \ (\mathrm{yyyy/mm/dd})$', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'ticklength', [0.005 0.025]);
% Set vertical plot to use its own min/max values
vel_min = min(vel_finite);
vel_max = max(vel_finite);
clim([vel_min vel_max]);
colormap(ax1, red_map);
% Add contour lines
hold on;
contour_levels = [4.5 5.5 6.5]; 
[C, h_contour] = contour(ax1, TT2, z_values, new_delta_sigma_vel, contour_levels, 'w-', 'LineWidth', 0.5);
clabel(C, h_contour, 'Color', 'w', 'FontSize', 10, 'LabelSpacing', 400);
hold on;
% Add vertical colorbar with fixed integer ticks (to avoid duplicates)
c1 = colorbar(ax1);
c1.Position = [0.89, 0.57, 0.03, 0.38]; % Position for top colorbar
% Create evenly spaced INTEGER ticks from min to max
min_int = ceil(vel_min);
max_int = floor(vel_max);
cbticks1 = min_int:1:max_int;
set(c1, 'Ticks', cbticks1);
c1.TickLabels = arrayfun(@(x) sprintf('%d', x), cbticks1, 'UniformOutput', false);
c1.FontSize = 12;
% Add delta^2 annotation next to top colorbar
annotation('textbox', [0.89, 0.885, 0.1, 0.1], ... % [x y width height]
    'String', '$\delta^2$', ... % LaTeX formatted text
    'Interpreter', 'latex', ... % Use LaTeX interpreter
    'FontSize', 14, ... % Font size
    'EdgeColor', 'none', ... % No border
    'FitBoxToText', 'on'); % Auto-adjust text box size
% Add vertical short-term and long-term plots
for ii = 1:size(new_delta_sigma_vel, 2)
    [~, min_index_vel(ii)] = min(new_delta_sigma_vel(:, ii));
end
hold on;
p1 = plot(TT2(1:5:end), z_values(min_index_vel(1:5:end)), 'o', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
[~, second_min_pos] = findSecondMinimum(new_delta_sigma_vel);
p2 = plot(TT2(1:5:195), z_values(second_min_pos(1:5:195)), 'o', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
leg1 = legend([p1, p2], {'Global minimum', char('Local minimum', '(secondary)')}, 'box', 'on', 'FontSize', 9);
title(leg1, 'Vertical', 'FontSize', 10,'Color','r');  
set(leg1, 'Position', [0.594666671852271         0.570410397493818         0.275999994814396         0.122499996771415]);
%% Plot horizontal change (bottom) with blue-green colormap
ax2 = axes(fg, 'Position', [0.12, 0.1, 0.75, 0.38]); % Bottom plot
h2 = imagesc(ax2, TT2, z_values, new_delta_sigma_hor);
set(ax2, 'YDir', 'normal');
xlim([TT2(1) TT2(end)]);
set(ax2, 'YScale', 'log');
set(ax2, 'XTick', [TT2(1) TT2(191) TT2(415)], 'XTickLabel', {datestr(TT(1), 'yyyy/mm/dd'), datestr(TT(191), 'yyyy/mm/dd'), datestr(TT(415), 'yyyy/mm/dd')});
yticks([1e1 1e2 320 1e3 4200]);
yticklabels({'10^1','10^2', '320','10^3','4200'});
ylim([8 7900]);
ylabel('$\zeta \ (\mathrm{s \cdot m^{\frac{1}{2}}})$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$\mathrm{Starting \ date} \ (\mathrm{yyyy/mm/dd})$', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'ticklength', [0.005 0.025]);
% Set horizontal plot to use its own min/max values
hor_min = min(hor_finite);
hor_max = max(hor_finite);
clim([hor_min hor_max]);
colormap(ax2, blue_map);
% Add contour lines
hold on;
contour_levels = [1.2 1.8 2.5];
[C, h_contour] = contour(ax2, TT2, z_values, new_delta_sigma_hor, contour_levels, 'w-', 'LineWidth', 0.5);
clabel(C, h_contour, 'Color', 'w', 'FontSize', 10, 'LabelSpacing', 400);
hold on;
% Add horizontal colorbar with fixed integer ticks (to avoid duplicates)
c2 = colorbar(ax2);
c2.Position = [0.89, 0.1, 0.03, 0.38]; % Position for bottom colorbar
% Create evenly spaced INTEGER ticks from min to max
cbticks2 = 1.5:0.5:2.5;
set(c2, 'Ticks', cbticks2);
c2.TickLabels =  {'1.5', '2.0', '2.5'};
c2.FontSize = 12;
% Add delta^2 annotation next to bottom colorbar
annotation('textbox', [0.89, 0.415, 0.1, 0.1], ... % [x y width height]
    'String', '$\delta^2$', ... % LaTeX formatted text
    'Interpreter', 'latex', ... % Use LaTeX interpreter
    'FontSize', 14, ... % Font size
    'EdgeColor', 'none', ... % No border
        'FitBoxToText', 'on'); % Auto-adjust text box size
% Add vertical short-term and long-term plots
for ii = 1:size(new_delta_sigma_hor, 2)
    [~, min_index_hor(ii)] = min(new_delta_sigma_hor(:, ii));
end
hold on;
p1 = plot(TT2(1:5:end), z_values(min_index_vel(1:5:end)), 'o', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
[second_min_vals, second_min_pos] = findSecondMinimum(new_delta_sigma_hor(1:800,:));
p2 = plot(TT2(1:5:195), z_values(second_min_pos(1:5:195)), 'o', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
leg2 = legend(ax2,[p1, p2], {'Global minimum', char('Local minimum', '(secondary)')}, 'box', 'on', 'FontSize', 9);
title(leg2, 'Horizontal', 'FontSize', 10,'Color','g');  
set(leg2, 'Position', [0.594666671852271         0.100410397493818         0.275999994814396         0.122499996771415]);

function [second_min_values, second_min_positions] = findSecondMinimum(data)
% FINDSECONDMINIMUM Finds the second smallest local minimum in each column of data
[rows, cols] = size(data);
second_min_values = nan(1, cols);
second_min_positions = nan(1, cols);
for col = 1:cols
    % Extract current column data for analysis
    column_data = data(:, col);
    % Initialize logical array to mark local minima positions
    % local_mins(i) = 1 if position i is a local minimum, 0 otherwise
    local_mins = zeros(rows, 1);
    for i = 2:(rows-1)
        if ~isnan(column_data(i)) && ~isnan(column_data(i-1)) && ~isnan(column_data(i+1))
            if column_data(i) < column_data(i-1) && column_data(i) <= column_data(i+1)
                local_mins(i) = 1;
            elseif column_data(i) <= column_data(i-1) && column_data(i) < column_data(i+1)
                local_mins(i) = 1;
            end
        end
    end
    %  Handle boundary conditions
    % Boundary points cannot be compared with neighbors on both sides, so they require special treatment
    % Check first point: is it a minimum compared to its right neighbor?
    if rows > 1 && ~isnan(column_data(1)) && ~isnan(column_data(2))
        if column_data(1) < column_data(2)
            local_mins(1) = 1;
        end
    end
    % Check last point: is it a minimum compared to its left neighbor?
    if rows > 1 && ~isnan(column_data(rows)) && ~isnan(column_data(rows-1))
        if column_data(rows) < column_data(rows-1)
            local_mins(rows) = 1;
        end
    end
    % Extract indices of all identified local minima
    min_indices = find(local_mins == 1);
    % Determine second minimum based on number of local minima found
    if length(min_indices) >= 2
        % Multiple local minima found
        min_values = column_data(min_indices);
        [sorted_mins, sorted_idx] = sort(min_values);
        second_min_position = min_indices(sorted_idx(2));
        second_min_value = sorted_mins(2);
        second_min_values(col) = second_min_value;
        second_min_positions(col) = second_min_position;
    elseif length(min_indices) == 1
        % Only one local minimum found
        second_min_values(col) = nan;
        second_min_positions(col) = nan;
    else
        % No local minima found (monotonic or flat data)
        valid_data = column_data(~isnan(column_data));
        if length(valid_data) >= 2
            [sorted_data, sorted_idx] = sort(column_data);
            valid_idx = ~isnan(sorted_data);
            sorted_data = sorted_data(valid_idx);
            sorted_idx = sorted_idx(valid_idx);
            if length(unique(sorted_data)) >= 2
                [unique_vals, unique_idx] = unique(sorted_data, 'first');
                second_min_values(col) = unique_vals(2);
                second_min_positions(col) = sorted_idx(unique_idx(2));
            end
        end
    end
end
end