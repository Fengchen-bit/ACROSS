%% Seismic Velocity Change (dv/v) vs Environmental Factors Correlation Analysis
% This script analyzes the relationship between seismic velocity changes and
% environmental parameters (precipitation and temperature). The analysis combines
% multiple seismic components into horizontal and vertical groups and compares
% them with meteorological data to identify potential causative mechanisms.
%
% Main objectives:
% 1. Calculate mean dv/v for all seismic components
% 2. Group components into horizontal and vertical wave propagation types
% 3. Apply polynomial detrending to identify long-term trends
% 4. Visualize correlations with precipitation (daily and monthly) and temperature
% 5. Create publication-ready figures showing temporal relationships
%
% Input files required:
% - dv_v_*.mat: Velocity change data for each seismic component
% - nan_day.mat: Days with NaN values to be excluded
% - rain_data.mat: Daily precipitation and temperature measurements
%
% Output: Two-panel figure showing dv/v trends and environmental forcing
%
% Author: FENG CHEN

clc
clear
close all
load('dv_v_Ur.mat');
load('nan_day.mat');
dv_mean=nan(length(t),6);
components = {'Ur', 'Rr', 'Tr', 'Ut', 'Rt', 'Tt'};
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

tf_num=6;
for ii=1:tf_num
    dv_mean_fill(:,ii) = fillmissing(dv_mean(:,ii), 'linear');
end
t_data=1:length(t);
dv_mean_fill_comb(:,1)=mean(dv_mean_fill(:,[2,3,5,6]),2);  % [Rr, Tr, Rt, Tt]
dv_mean_fill_comb(:,2)=mean(dv_mean_fill(:,[1,4]),2); % [Ur, Ut]

% Third-order polynomial fitting
for tf_index=1:2
poly_model = polyfit(t_data', dv_mean_fill_comb(:,tf_index), 3);
fitted_dv_v(:,tf_index) = polyval(poly_model, t_data');
end
fitted_dv_v_hor=fitted_dv_v(:,1);
fitted_dv_v_vel=fitted_dv_v(:,2);

fg=figure('Position',[50,100,1300,700]);
ax=gobjects(1,4);
ax(1)=axes(fg,'Position',[0.17,0.2,0.75,0.25]);
hold(ax(1),'on');

scatter(ax(1),datenum(t),(dv_mean(:,1)+dv_mean(:,4))/2,15,'r','filled'); % veritical
scatter(ax(1),datenum(t),(dv_mean(:,2)+dv_mean(:,3)+dv_mean(:,5)+dv_mean(:,6))/4,15,'g','filled'); % horizental
g1=plot(ax(1),datenum(t),fitted_dv_v_hor,'Color',[0.5 0.9 0.5],'LineWidth',3); %b
g1.Color(4) = 0.8; 
g2=plot(ax(1),datenum(t),fitted_dv_v_vel,'Color',[1 0.7 0.7],'LineWidth',3); %b
g2.Color(4) = 0.8; 

xlim([datenum(t(1)) datenum(t(end))]);
box on
set(ax(1),'Xtick',[datenum(t(19)) datenum(t(49)) datenum(t(80)) datenum(t(111)) datenum(t(139))...
    datenum(t(170)) datenum(t(200)) datenum(t(231)) datenum(t(258)) datenum(t(261)) datenum(t(264)) datenum(t(292))]...
    ,'Xticklabel',{},'FontWeight','normal','FontSize',13);
gx = gca; 
gx.XAxisLocation = 'bottom'; 
gx.XAxis.TickDirection = 'out'; 
xlabelObj = xlabel('Calendar time (month:year)');  
ylabel('dv/v (%)');
ylim([-0.75 0.7])
pos = get(xlabelObj, 'Position');  
set(xlabelObj, 'Position', [pos(1), pos(2)-0.13, pos(3)]); 
% Create custom time axis with month/year labels
Tdata=t(1:end);
years = year(Tdata);
months = month(Tdata);
% Calculate first occurrence of each month in each year for labeling
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
firstOccurrence_month=reshape(firstOccurrence_month',1,[]);
firstOccurrence_month(find(isnan(firstOccurrence_month)))=[];
firstOccurrence_month=firstOccurrence_month+14;
% Create auxiliary axis for custom time labels
axPos = get(ax(1), 'Position');
ax1 = axes('Position', [axPos(1), axPos(2)-0.03, axPos(3), 0.03], 'Visible', 'off');
yeardata=datenum(Tdata(firstOccurrence_month));
date=nan(length(Tdata),1);
bar(ax1,datenum(Tdata),date,1,'white');
for i = 1:length(firstOccurrence_month)
    text(datenum(Tdata(firstOccurrence_month(i))), 1, datestr(Tdata(firstOccurrence_month(i)),'mm'), 'Parent', ax1, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center','FontSize',10);
end  
years = unique(year(Tdata));
firstOccurrence_year2 = zeros(size(years));  
for i = 1:length(years)
    firstOccurrence_year2(i) = find(year(Tdata) == years(i), 1, 'first');
end
if firstOccurrence_year2(1)==1
firstOccurrence_year2(1)=[];
end
years = year(Tdata);
midpoint = nan(max(years) - min(years) + 1, 1); 
for i = min(years):max(years)
    idxs = find(years == i);  
    if ~isempty(idxs)
        midpoint(i - min(years) + 1) = idxs(round(length(idxs)/2)); 
    end
end
date(firstOccurrence_year2)=1.5;
hold on
bar(ax1,datenum(Tdata),date,1,'k');
years = unique(year(Tdata));
for i = 1:length(midpoint)
    text(datenum(Tdata(midpoint(i))), 0.1, num2str(years(i)), 'Parent', ax1, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center','FontSize',12);
end  
xlim([datenum(Tdata(1)) datenum(Tdata(end))])
set(gca,'Visible','off')
box off
fig = fg;
fig.Color = 'white';  
gx = gca;  
gx.XColor = 'white'; 
gx.YColor = 'white';  

line(ax(1), [ax(1).XLim(1) ax(1).XLim(2)], [ax(1).YLim(1) ax(1).YLim(1)], 'Color', 'k');
line(ax(1), [ax(1).XLim(1) ax(1).XLim(1)], [ax(1).YLim(1) ax(1).YLim(2)], 'Color', 'k'); 
line(ax(1), [ax(1).XLim(2) ax(1).XLim(2)], [ax(1).YLim(1) ax(1).YLim(2)], 'Color', 'k'); 
% Alternative simpler method:
% Set the axes box to 'off' and then set the top color to 'none'
ax(1).Box = 'off';
ax(1).XAxis.Color = 'k';  % Keep the bottom x-axis visible
ax(1).YAxis.Color = 'k';  % Keep the y-axis visible
% Remove the top edge
ax(1).TickDir = 'out';  % Make ticks point outward (optional)
%% Panel 2: Environmental factor (precipitation and temperature)
% Load meteorological data
load('rain_data.mat','Rain_data','Temperature_data');
month_Rain = [30 43.5 12 95.5 74 261.5 227.5 247.5 149.5 406 327];
Tdata_Month = Tdata([1,21,51,82,113,141,172,202,233,263,294]+15);

ax(2)=axes(fg,'position',[.17 .45 .75 .20]);
Rain_data(1:2)=[];
Rain_data(length(Tdata)+1:end) = [];

yyaxis left
p1=bar(ax(2), datenum(t), -Rain_data, 1, 'b');
scale_factor = 150/400;  % Scale monthly values to fit in daily range
hold on
p2=plot(ax(2), datenum(Tdata_Month), -month_Rain* scale_factor, '-o', 'Color', [0.3010 0.7450 0.9330], 'Linewidth', 2, 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
ylim([-180 0]);
% Position x-axis at the top
set(ax(2), 'XAxisLocation', 'top', 'FontWeight', 'bold', 'FontSize', 13);
% Invert y-axis labels to show positive values
set(ax(2), 'ytick', [-150 -100 -50 0], 'yticklabel', {'150', '100', '50', '0'}, 'FontWeight', 'normal', 'FontSize', 10, 'YColor', 'b')
yl = ylabel('Precipitation (mm)','Color','b');
set(yl, 'Position', get(yl,'Position') + [0 20 0]);

%% Plot monthly precipitation trend (light blue line with markers)
hold on
yyaxis right
Temperature_data(1:2)=[];
Temperature_data(length(Tdata)+1:end)=[];
p3=scatter(ax(2),datenum(t),Temperature_data,15,'m','filled');
ylim([-1.5 30]);
set(ax(2),'ytick',[0 10 20 30],'yticklabel',{'0' ,'10' ,'20','30'},'FontWeight','bold','FontSize',10, 'YColor', 'm')
set(ax(2),'Xtick',[],'Xticklabel',{},'FontWeight','normal','FontSize',10);
xlim([datenum(t(1)) datenum(t(end))]);
bx = gca;
bx.TickLength = [0.005 0.01];
box on
ylabel('Temperature (C^{o})','Color', 'm');
ax_pos = get(ax(2), 'Position');
% Create the second axis - adjusted position to be left of the main axis
ax2 = axes('Position', [ax_pos(1)-0.044, ax_pos(2), ax_pos(3), ax_pos(4)], 'yticklabel', {'400', '300', '200', '100', '0'}, 'Color', 'none');
% Set the y-limits scaled to match the visual representation
set(ax2, 'YLim', [-450 0], 'YTick', [-400 -300 -200 -100 0], 'YAxisLocation', 'left', ...
    'XTick', [], 'XColor', 'none', 'Color', 'none', 'FontWeight', 'normal', 'FontSize', 10, ...
    'YColor', [0.3010 0.7450 0.9330], 'Box', 'off', 'TickLength', [0.01 0.025]);
% Make sure tick marks are visible
ax2.TickDir = 'out';
% Move the ylabel further to the left
ylabel(ax2, 'Monthly precipitation (mm)', 'Color', [0.3010 0.7450 0.9330]);

set(ax(2),'Xtick',[],'Xticklabel',{});
ax(2).Box = 'off';  % Turn off the box around the plot
% Create custom box with only three sides (left, right, and top)
line(ax(2), [ax(2).XLim(1) ax(2).XLim(2)], [ax(2).YLim(2) ax(2).YLim(2)], 'Color', 'k'); % Top
line(ax(2), [ax(2).XLim(1) ax(2).XLim(1)], [ax(2).YLim(1) ax(2).YLim(2)], 'Color', 'k'); % Left
line(ax(2), [ax(2).XLim(2) ax(2).XLim(2)], [ax(2).YLim(1) ax(2).YLim(2)], 'Color', 'k'); % Right
ax(2).Box = 'off';
