%% Simulation for seismic velocity change
% This script compares measured seismic velocity changes with synthetic velocity changes.

% Main objectives:
% 1. Calculate mean dv/v for horizontal and vertical components
% 2. Perform polynomial trend fitting to separate long-term and short-term variations
% 3. Compare observed dv/v with synthetic velocity changes from poroelastic model
% 4. Visualize correlations between precipitation, pore pressure, and dv/v
%
% Input files required:
% - dv_v_*.mat: Velocity change data for each component of transfer
% functin
% - nan_day.mat: Days with NaN values
% - Relation_PPC_z_Hor.mat: Simulation result for horizontal component
% - Relation_PPC_z_Vel.mat: Simulation result for vertical component
%
% Output: Multi-panel figure showing model validation and environmental correlations
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
%% Fill missing values using linear interpolation
% This step ensures continuous time series for polynomial fitting
for ii=1:length(components)
    dv_mean_fill(:,ii) = fillmissing(dv_mean(:,ii), 'linear');
end

dv_mean_fill_comb(:,1)=mean(dv_mean_fill(:,[2,3,5,6]),2);% Horizontal components: Rr, Tr, Rt, Tt
dv_mean_fill_comb(:,2)=mean(dv_mean_fill(:,[1,4]),2);% Vertical components: Ur, Ut
t_data=1:length(t);

for tf_index=1:2
    poly_model = polyfit(t_data', dv_mean_fill_comb(:,tf_index), 3);
    fitted_dv_v(:,tf_index) = polyval(poly_model, t_data');
end

fitted_dv_v_hor=fitted_dv_v(:,1);
fitted_dv_v_vel=fitted_dv_v(:,2);

fg=figure('Position',[50,100,1800,1300]);
ax=gobjects(1,7);
st=[.025 .08];
bd2=[.145 .6];  
bd1=[.35 .15]; 
% fitting figure
ax(1)=axes(fg,'Position',[st,bd2]);
% dv/v in hor(short-term)
ax(2)=axes(fg,'Position',[st(1)+bd2(1)+0.032,st(2)+0.05,bd1]);
ax(3)=axes(fg,'Position',[st(1)+bd2(1)+0.032,st(2)+bd1(2)+0.05,bd1(1),bd1(2)*0.5]);
% dv/v in vet(Long-term)
ax(4)=axes(fg,'Position',[st(1)+bd2(1)+0.032,st(2)+2*bd1(2)+0.052,bd1]);
ax(5)=axes(fg,'Position',[st(1)+bd2(1)+0.032,st(2)+2*bd1(2)+0.052+bd1(2),bd1(1),bd1(2)*0.5]);

% dv/v in vet(short-term)
ax(6)=axes(fg,'Position',[st(1)+bd2(1)+0.09+bd1(1),st(2)+0.05,bd1]);
ax(7)=axes(fg,'Position',[st(1)+bd2(1)+0.09+bd1(1),st(2)+bd1(2)+0.05,bd1(1),bd1(2)*0.5]);
% dv/v in hor(Long-term)
ax(8)=axes(fg,'Position',[st(1)+bd2(1)+0.09+bd1(1),st(2)+2*bd1(2)+0.052,bd1]);
ax(9)=axes(fg,'Position',[st(1)+bd2(1)+0.09+bd1(1),st(2)+2*bd1(2)+0.052+bd1(2),bd1(1),bd1(2)*0.5]);


%% plot select z
% plot horizental short-term veriation
load('Relation_PPC_z_Hor.mat','delta_sigma_hor','z_values','tmp_dv_v_syn','new_fitted_dv_v_hor');
dv_v_syn_hor=nan(length(t),size(tmp_dv_v_syn,2),size(tmp_dv_v_syn,3));
dv_v_syn_hor(~isnan(dv_mean(:,1)),:,:)=tmp_dv_v_syn;
g1=semilogx(ax(1),z_values,delta_sigma_hor(:,1),'Color','g','LineWidth',2);
g1.Color(4) = 0.5;
hold(ax(1),'on');
% plot vertical short-term veriation
load('Relation_PPC_z_Vel.mat','delta_sigma_vel','tmp_dv_v_syn','new_fitted_dv_v_vel');
dv_v_syn_vel=nan(length(t),size(tmp_dv_v_syn,2),size(tmp_dv_v_syn,3));
dv_v_syn_vel(~isnan(dv_mean(:,1)),:,:)=tmp_dv_v_syn;
g2=semilogx(ax(1),z_values,delta_sigma_vel(:,1),'Color','r','LineWidth',2);
g2.Color(4) = 0.5;
ylabel(ax(1),'\sigma^2','FontSize',15);
xlabel(ax(1),'$\zeta \ (\mathrm{s \cdot m^{\frac{1}{2}}})$', 'Interpreter', 'latex', 'FontSize', 15);
ylim(ax(1),[0 8]);
xlim(ax(1),[8 7900]);
xticks(ax(1),[1e1 1e2 310 1e3 4200]);
xticklabels(ax(1),{'10^1','10^2', '310','10^3','4200'});
legend(ax(1), [g2,g1], ...
    { 'Vertical','Horizontal'}, ...
    'Box', 'off', ...
    'FontSize', 20,'Position', [0.07, 0.3, 0.1, 0.1]);

%% Panel 2: Short-term Vertical dv/v Comparison
% Find optimal z for vertical components (minimum σ²)
[~, min_index_vel] = min(delta_sigma_vel(:, 1));
scatter(ax(2),datenum(t),new_fitted_dv_v_vel,11,'MarkerEdgeColor','r','MarkerFaceColor','r');
hold(ax(2),'on');
plot(ax(2),datenum(t),fitted_dv_v_vel,'Color',[1 0.7 0.7],'LineWidth',3);
scatter(ax(2),datenum(t),dv_v_syn_vel(:,min_index_vel,1),11,'MarkerEdgeColor','k','MarkerFaceColor','k');
box(ax(2),'on');
z_integer = round(z_values(min_index_vel));
text_string = sprintf('$\\zeta = %d \\, \\mathrm{(s \\cdot m^{\\frac{1}{2}})}$', z_integer);
text(ax(2),datenum(t(140)),0.5, text_string, 'Interpreter', 'latex', 'FontSize', 18);
xlim(ax(2),[datenum(t(1)) datenum(t(end))]);
set(ax(2),'Xtick',[datenum(t(19)) datenum(t(49)) datenum(t(80)) datenum(t(111)) datenum(t(139))...
    datenum(t(170)) datenum(t(200)) datenum(t(231)) datenum(t(258)) datenum(t(261)) datenum(t(264)) datenum(t(292))]...
    ,'Xticklabel',{},'FontWeight','normal','FontSize',13);
ax(2).TickLength = [0.005 0.01];
ax(2).FontSize = 11;
axes(ax(2));
gx = gca;
gx.XAxisLocation = 'bottom';
gx.XAxis.TickDirection = 'out';
y1=ylabel('dv/v (%)','FontSize',11);
ylim(ax(2),[-0.75 0.7]);
set(ax(2),'Ytick',[-0.5 0 0.5]);
% ylim([-1 1]);
xlabelObj = xlabel('Calendar time (month:year)');
set(xlabelObj, 'FontSize', 14);
pos = get(xlabelObj, 'Position');
set(xlabelObj, 'Position', [pos(1), pos(2)-0.3, pos(3)]);

% Create custom time axis with month/year labels for Panel 2
Tdata=t(1:end);
years = year(Tdata);
months = month(Tdata);
% Calculate first occurrence of each month in each year
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
axPos = get(ax(2), 'Position');
ax1 = axes('Position', [axPos(1), axPos(2)-0.03, axPos(3), 0.03], 'Visible', 'off');
yeardata=datenum(Tdata(firstOccurrence_month));
date=nan(length(Tdata),1);
bar(ax1,datenum(Tdata),date,1,'white');
% Add month labels
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
cx = gca;
cx.XColor = 'white';
cx.YColor = 'white';

line(ax(2), [ax(2).XLim(1) ax(2).XLim(2)], [ax(2).YLim(1) ax(2).YLim(1)], 'Color', 'k');
line(ax(2), [ax(2).XLim(1) ax(2).XLim(1)], [ax(2).YLim(1) ax(2).YLim(2)], 'Color', 'k');
line(ax(2), [ax(2).XLim(2) ax(2).XLim(2)], [ax(2).YLim(1) ax(2).YLim(2)], 'Color', 'k');
ax(2).Box = 'off';
ax(2).XAxis.Color = 'k';
ax(2).YAxis.Color = 'k';
ax(2).TickDir = 'out';

%% Panel 3: Precipitation and Pore Pressure for Short-term Vertical
% Load precipitation and model parameters
load('Relation_PPC_z_Vel.mat','delta_p','Rain_data');
z = z_values(min_index_vel);
% Calculate pore pressure changes using 1D diffusion model
n = length(delta_p); % day number
P = zeros(n, 1);
delta_tt = 86400; % Time increment in seconds (1 day)
for tt = 2:n
    tmp_P=0;
    for  ii=1:tt
        tmp_P=tmp_P+delta_p(ii) * erfc(z/sqrt(delta_tt*(tt-ii)));
    end
    P(tt)=tmp_P;
end
axes(ax(3));
yyaxis left
bar(ax(3), datenum(t), -Rain_data, 1, 'b');
ylim([-180 0]);
set(ax(3), 'XAxisLocation', 'top', 'FontWeight', 'bold', 'FontSize', 13);
% Invert y-axis labels to show positive values
set(ax(3), 'ytick', [-150 -100 -50 0], 'yticklabel', {'150', '100', '50', '0'}, 'FontWeight', 'normal', 'FontSize', 10)
yl = ylabel('Precipitation (mm)','FontSize',11);
set(yl, 'Position', get(yl,'Position') + [0 20 0]);
hold on
bx = gca;
bx.TickLength = [0.005 0.01];
ax(3).Box = 'off';  % Turn off the box around the plot
% Create custom box with only three sides (left, right, and top)
line(ax(3), [ax(3).XLim(1) ax(3).XLim(2)], [ax(3).YLim(2) ax(3).YLim(2)], 'Color', 'k');
line(ax(3), [ax(3).XLim(1) ax(3).XLim(1)], [ax(3).YLim(1) ax(3).YLim(2)], 'Color', 'k');
line(ax(3), [ax(3).XLim(2) ax(3).XLim(2)], [ax(3).YLim(1) ax(3).YLim(2)], 'Color', 'k');
ax(3).Box = 'off';
% Since your plot has the x-axis at the top, ensure it remains visible
ax(3).XAxis(1).Color = 'k';  % Keep the top x-axis visible
ax(3).YAxis(2).Color = 'k';  % Keep the y-axis visible
hold on
% Right y-axis: Plot pore pressure changes (orange)
yyaxis right
plot(ax(3), datenum(t), -P(end-321:end),'Color',[1, 0.5, 0], 'LineWidth', 1.5);
ax(3).YAxis(2).Color = [1, 0.5, 0];
% Set y-axis ticks with negative values but positive labels
set(ax(3),'Xtick',[],'Xticklabel',{});
xlim(ax(3),[datenum(t(1)) datenum(t(end))]);
% Update the ylabel
y2 = ylabel('Pore pressure changes (Pa)', 'FontSize', 10, 'Color', [1, 0.5, 0]);
ylim(ax(3),[-450 0]);
set(ax(3), 'ytick', [-450 -200 0], 'yticklabel', {'450', '200', '0'}, 'FontWeight', 'normal', 'FontSize', 10);

%% Panel 4: Long-term Vertical dv/v Comparison
% Find second minimum for long-term analysis
[~, min_index_vel_2] = findSecondMinimum(delta_sigma_vel(:, 1));
scatter(ax(4),datenum(t),new_fitted_dv_v_vel,11,'MarkerEdgeColor','r','MarkerFaceColor','r');
hold(ax(4),'on');
plot(ax(4),datenum(t),fitted_dv_v_vel,'Color',[1 0.7 0.7],'LineWidth',3);
scatter(ax(4),datenum(t),dv_v_syn_vel(:,min_index_vel_2,1),11,'MarkerEdgeColor','k','MarkerFaceColor','k');
box(ax(4),'on')
z_integer = round(z_values(min_index_vel_2));
text_string = sprintf('$\\zeta = %d \\, \\mathrm{(s \\cdot m^{\\frac{1}{2}})}$', z_integer);
text(ax(4),datenum(t(140)),0.2, text_string, 'Interpreter', 'latex', 'FontSize', 18);
xlim(ax(4),[datenum(t(1)) datenum(t(end))]);
set(ax(4),'Xtick',[datenum(t(19)) datenum(t(49)) datenum(t(80)) datenum(t(111)) datenum(t(139))...
    datenum(t(170)) datenum(t(200)) datenum(t(231)) datenum(t(261)) datenum(t(292))]...
    ,'Xticklabel',{},'FontWeight','normal','FontSize',13);
ax(4).TickLength = [0.005 0.01];
ax(4).FontSize = 11;
axes(ax(4));
gx = gca;
gx.XAxisLocation = 'bottom';
gx.XAxis.TickDirection = 'out';
ylabel('dv/v (%)','FontSize',11);
ylim([-0.3 0.3]);
set(ax(4),'Ytick',[-0.2 0 0.2]);
Tdata=t(1:end);
years = year(Tdata);
months = month(Tdata);
firstOccurrence_month = nan(max(years) - min(years) + 1, 12);
for i = min(years):max(years)
    for j = 1:12
        idx = find(years == i & months == j, 1, 'first');
        if ~isempty(idx)
            firstOccurrence_month(i - min(years) + 1, j) = idx;
        end
    end
end
firstOccurrence_month=reshape(firstOccurrence_month',1,[]);
firstOccurrence_month(find(isnan(firstOccurrence_month)))=[];
firstOccurrence_month=firstOccurrence_month+14;
axPos = get(ax(4), 'Position');
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
cx = gca;
cx.XColor = 'white';
cx.YColor = 'white';
line(ax(4), [ax(4).XLim(1) ax(4).XLim(2)], [ax(4).YLim(1) ax(4).YLim(1)], 'Color', 'k'); % 底部
line(ax(4), [ax(4).XLim(1) ax(4).XLim(1)], [ax(4).YLim(1) ax(4).YLim(2)], 'Color', 'k'); % 左侧
line(ax(4), [ax(4).XLim(2) ax(4).XLim(2)], [ax(4).YLim(1) ax(4).YLim(2)], 'Color', 'k'); % 右侧
% Alternative simpler method:
% Set the axes box to 'off' and then set the top color to 'none'
ax(4).Box = 'off';
ax(4).XAxis.Color = 'k';  % Keep the bottom x-axis visible
ax(4).YAxis.Color = 'k';  % Keep the y-axis visible
% Remove the top edge
ax(4).TickDir = 'out';  % Make ticks point outward (optional)

%% Panel 5: Precipitation and Pore Pressure for Long-term Vertical
% Calculate pore pressure for long-term z
z = z_values(min_index_vel_2);
P = zeros(n, 1);
for tt = 2:n
    tmp_P=0;
    for  ii=1:tt
        tmp_P=tmp_P+delta_p(ii) * erfc(z/sqrt(delta_tt*(tt-ii)));
    end
    P(tt)=tmp_P;
end
axes(ax(5));
yyaxis left
bar(ax(5), datenum(t), -Rain_data, 1, 'b');
ylim([-180 0]);
set(ax(5), 'XAxisLocation', 'top', 'FontWeight', 'bold', 'FontSize', 13);
set(ax(5), 'ytick', [-150 -100 -50 0], 'yticklabel', {'150', '100', '50', '0'}, 'FontWeight', 'normal', 'FontSize', 10)
yl = ylabel('Precipitation (mm)','FontSize',11);
set(yl, 'Position', get(yl,'Position') + [0 20 0]);
hold on
% Adjust tick properties
bx = gca;
bx.TickLength = [0.005 0.01];
ax(5).Box = 'off';  % Turn off the box around the plot
% Create custom box with only three sides (left, right, and top)
line(ax(5), [ax(5).XLim(1) ax(5).XLim(2)], [ax(5).YLim(2) ax(5).YLim(2)], 'Color', 'k');
line(ax(5), [ax(5).XLim(1) ax(5).XLim(1)], [ax(5).YLim(1) ax(5).YLim(2)], 'Color', 'k');
line(ax(5), [ax(5).XLim(2) ax(5).XLim(2)], [ax(5).YLim(1) ax(5).YLim(2)], 'Color', 'k');
ax(5).Box = 'off';
% Since your plot has the x-axis at the top, ensure it remains visible
ax(5).XAxis(1).Color = 'k';  % Keep the top x-axis visible
ax(5).YAxis(2).Color = 'k';  % Keep the y-axis visible
hold on
yyaxis right
plot(ax(5), datenum(t), -P(end-321:end),'Color',[1, 0.5, 0], 'LineWidth', 1.5);
ax(5).YAxis(2).Color = [1, 0.5, 0];
% Set y-axis ticks with negative values but positive labels
set(ax(5),'Xtick',[],'Xticklabel',{});
xlim(ax(5),[datenum(t(1)) datenum(t(end))]);
% Update the ylabel
y2 = ylabel('Pore pressure changes (Pa)', 'FontSize', 10, 'Color', [1, 0.5, 0]);
ylim(ax(5),[-35 -20]);
set(ax(5), 'ytick', [-30 -20], 'yticklabel', {'30', '20'}, 'FontWeight', 'normal', 'FontSize', 10)
%% Panels 6-7: Horizontal Component Analysis
% Similar structure to vertical analysis but for horizontal components
[~, min_index_hor] = min(delta_sigma_hor(:, 1));
scatter(ax(6),datenum(t),new_fitted_dv_v_hor,11,'MarkerEdgeColor','g','MarkerFaceColor','g');
hold(ax(6),'on');
plot(ax(6),datenum(t),fitted_dv_v_hor,'Color',[0.5 0.9 0.5],'LineWidth',3);
scatter(ax(6),datenum(t),dv_v_syn_hor(:,min_index_hor,1),11,'MarkerEdgeColor','k','MarkerFaceColor','k');
% box on
z_integer = round(z_values(min_index_hor));
text_string = sprintf('$\\zeta = %d \\, \\mathrm{(s \\cdot m^{\\frac{1}{2}})}$', z_integer);
text(ax(6),datenum(t(140)),0.5, text_string, 'Interpreter', 'latex', 'FontSize', 18);
xlim(ax(6),[datenum(t(1)) datenum(t(end))]);
set(ax(6),'Xtick',[datenum(t(19)) datenum(t(49)) datenum(t(80)) datenum(t(111)) datenum(t(139))...
    datenum(t(170)) datenum(t(200)) datenum(t(231)) datenum(t(258)) datenum(t(261)) datenum(t(264)) datenum(t(292))]...
    ,'Xticklabel',{},'FontWeight','normal','FontSize',13);
ax(6).TickLength = [0.005 0.01];
ax(6).FontSize = 11;
axes(ax(6));
gx = gca;
gx.XAxisLocation = 'bottom';
gx.XAxis.TickDirection = 'out';
ylabel('dv/v (%)','FontSize',11);
ylim([-0.75 0.7]);
set(ax(6),'Ytick',[-0.5 0 0.5]);
xlabelObj = xlabel('Calendar time (month:year)');
set(xlabelObj, 'FontSize', 14);
pos = get(xlabelObj, 'Position');
set(xlabelObj, 'Position', [pos(1), pos(2)-0.3, pos(3)]);
Tdata=t(1:end);
years = year(Tdata);
months = month(Tdata);
firstOccurrence_month = nan(max(years) - min(years) + 1, 12);
for i = min(years):max(years)
    for j = 1:12
        idx = find(years == i & months == j, 1, 'first');
        if ~isempty(idx)
            firstOccurrence_month(i - min(years) + 1, j) = idx;
        end
    end
end
firstOccurrence_month=reshape(firstOccurrence_month',1,[]);
firstOccurrence_month(find(isnan(firstOccurrence_month)))=[];
firstOccurrence_month=firstOccurrence_month+14;
axPos = get(ax(6), 'Position');
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
    text(datenum(Tdata(midpoint(i))), 0.45, num2str(years(i)), 'Parent', ax1, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center','FontSize',12);
end
xlim([datenum(Tdata(1)) datenum(Tdata(end))])
set(gca,'Visible','off')
box off
fig = fg;
fig.Color = 'white';
cx = gca;
cx.XColor = 'white';
cx.YColor = 'white';
box(ax(6),'on');
line(ax(6), [ax(6).XLim(1) ax(6).XLim(2)], [ax(6).YLim(1) ax(6).YLim(1)], 'Color', 'k'); % 底部
line(ax(6), [ax(6).XLim(1) ax(6).XLim(1)], [ax(6).YLim(1) ax(6).YLim(2)], 'Color', 'k'); % 左侧
line(ax(6), [ax(6).XLim(2) ax(6).XLim(2)], [ax(6).YLim(1) ax(6).YLim(2)], 'Color', 'k'); % 右侧
% Alternative simpler method:
% Set the axes box to 'off' and then set the top color to 'none'
ax(6).Box = 'off';
ax(6).XAxis.Color = 'k';  % Keep the bottom x-axis visible
ax(6).YAxis.Color = 'k';  % Keep the y-axis visible
% Remove the top edge
ax(6).TickDir = 'out';  % Make ticks point outward (optional)
%% plot month rain figure
z = z_values(min_index_hor);
P = zeros(n, 1);
for tt = 2:n
    tmp_P=0;
    for  ii=1:tt
        tmp_P=tmp_P+delta_p(ii) * erfc(z/sqrt(delta_tt*(tt-ii)));
    end
    P(tt)=tmp_P;
end
axes(ax(7));
yyaxis left
bar(ax(7), datenum(t), -Rain_data, 1, 'b');
ylim([-180 0]);
set(ax(7), 'XAxisLocation', 'top', 'FontWeight', 'bold', 'FontSize', 13);
% Invert y-axis labels to show positive values
set(ax(7), 'ytick', [-150 -100 -50 0], 'yticklabel', {'150', '100', '50', '0'}, 'FontWeight', 'normal', 'FontSize', 10)
yl = ylabel('Precipitation (mm)','FontSize',11);
set(yl, 'Position', get(yl,'Position') + [0 20 0]);
hold on
% Adjust tick properties
bx = gca;
bx.TickLength = [0.005 0.01];
ax(7).Box = 'off';  % Turn off the box around the plot
% Create custom box with only three sides (left, right, and top)
line(ax(7), [ax(7).XLim(1) ax(7).XLim(2)], [ax(7).YLim(2) ax(7).YLim(2)], 'Color', 'k');
line(ax(7), [ax(7).XLim(1) ax(7).XLim(1)], [ax(7).YLim(1) ax(7).YLim(2)], 'Color', 'k');
line(ax(7), [ax(7).XLim(2) ax(7).XLim(2)], [ax(7).YLim(1) ax(7).YLim(2)], 'Color', 'k');
ax(7).Box = 'off';
% nce your plot has the x-axis at the top, ensure it remains visible
ax(7).XAxis(1).Color = 'k';  % Keep the top x-axis visible
ax(7).YAxis(2).Color = 'k';  % Keep the y-axis visible
hold on

yyaxis right
plot(ax(7), datenum(t), -P(end-321:end),'Color',[1, 0.5, 0], 'LineWidth', 1.5);
ax(7).YAxis(2).Color = [1, 0.5, 0];
% Set y-axis ticks with negative values but positive labels
set(ax(7),'Xtick',[],'Xticklabel',{});
xlim(ax(7),[datenum(t(1)) datenum(t(end))]);
% Update the ylabel
y2 = ylabel('Pore pressure changes (Pa)', 'FontSize', 10, 'Color', [1, 0.5, 0]);
ylim(ax(7),[-450 0]);
set(ax(7), 'ytick', [-450 -200 0], 'yticklabel', {'450', '200', '0'}, 'FontWeight', 'normal', 'FontSize', 10)

%% horizontal
[~, min_index_hor_2] = findSecondMinimum(delta_sigma_hor(:, 1));
scatter(ax(8),datenum(t),new_fitted_dv_v_hor,11,'MarkerEdgeColor','g','MarkerFaceColor','g'); %b
hold(ax(8),'on');
plot(ax(8),datenum(t),fitted_dv_v_hor,'Color',[0.5 0.9 0.5],'LineWidth',3); %b
scatter(ax(8),datenum(t),dv_v_syn_hor(:,min_index_hor_2,1),11,'MarkerEdgeColor','k','MarkerFaceColor','k'); %b
box(ax(8),'on');
z_integer = round(z_values(min_index_hor_2));
text_string = sprintf('$\\zeta = %d \\, \\mathrm{(s \\cdot m^{\\frac{1}{2}})}$', z_integer);
text(ax(8),datenum(t(140)),0.2, text_string, 'Interpreter', 'latex', 'FontSize', 18);
xlim(ax(8),[datenum(t(1)) datenum(t(end))]);
set(ax(8),'Xtick',[datenum(t(19)) datenum(t(49)) datenum(t(80)) datenum(t(111)) datenum(t(139))...
    datenum(t(170)) datenum(t(200)) datenum(t(231)) datenum(t(261)) datenum(t(292))]...
    ,'Xticklabel',{},'FontWeight','normal','FontSize',13);
ax(8).TickLength = [0.005 0.01];
ax(8).FontSize = 11;
axes(ax(8));
gx = gca;
gx.XAxisLocation = 'bottom';
gx.XAxis.TickDirection = 'out';
ylabel('dv/v (%)','FontSize',11);
ylim([-0.3 0.3]);
set(ax(8),'Ytick',[-0.2 0 0.2]);
Tdata=t(1:end);
years = year(Tdata);
months = month(Tdata);
firstOccurrence_month = nan(max(years) - min(years) + 1, 12);
for i = min(years):max(years)
    for j = 1:12
        idx = find(years == i & months == j, 1, 'first');
        if ~isempty(idx)
            firstOccurrence_month(i - min(years) + 1, j) = idx;
        end
    end
end
firstOccurrence_month=reshape(firstOccurrence_month',1,[]);
firstOccurrence_month(find(isnan(firstOccurrence_month)))=[];
firstOccurrence_month=firstOccurrence_month+14;

axPos = get(ax(8), 'Position');
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
    text(datenum(Tdata(midpoint(i))), 0.45, num2str(years(i)), 'Parent', ax1, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center','FontSize',12);
end
xlim([datenum(Tdata(1)) datenum(Tdata(end))])
set(gca,'Visible','off')
box off
fig = fg;
fig.Color = 'white';
cx = gca;
cx.XColor = 'white';
cx.YColor = 'white';
box(ax(8),'on');

line(ax(8), [ax(8).XLim(1) ax(8).XLim(2)], [ax(8).YLim(1) ax(8).YLim(1)], 'Color', 'k'); % 底部
line(ax(8), [ax(8).XLim(1) ax(8).XLim(1)], [ax(8).YLim(1) ax(8).YLim(2)], 'Color', 'k'); % 左侧
line(ax(8), [ax(8).XLim(2) ax(8).XLim(2)], [ax(8).YLim(1) ax(8).YLim(2)], 'Color', 'k'); % 右侧
% Alternative simpler method:
% Set the axes box to 'off' and then set the top color to 'none'
ax(8).Box = 'off';
ax(8).XAxis.Color = 'k';  % Keep the bottom x-axis visible
ax(8).YAxis.Color = 'k';  % Keep the y-axis visible
% Remove the top edge
ax(8).TickDir = 'out';  % Make ticks point outward (optional)

%% plot month rain figure
z = z_values(min_index_hor_2);
P = zeros(n, 1);
for tt = 2:n
    tmp_P=0;
    for  ii=1:tt
        tmp_P=tmp_P+delta_p(ii) * erfc(z/sqrt(delta_tt*(tt-ii)));
    end
    P(tt)=tmp_P;
end
axes(ax(9));

yyaxis left
bar(ax(9), datenum(t), -Rain_data, 1, 'b');
ylim([-180 0]);
set(ax(9), 'XAxisLocation', 'top', 'FontWeight', 'bold', 'FontSize', 13);
% Invert y-axis labels to show positive values
set(ax(9), 'ytick', [-150 -100 -50 0], 'yticklabel', {'150', '100', '50', '0'}, 'FontWeight', 'normal', 'FontSize', 10)
yl = ylabel('Precipitation (mm)','FontSize',11);
set(yl, 'Position', get(yl,'Position') + [0 20 0]);
hold on
% Adjust tick properties
bx = gca;
bx.TickLength = [0.005 0.01];
ax(9).Box = 'off';  % Turn off the box around the plot
% Create custom box with only three sides (left, right, and top)
line(ax(9), [ax(9).XLim(1) ax(9).XLim(2)], [ax(9).YLim(2) ax(9).YLim(2)], 'Color', 'k');
line(ax(9), [ax(9).XLim(1) ax(9).XLim(1)], [ax(9).YLim(1) ax(9).YLim(2)], 'Color', 'k');
line(ax(9), [ax(9).XLim(2) ax(9).XLim(2)], [ax(9).YLim(1) ax(9).YLim(2)], 'Color', 'k');
ax(9).Box = 'off';
% Since your plot has the x-axis at the top, ensure it remains visible
ax(9).XAxis(1).Color = 'k';  % Keep the top x-axis visible
ax(9).YAxis(2).Color = 'k';  % Keep the y-axis visible
hold on

yyaxis right
plot(ax(9), datenum(t), -P(end-321:end),'Color',[1, 0.5, 0], 'LineWidth', 1.5);
ax(9).YAxis(2).Color = [1, 0.5, 0];
% Set y-axis ticks with negative values but positive labels
set(ax(9),'Xtick',[],'Xticklabel',{});
xlim(ax(9),[datenum(t(1)) datenum(t(end))]);
% Update the ylabel
y2 = ylabel('Pore pressure changes (Pa)', 'FontSize', 10, 'Color', [1, 0.5, 0]);
ylim(ax(9),[-35 -20]);
set(ax(9), 'ytick', [-30 -20], 'yticklabel', {'30', '20'}, 'FontWeight', 'normal', 'FontSize', 10)


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