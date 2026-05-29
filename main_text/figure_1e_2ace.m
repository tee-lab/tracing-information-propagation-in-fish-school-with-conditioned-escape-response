close all
clear
clc

%%

% load one file to find the length of the video

file_name = dir('*.csv');
file_name = file_name.name;

tot_time = readmatrix(file_name);
tot_time = tot_time(end,1) + 1; % because time starts from t = 0. 

file_name = file_name(1:end-5);

dt = 0.04; % del t between 2 frames. 
n = 5; % no.of fish
col_posx = 2; % column with x coordinates
col_posy = 3; % column with y coordinates
smooth_window = 6;

t_st = 1;
t_skip = 1;
t_end = tot_time;
t_plt = (t_st:t_skip:t_end)*dt;

time_light_on = 53*dt; % time when light was turned on

time_scale = 2; % 1 to consider t = 0 when first fish crosses the border 
% and 2 to consider t = 0 when light is turned on.

% tank compartment distance
tank_min_x = 22.5; % 21.8. 
tank_max_x = 28.5; % 29.8;
tank_middle = (tank_max_x + tank_min_x)/2;

speed_plt_color = ["#6929c4", "#1192e8", "#005d5d", "#9f1853", "#d2a106"];

two_d_color = "#A52A2A";
ini_color = "#0096FF";
relax_color = "#097969";
font_size = 25;
lw_plot = 2;
lw_axis = 2;
lw_xline = 4;
x_min_lim = -0.5;
x_max_lim = 3;

%% collecting positions of all fish from tracked data files

pos_t = nan(n,2,tot_time); % store position data

% intercept if there are missing data. ideally we shouldnt have any missing
% data in the new tracks by Fathimath.
pos_interp = nan(n,2,tot_time); 
t_interp = 1:tot_time;
mthd = 'spline'; % method to interpolate.

for i = 1:n 
    
    fname = strcat(file_name,num2str(i),'.csv'); % name of the file for 'i'th fish
    fno = readmatrix(fname); % load the file
    tot_time_temp = fno(end,1); % no.of frames (or time)
    
    % in case time is greater than the total time we loaded above, set the
    % tot_time to new tot_time. and pos_t with new length. 
    if tot_time_temp > tot_time    
        pos_t_temp = pos_t;
        pos_t = nan(n,2,tot_time_temp);
        pos_t(:,:,1:tot_time) = pos_t_temp;
        tot_time = tot_time_temp;  
    end
    
    % if the first frame for this fish is not t = 0, then fill NaN till
    % t = T - 1, where T is the first frame for this fish
    % ideally this loop shouldnt run for the new tracked videos by
    % Fathimath. 
    if fno(1,1) > 0 
        fno_temp = nan(size(fno,1)+fno(1,1), size(fno,2));
        fno_temp(:,1) = 0:(tot_time-1);
        fno_temp(fno(1,1)+1:end,:) = fno;
        fno = fno_temp;
    end
    
    pos_x = fno(:,col_posx); % x coordinates for _i_ th fish
    pos_y = fno(:,col_posy); % y coordinates for _i_ th fish
    inf_id = (isinf(pos_x) | isinf(pos_y)); % identify Inf in either x or y coordinates (again should be empty set for new tracked videos) 
    pos_x(inf_id) = nan;
    pos_y(inf_id) = nan;
    
    % store position
    pos_t(i,1,1:tot_time) = pos_x;
    pos_t(i,2,1:tot_time) = pos_y;

    t_interp = 1:tot_time;
    
    % again we there shouldnt be missing values in beedance tracks
    pos_interp_x = squeeze(pos_t(i,1,:));
    pos_interp_y = squeeze(pos_t(i,2,:));
    xq_nan = t_interp(isnan(pos_interp_x));
    yq_nan = t_interp(isnan(pos_interp_y));
    pos_interp_x(xq_nan) = interp1(t_interp(~isnan(pos_interp_x))*dt, pos_interp_x(~isnan(pos_interp_x)), xq_nan*dt, mthd);
    pos_interp_y(yq_nan) = interp1(t_interp(~isnan(pos_interp_y))*dt, pos_interp_y(~isnan(pos_interp_y)), yq_nan*dt, mthd);
    pos_interp(i,1,:) = pos_interp_x;
    pos_interp(i,2,:) = pos_interp_y;

end

pos_t(:,1,:) = smoothdata(pos_t(:,1,:), 3, 'gaussian', smooth_window);
pos_t(:,2,:) = smoothdata(pos_t(:,2,:), 3, 'gaussian', smooth_window);

% calculate and store velocity
del_step = 3;
vel_t = nan(n,2,tot_time);
vel_x = (pos_t(:,1,del_step+1:end) - pos_t(:,1,1:(end-del_step)))/(dt*del_step); % vel_x
vel_y = (pos_t(:,2,del_step+1:end) - pos_t(:,2,1:(end-del_step)))/(dt*del_step); % vel_y
vel_t(:,1,1:end-del_step) = vel_x;
vel_t(:,2,1:end-del_step) = vel_y;

% calculate speed.
speed_t = squeeze(vecnorm(vel_t,2,2));

% this should be same as vel_t and speed_t for beedance tracks
vel_interp = nan(n,2,tot_time);
vel_x_interp = (pos_interp(:,1,2:end) - pos_interp(:,1,1:(end-1)))/dt;
vel_y_interp = (pos_interp(:,2,2:end) - pos_interp(:,2,1:(end-1)))/dt;
vel_interp(:,1,2:end) = vel_x_interp;
vel_interp(:,2,2:end) = vel_y_interp;
speed_interp = squeeze(vecnorm(vel_interp,2,2));

%% identifying when the conditioned fish first crosses the boarder

pos_fish = pos_t; % input what pos data to use. 'pos_t', 'pos_interp'
vel_fish = vel_t; % input what vel data to use: 'vel_t', 'vel_interp'
speed_fish = speed_t; % input what speed data to use: 'speed_cal', 'speed_interp'

pos_ini_x = pos_fish(:,1,1); % position of fish at t = 0
min_pos_ini_x = min(pos_ini_x, [], 'omitmissing'); % min pos at t = 0
max_pos_ini_x = max(pos_ini_x, [], 'omitmissing'); % max pos at t = 0
pos_x = squeeze(pos_fish(:,1,:)); % x position of all fish across all times 
% (as y position doesnt matter, only when they crossed borders matters)

if max_pos_ini_x < tank_middle % if max initial x pos should also be < tank border

    t_cross_min = pos_x > tank_min_x; % time at which fish cross the first border
    t_cross_max = pos_x > tank_max_x; % time at which fish cross the second border
    t_cross_centre = pos_x > tank_middle; % time at which fish cross the middle of the box
        
    fcross_tmin = nan(n,1); % store the time at which fish 1st cross the 1st border
    fcross_tmax = nan(n,1); % 2nd border 
    fcross_tcentre = nan(n,1); % centre of the box

    for i = 1:n

        fcross_tmin(i) = find(t_cross_min(i,:), 1);
        fcross_tmax(i) = find(t_cross_max(i,:), 1);
        fcross_tcentre(i) = find(t_cross_centre(i,:),1);

    end

elseif min_pos_ini_x > tank_middle % if even min initial x pos also be > tank border

    t_cross_min = pos_x < tank_max_x;
    t_cross_max = pos_x < tank_min_x;
    t_cross_centre = pos_x < tank_middle;
        
    fcross_tmin = nan(n,1);
    fcross_tmax = nan(n,1);
    fcross_tcentre = nan(n,1);

    for i = 1:n

        fcross_tmin(i) = find(t_cross_min(i,:), 1);
        fcross_tmax(i) = find(t_cross_max(i,:), 1);
        fcross_tcentre(i) = find(t_cross_centre(i,:),1);

    end

end

% time when 1st fish crossed the 1st and 2nd border, and centre of the box
tcross_min = min(fcross_tmin);
tcross_max = max(fcross_tmax);
tcross_centre_min = min(fcross_tcentre)*dt;
tcross_centre_max = max(fcross_tcentre)*dt;

if time_scale == 1

    t_plt = (t_plt - tcross_centre_min)/(tcross_centre_max - tcross_centre_min);

elseif time_scale == 2

    t_plt = (t_plt - time_light_on)/(tcross_centre_max - time_light_on);

end

[~, rank_id] = sort(fcross_tmin, 'ascend');
[~, rank_id_centre] = sort(fcross_tcentre, 'ascend');

tcross_min_diff = fcross_tmin(rank_id(2:end)) - fcross_tmin(rank_id(1:end-1));
tcross_max_diff = fcross_tmax(rank_id(2:end)) - fcross_tmax(rank_id(1:end-1));
tcross_centre_diff = fcross_tcentre(rank_id_centre(2:end)) - fcross_tcentre(rank_id_centre(1:end-1));

%% group polarisation and cohesion (dispersion)

vel_norm = vel_fish./(vecnorm(vel_fish,2,2)+eps);
vn_x = squeeze(vel_norm(:,1,:)); % vx
vn_y = squeeze(vel_norm(:,2,:)); % vy
mx = mean(vn_x,1,'omitmissing'); % mx
my = mean(vn_y,1,'omitmissing'); % my
m = sqrt(mx.^2 + my.^2); % m

m = smoothdata(m, 'gaussian', smooth_window);
mx = smoothdata(abs(mx), 'gaussian', smooth_window);
my = smoothdata(abs(my), 'gaussian', smooth_window);

plt_count = 0;

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(t_plt, m(t_st:t_skip:t_end), 'LineWidth', lw_plot, 'Color', two_d_color)
hold on
plot(t_plt, mx(t_st:t_skip:t_end), 'LineWidth', lw_plot, 'Color', ini_color)
hold on
plot(t_plt, my(t_st:t_skip:t_end), 'LineWidth', lw_plot, 'Color', relax_color)
hold on
xline(0, '--r', 'LineWidth', lw_xline)
hold on
xline(1, '--r', 'LineWidth', lw_xline)
hold off

% legend({'P', 'P_x', 'P_y'}, 'Location', 'best')
% legend('boxoff')

set(gca, 'XLim', [x_min_lim x_max_lim], 'YLim', [0, 1], 'YTick', 0:.2:1, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('')
ylabel('Polarisation')

% exportgraphics(gca, 'pol_ts_samp_data.pdf', 'ContentType', 'vector')

gc_x = mean(pos_fish(:,1,:), 1, 'omitmissing'); % group centre - x 
gc_y = mean(pos_fish(:,2,:), 1, 'omitmissing'); % group centre - y

% position in group centre reference frame
pos_gc(:,1,:) = pos_fish(:,1,:) - gc_x;
pos_gc(:,2,:) = pos_fish(:,2,:) - gc_y;
dist_gc = vecnorm(pos_gc,2,2);
mean_dist_gc = squeeze(mean(dist_gc,1,'omitmissing'));

% smoothening the data
mean_dist_gc = smoothdata(mean_dist_gc, 'gaussian', smooth_window);

% dispersion only along x axis

pos_gc_x = squeeze(pos_gc(:,1,:));
pos_gc_x = abs(pos_gc_x);
mean_gc_x = mean(pos_gc_x, 1, 'omitmissing');
mean_gc_x = smoothdata(mean_gc_x, 'gaussian', smooth_window);

% dispersion along y axis
pos_gc_y = squeeze(pos_gc(:,2,:));
pos_gc_y = abs(pos_gc_y);
mean_gc_y = mean(pos_gc_y, 1);
mean_gc_y = smoothdata(mean_gc_y, 'gaussian', smooth_window);

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(t_plt, mean_dist_gc(t_st:t_skip:t_end), 'LineWidth', lw_xline, 'Color', two_d_color)
hold on
plot(t_plt, mean_gc_x(t_st:t_skip:t_end), 'LineWidth', lw_xline, 'Color', ini_color)
hold on
plot(t_plt, mean_gc_y(t_st:t_skip:t_end), 'LineWidth', lw_xline, 'Color', relax_color)
hold on
xline(0, '--r', 'LineWidth', lw_xline)
hold on
xline(1, '--r', 'LineWidth', lw_xline)
hold off

set(gca, 'XLim', [x_min_lim x_max_lim], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('')
ylabel('Dispersion (cm)')
% legend({'D', 'D_x', 'D_y'}, 'Location', 'best')
% legend('boxoff')

% exportgraphics(gca, 'disp_ts_sample_data.pdf', 'ContentType', 'vector')

%% speed time series

speed_rank = speed_fish(rank_id,:);

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];
t = tiledlayout(5, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

for i = 1:n
    
    spd_plt_temp = speed_rank(i,:);
    spd_plt_temp = smoothdata(spd_plt_temp, "gaussian", 10);
    
    nexttile
    plot(t_plt, spd_plt_temp(t_st:t_skip:t_end), 'LineWidth', lw_plot, ...
        'Color', speed_plt_color(i))
    hold on
    xline(0, '--r', 'LineWidth', lw_xline)
    hold on
    xline(1, '--r', 'LineWidth', lw_xline)
    hold off

    % legend({strcat("CR - ", num2str(i))}, 'Location', 'best', 'FontSize', 20)
    % legend('Box', 'off')

    set(gca, 'XLim', [x_min_lim x_max_lim], 'YLim', [0, 8], 'YTick', 0:4:8, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', 15, 'FontName', 'Helvetica')

    if i < n
        % Remove x-tick labels for all but the bottom plot
        set(gca, 'XTickLabel', []);
    end

    % if i > 1
    %     % For all plots EXCEPT the top one, hide the top label (8)
    %     yticklabels({'0', '4', ''});
    % else
    %     % Keep all labels for the top plot
    %     yticklabels({'0', '4', '8'});
    % end

end

xlabel(t, '', 'FontSize', font_size, 'FontName', 'Helvetica')
ylabel(t, 'Speed (cm/s)', 'FontSize', font_size, 'FontName', 'Helvetica')

% exportgraphics(fig, 'speed_ts_sample_data.pdf', 'ContentType', 'vector')

%% x time series

pos_rank = pos_t(rank_id,:,:);

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

for i = 1:n
    
    pos_plt_temp = squeeze(pos_rank(i,1,:));
    pos_plt_temp = smoothdata(pos_plt_temp, "gaussian", 10);
    
    plot(t_plt, pos_plt_temp(t_st:t_skip:t_end), 'LineWidth', lw_plot)
    hold on

end

yline(tank_min_x, '--k', 'LineWidth', lw_xline)
hold on
yline(tank_max_x, '--k', 'LineWidth', lw_xline)
hold on
xline(0, '--r', 'LineWidth', lw_xline)
hold on
xline(1, '--r', 'LineWidth', lw_xline)
hold off

% legend({'CR-1', 'CR-2', 'CR-3', 'CR-4', 'CR-5'}, 'Location', 'best')
% legend('Box', 'off')

set(gca, 'XLim', [x_min_lim x_max_lim], 'XTick', -1:1:x_max_lim, ...
    'YLim', [0, 50], 'YTick', 0:10:50, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('t_n')
ylabel('x (cm)')

% exportgraphics(fig, 'pos_x_ts_sample_data.pdf', 'ContentType', 'vector')