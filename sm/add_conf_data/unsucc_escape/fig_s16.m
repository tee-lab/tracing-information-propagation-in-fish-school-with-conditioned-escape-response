close all
clear
clc

%% load data

all_par_cross_data = load("partial_cross_data.mat"); % load all pos, vel data
load("partial_cross_data.mat", "dt", "n", "no_exp", "tank_max_x", "tank_min_x", "tank_middle", "del_step")

t_cut_min = -0.5;
t_cut_max = 10;

ini_color = "#33b1ff";
esp_color = "#fa4d56";
relax_color = "#198038";
speed_plt_color = ["#6929c4", "#1192e8", "#005d5d", "#9f1853", "#d2a106"];
font_size = 30;
lw_plot = 2;
lw_axis = 2;
lw_xline = 4;
tl_major = 0.02;
tl_minor = 0.0005;

%% storing all data

speed_all = []; % speed 
speed_ini_all = []; % store all speeds in initial phase
speed_escape_all = []; % store all speed in escape phase.
speed_relax_all = []; % store all speeds in relax phase

%% Group properties analysis

for e = 1:no_exp
    
    % fish data for this given experiment
    pos_fish = all_par_cross_data.(strcat('pos_t_ex_', num2str(e))); % input what pos data to use. 
    vel_fish = all_par_cross_data.(strcat('vel_t_ex_', num2str(e))); % input what vel data to use
    speed_fish = all_par_cross_data.(strcat('speed_t_ex_', num2str(e))); % input what speed data to use
    time_light_on = all_par_cross_data.(strcat('frame_light_on_ex_', num2str(e))); % when was the light turned on
    time_light_on = time_light_on*dt; % in seconds

    vel_fish = vel_fish(:,:,1:end-del_step);
    speed_fish = speed_fish(:,1:end-del_step);
    
    tot_time = size(vel_fish,3); % total number of frames
    t_st = 1; % frame start
    t_skip = 1; % frames to skip
    t_end = tot_time; % last frame
    t_plt = (t_st:t_skip:t_end)*dt; % time in sec

    % identifying when the conditioned fish first crosses the boarder
    pos_ini_x = pos_fish(:,1,1); % position of fish at t = 0
    min_pos_ini_x = min(pos_ini_x); % min pos at t = 0
    max_pos_ini_x = max(pos_ini_x); % max pos at t = 0
    pos_x = squeeze(pos_fish(:,1,:)); % x position of all fish across all times
    % (as y position doesnt matter, only when they crossed borders matters)

    if max_pos_ini_x < tank_middle % if max initial x pos should also be < tank border

        t_cross_min = pos_x > tank_min_x; % time at which fish cross the first border
       
        fcross_tmin = nan(n,1); % store the time at which fish 1st cross the 1st border

        for i = 1:n
            
            idx = [find(t_cross_min(i,:), 1), NaN];
            fcross_tmin(i) = idx(1);

        end

    elseif min_pos_ini_x > tank_middle % if even min initial x pos also be > tank border

        t_cross_min = pos_x < tank_max_x;

        fcross_tmin = nan(n,1);

        for i = 1:n
            
            idx = [find(t_cross_min(i,:), 1), NaN];
            fcross_tmin(i) = idx(1);

        end

    end

    [~, rank_id] = sort(fcross_tmin, 'ascend');

    t_plt = t_plt - time_light_on;

    % speed analysis
    speed_rank = speed_fish(rank_id,:);
    
    spd_all_temp = speed_rank(:,t_plt >= t_cut_min & t_plt <= t_cut_max);
    spd_all_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', spd_all_temp'];
    speed_all = cat(1, speed_all, spd_all_temp);

end

%% plotting average plots

% plotting all speed

no_time_edges = 51;

plt_count = 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

min_spd = min(min(speed_all(:,2:end)));
max_spd = max(max(speed_all(:,2:end)));

[~, spd_sort_id] = sort(speed_all(:,1));
speed_all = speed_all(spd_sort_id,:);

spd_cf_ts_ideal = readmatrix('spd_cf_ts_mt.csv');
errorbar(spd_cf_ts_ideal(:,1), spd_cf_ts_ideal(:,2), spd_cf_ts_ideal(:,3), ...
    '-o', 'LineWidth', 2, 'Color', speed_plt_color(2))

hold on  

spd_cf_ts_add = readmatrix('spd_cf_ts_all_cross_conf.csv');
errorbar(spd_cf_ts_add(:,1), spd_cf_ts_add(:,2), spd_cf_ts_add(:,3), ...
    '-o', 'LineWidth', 2, 'Color', speed_plt_color(3))

hold on

for i = 2:2    

    t_edges = linspace(t_cut_min, t_cut_max, no_time_edges);
    spd_edges = linspace(min_spd, max_spd, 31);
    [histcount_spd, ~, ~, t_bins, spd_bins] = histcounts2(speed_all(:,1), ...
        speed_all(:,i), t_edges, spd_edges);

    t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
    spd_edges = spd_edges(1:end-1) + (spd_edges(2) - spd_edges(1))/2;

    mean_spd_i = nan(1,size(histcount_spd,1));
    std_spd_i = nan(1,size(histcount_spd,1));
    se_spd_i = nan(1,size(histcount_spd,1));
    median_spd_i = nan(1,size(histcount_spd,1));
    for j = 1:size(histcount_spd,1)
        t_id = t_bins == j;
        spd_temp = speed_all(t_id,i);
        mean_spd_i(j) = mean(spd_temp);
        std_spd_i(j) = std(spd_temp);
        se_spd_i(j) = std(spd_temp)/sqrt(length(spd_temp));
        median_spd_i(j) = median(spd_temp);
    end

    % plot(t_edges, mean_spd_i, '-o')
    errorbar(t_edges, mean_spd_i, se_spd_i, '-o', 'LineWidth', 2, ...
        'Color', speed_plt_color(i-1))
    hold on

end

xline(0, '--r', 'LineWidth', lw_xline)
hold off

legend({'All fish cross (considered, N = 53)', 'All fish cross (confounded, N = 25)', ...
    'Partial crossing (N = 22)'}, ...
    'Location', 'southeast', 'FontSize', 23)
legend('box','off')

set(gca, 'XLim', [-0.5 5], 'YLim', [0,16], 'YTick', 0:3:20,...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [tl_major, tl_minor],...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')
ax = gca;
ax.XAxis.MinorTickValues = -1:0.5:10;  % Minor X-ticks every 0.1
ax.YAxis.MinorTickValues = 0:1.5:30;

xlabel('$\tilde{t} - \tilde{t}_0$', 'Interpreter', 'latex')
ylabel('Mean speed (cm/s)')

% exportgraphics(gca, 'speed_ts_comparision.pdf', 'ContentType', 'vector')