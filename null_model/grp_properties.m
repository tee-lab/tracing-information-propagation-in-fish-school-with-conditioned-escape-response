close all
clear
clc

%% Measuring group properties like polarisation, cohesion, etc.

% tic

fname = dir('*.mat');
fname = fname.name;
load(fname) % load data
fname = fname(1:end-4);

%%

ini_color = "#33b1ff";
esp_color = "#fa4d56";
relax_color = "#198038";
speed_plt_color = ["#6929c4", "#1192e8", "#005d5d", "#9f1853", "#d2a106"];

font_size = 25;
lw_axis = 2;
lw_plot = 2;
lw_xline = 3;
mr_size = 10;
label_fs = 25;

%%

no_it = no_it*no_exp; % total no.of iterations (to calculate group properties)
en_start = reshape(en_start, n, no_it); % start time
en_centre = reshape(en_centre, n, no_it); % centre crossing time
en_end = reshape(en_end, n, no_it); % crossing time of 2nd border.
rank_order_atk = reshape(rank_order_atk, n, no_it); % rank of crossing
pos_t = reshape(pos_t, n, 2, n_iter, no_it); % position data
s_t = reshape(s_t, n, n_iter, no_it); % speed data
theta_t = reshape(theta_t, n, n_iter, no_it); % orientation data
cg_time = 5; % smoothing window

%% Encounter time

% Escape phase time

esp_time_st = zeros(no_it,1); % time at which fish crosses the first barrier
esp_time_end = zeros(no_it,1); % time at which fish crosses the second barrier

for i = 1:no_it
    if min(en_start(:,i)) > 0
        esp_time_st(i) = min(en_start(:,i))*dt;
        esp_time_end(i) = max(en_end(:,i))*dt;
    else
        esp_time_st(i) = 0;
        esp_time_end(i) = 0;
    end
end

esp_time = esp_time_end - esp_time_st;

% en_start_min = min(en_start,[],1);
en_end_max = max(en_end,[],1); % time at which the last fish crosses the last barrier

t_plt = (1:n_iter);
t_plt = repmat(t_plt', no_it, 1);
% en_st_time_rep = repmat(en_start_min, n_iter, 1);
% en_st_time_rep = en_st_time_rep(:);
en_st_time_rep = min(t_atk);
en_end_time_rep = repmat(en_end_max, n_iter, 1);
en_end_time_rep = en_end_time_rep(:);
t_plt = (t_plt - en_st_time_rep)./(en_end_time_rep - en_st_time_rep); % normalised time
% 0 = time of light on
% 1 = time at which the last fish crosses the barrier. 

%% Grp polarisation

grp_pol = zeros(n_iter, no_it);

for i = 1:no_it

    theta_temp = squeeze(theta_t(:,:,i));
    vel_x = mean(cos(theta_temp),1);
    vel_y = mean(sin(theta_temp),1);
    pol_tem = sqrt(vel_x.^2 + vel_y.^2);
    grp_pol(:,i) = pol_tem.';

end

pol = [t_plt, grp_pol(:)];

no_pol_edges = 21;
no_time_edges = 201;
[~, pol_sort_id] = sort(pol(:,1));
pol = pol(pol_sort_id, :);

t_edges = linspace(min(t_plt), max(t_plt), no_time_edges);
pol_edges = linspace(0, 1, no_pol_edges);
[histcount_pol, ~, ~, t_bins, pol_bins] = histcounts2(pol(:,1), pol(:,2), t_edges, pol_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
pol_edges = pol_edges(1:end-1) + (pol_edges(2) - pol_edges(1))/2;

mean_pol = nan(1,size(histcount_pol,1));
std_pol = nan(1,size(histcount_pol,1));
se_pol = nan(1,size(histcount_pol,1));
median_pol = nan(1,size(histcount_pol,1));

for i = 1:size(histcount_pol,1)
    t_id = t_bins == i;
    pol_temp = pol(t_id,2);
    mean_pol(i) = mean(pol_temp);
    std_pol(i) = std(pol_temp);
    se_pol(i) = std(pol_temp)/sqrt(length(pol_temp));
    median_pol(i) = median(pol_temp);
end

plt_count = 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(t_edges, mean_pol, 'o-', 'Color', ini_color, 'LineWidth', lw_plot, ...
    'MarkerFaceColor', ini_color)
hold on
xline(0, '--r', 'LineWidth', lw_plot)
hold on
xline(1, '--r', 'LineWidth', lw_plot)

set(gca, 'XLim', [-0.5 3], 'YLim', [0,1], 'YTick', 0:0.2:1, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('t_n', 'FontSize', label_fs)
ylabel('Polarisation', 'FontSize', label_fs)

%% Group speed dynamics

median_grp_spd = zeros(n_iter, no_it);
speed_rank = nan(size(s_t));

for i = 1:no_it
    
    speed_rank(:,:,i) = s_t(rank_order_atk(:,i),:,i);
    spd_temp = squeeze(s_t(rank_order_atk(:,i),:,i));
    median_grp_spd_temp = median(spd_temp,1);
    median_grp_spd(:,i) = median_grp_spd_temp.';

end

speed_rank = reshape(speed_rank, n, no_it*n_iter);
speed_rank = speed_rank';
speed_all = [t_plt, speed_rank];

min_spd = min(min(speed_all(:,2:end)));
max_spd = max(max(speed_all(:,2:end)));

[~, spd_sort_id] = sort(speed_all(:,1));
speed_all = speed_all(spd_sort_id, :);

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

mean_spd_data = [];
se_spd_data = [];
std_spd_data = [];
median_spd_data = [];

for i = 2:n+1

    t_edges = linspace(min(t_plt), max(t_plt), no_time_edges);
    spd_edges = linspace(min_spd, max_spd, 31);
    [histcount_spd, ~, ~, t_bins, spd_bins] = histcounts2(speed_all(:,1), speed_all(:,i), t_edges, spd_edges);

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

    mean_spd_data = cat(2, mean_spd_data, mean_spd_i');
    median_spd_data = cat(2, median_spd_data, median_spd_i');
    std_spd_data = cat(2, std_spd_data, std_spd_i');
    se_spd_data = cat(2, se_spd_data, se_spd_i');

    plot(t_edges, mean_spd_i, '-o', 'LineWidth', lw_plot, ...
        'Color', speed_plt_color(i-1))
    hold on

end

xline(0, '--r', 'LineWidth', lw_plot)
hold on
xline(1, '--r', 'LineWidth', lw_plot)

legend({'CR-1', 'CR-2', 'CR-3', 'CR-4', 'CR-5'}, 'Location', 'best')
legend('box','off')
set(gca, 'XLim', [-0.5 2], 'YLim', [0, 5], 'YTick', 0:.5:5, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('t_n', 'FontSize', label_fs)
ylabel('Speed (cm/s)', 'FontSize', label_fs)

%% Group cohesion

grp_cohesion = zeros(n_iter, no_it);

for i = 1:no_it

    pos_it = pos_t(:,:,:,i);
        
    gc_x = mean(pos_it(:,1,:), 1);
    gc_y = mean(pos_it(:,2,:), 1);

    pos_gc(:,1,:) = pos_it(:,1,:) - gc_x;
    pos_gc(:,2,:) = pos_it(:,2,:) - gc_y;

    dist_gc = vecnorm(pos_gc,2,2);
    mean_dist_gc = squeeze(mean(dist_gc,1));

    grp_cohesion(:,i) = mean_dist_gc;


end

grp_coh = [t_plt, grp_cohesion(:)];

[~, gc_sort_id] = sort(grp_coh(:,1));
grp_coh = grp_coh(gc_sort_id, :);

t_edges = linspace(min(t_plt), max(t_plt), no_time_edges);
min_gc = min(grp_coh(:,2));
max_gc = max(grp_coh(:,2));
gc_edges = linspace(min_gc, max_gc, 31);
[histcount_gc, ~, ~, t_bins, gc_bins] = histcounts2(grp_coh(:,1), grp_coh(:,2), t_edges, gc_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
gc_edges = gc_edges(1:end-1) + (gc_edges(2) - gc_edges(1))/2;

mean_gc = nan(1,size(histcount_gc,1));
std_gc = nan(1,size(histcount_gc,1));
se_gc = nan(1,size(histcount_gc,1));
median_gc = nan(1,size(histcount_gc,1));
for i = 1:size(histcount_gc,1)
    t_id = t_bins == i;
    gc_temp = grp_coh(t_id,2);
    mean_gc(i) = mean(gc_temp);
    std_gc(i) = std(gc_temp);
    se_gc(i) = std(gc_temp)/sqrt(length(gc_temp));
    median_gc(i) = median(gc_temp);
end

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(t_edges, mean_gc, 'o-', 'Color', '#A52A2A', 'LineWidth', lw_plot, ...
    'MarkerFaceColor', '#A52A2A')

set(gca, 'XLim', [-0.5 3], 'YLim', [0, 10], 'LineWidth', lw_axis, ...
    'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('t_n', 'FontSize', label_fs)
ylabel('Dispersion (cm)', 'FontSize', label_fs)

%% time crossing

[fcross_tmin, rank_id_min] = sort(en_start, 1, 'ascend'); % arranging start in ascending order.
% [fcross_tcentre, rank_id_centre] = sort(en_centre, 1, 'ascend');
% [fcross_tmax, rank_id_max] = sort(en_end, 1, 'ascend');

c_id = repmat((2:n)', no_it, 1);
tcross_min_diff = fcross_tmin(2:end,:) - fcross_tmin(1:(end-1),:);
tmind_all = [c_id, tcross_min_diff(:)];
% tcross_centre_diff = fcross_tcentre(2:end,:) - fcross_tcentre(1:(end-1),:);
% tcd_all = [c_id, tcross_centre_diff(:)];
% tcross_max_diff = fcross_tmax(2:end,:) - fcross_tmax(1:(end-1),:);
% tmaxd_all = [c_id, tcross_max_diff(:)];

%  Time difference plots

mean_td_min = nan(1,n-1);
std_td_min = nan(1,n-1);
se_td_min = nan(1,n-1);
median_td_min = nan(1,n-1);
en_time_diff = nan(n, no_it);

for i = 2:n
    td_id = tmind_all(:,1) == i;
    td_temp = tmind_all(td_id,2);
    en_time_diff(i,:) = td_temp;
    mean_td_min(i-1) = mean(td_temp);
    std_td_min(i-1) = std(td_temp);
    se_td_min(i-1) = (1.96*std(td_temp))/sqrt(length(td_temp));
    median_td_min(i-1) = median(td_temp);
end

% mean_td_centre = nan(1,n-1);
% std_td_centre = nan(1,n-1);
% se_td_centre = nan(1,n-1);
% median_td_centre = nan(1,n-1);

% for i = 2:n
%     td_id = tcd_all(:,1) == i;
%     td_temp = tcd_all(td_id,2);
%     mean_td_centre(i-1) = mean(td_temp);
%     std_td_centre(i-1) = std(td_temp);
%     se_td_centre(i-1) = (1.96*std(td_temp))/sqrt(length(td_temp));
%     median_td_centre(i-1) = median(td_temp);
% end
% 
% mean_td_max = nan(1,n-1);
% std_td_max = nan(1,n-1);
% se_td_max = nan(1,n-1);
% median_td_max = nan(1,n-1);
% 
% for i = 2:n
%     td_id = tmaxd_all(:,1) == i;
%     td_temp = tmaxd_all(td_id,2);
%     mean_td_max(i-1) = mean(td_temp);
%     std_td_max(i-1) = std(td_temp);
%     se_td_max(i-1) = (1.96*std(td_temp))/sqrt(length(td_temp));
%     median_td_max(i-1) = median(td_temp);
% end

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

st_id = 2;

errorbar((st_id:n)-0.1, mean_td_min*dt, se_td_min*dt, 'o', ...
        'LineWidth', lw_plot, 'Color', ini_color, 'MarkerSize', mr_size, ...
        'MarkerFaceColor', ini_color)
% hold on
% errorbar((2:n), mean_td_centre, se_td_centre, 'o', 'LineWidth', 1, 'Color', "#4e91fd")
% hold on
% errorbar((2:n)+0.1, mean_td_max, se_td_max, 'o', 'LineWidth', 1, 'Color', "magenta")

% legend({'barrier-start', 'barrier-centre', 'barrier-last'}, 'Location', 'best')
% legend('Box', 'off')

set(gca, 'XLim', [st_id-0.4 n+0.2], 'XTick', 1:n, 'YLim', [0,3], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', label_fs)
ylabel('Time since previous fish crossed (s)', 'FontSize', label_fs)

%%

dist_min_wall = nan(n, no_it); % distance to the closest wall (only along x axis)
dist_cent_wall = nan(n, no_it); % distance to centre of closest wall (2-D)
dist_to_conditioned_fish_x = nan(n-1, no_it); % distance to conditioned fish (x-asis only)
dist_to_conditioned_fish = nan(n-1, no_it); % distance to conditioned fish

for e = 1:no_it
    
    % fish data for this given experiment
    pos_it = pos_t(:,:,:,e); % input what pos data to use.
    time_light_on = min(t_atk); % when was the light turned on
    rank_id_min = rank_order_atk(:,e); % rank of crossing
    
    % position when light was turned on in the order of crossing
    pos_fish_light_on = pos_it(rank_id_min,:,time_light_on);
    % distance to conditioned fish
    dist_to_cond = pos_fish_light_on(2:end,:) - pos_fish_light_on(1,:);
    dist_to_conditioned_fish_x(:,e) = abs(dist_to_cond(:,1));
    dist_to_conditioned_fish(:,e) = vecnorm(dist_to_cond,2,2);
        
    % check which side of the tank fish are in
    if pos_fish_light_on(:,1) < (box_len/2)

        dist_x = mini_box_len - pos_fish_light_on(:,1); % x_min_wall - x

    elseif pos_fish_light_on(:,1) > (box_len/2)

        dist_x = pos_fish_light_on(:,1) - (box_len - mini_box_len); % x - x_max_wall

    end

    dist_y = (box_width/2) - pos_fish_light_on(:,2); % distance to centre of the tank

    dist_min_wall(:,e) = dist_x;
    dist_cent_wall(:,e) = sqrt(dist_x.^2 + dist_y.^2);

end

%% distribution of distance to conditioned fish

mean_dist_to_cond = mean(dist_to_conditioned_fish, 2); % mean dist to conditioned fish
std_dist_to_cond = std(dist_to_conditioned_fish, 0, 2); % sd
se_dist_to_cond = std_dist_to_cond/sqrt(no_it); % se

mean_dist_to_cond_x = mean(dist_to_conditioned_fish_x,2);
std_dist_to_cond_x = std(dist_to_conditioned_fish_x, 0, 2);
se_dist_to_cond_x = std_dist_to_cond_x/sqrt(no_it);

plt_count = plt_count + 1;
fig = figure(plt_count); 
fig.Position = [300, 1200, 800, 700];

errorbar((2:n)-0.075, mean_dist_to_cond_x, se_dist_to_cond_x, ...
    'LineStyle', "none", "Marker", "s", ...
    "Color", ini_color, 'LineWidth', lw_plot, 'MarkerSize', ...
    10, 'MarkerFaceColor', ini_color)
hold on
errorbar((2:n)-0.075, mean_dist_to_cond, se_dist_to_cond, 'LineStyle', ...
    "none", "Marker", "d", ...
    "Color", ini_color, 'LineWidth', lw_plot, 'MarkerSize', 10, ...
    'MarkerFaceColor', ini_color)

legend({'x^i_{C}', 'r^i_{C}'}, 'Location', 'best')
legend('Box', 'off')

set(gca, 'XLim', [1.7 5+0.3], 'XTick', 1:n, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel("Distance to conditioned fish (cm)", 'FontSize', 23)


%% distribution of distance to min tank x

% plt_count = plt_count + 1;
% figure(plt_count)
% 
% mean_dist_min_wall_x = mean(dist_min_wall, 2); % mean distance to closest x-wall
% std_dist_min_wall_x = std(dist_min_wall, 0, 2);
% se_dist_min_wall_x = std_dist_min_wall_x/sqrt(no_it);
% 
% for i = 1:n
% 
%     scatter(i, dist_min_wall(i,:), 10, "k", 'filled')
%     hold on
% 
% end
% ylabel("Distance to tank wall")
% xlim([0.5 5.5])
% 
% plt_count = plt_count + 1;
% figure(plt_count)
% 
% errorbar(1:n, mean_dist_min_wall_x, se_dist_min_wall_x, 'LineStyle', "none", "Marker", "o", ...
%     "Color", "k")
% xlim([0.5 5.5])
% ylabel("Distance to tank wall (x)")
% 
% %% distribution of distance to tank centre
% 
% plt_count = plt_count + 1;
% figure(plt_count)
% 
% mean_dist_to_cent = mean(dist_cent_wall, 2); % distance to centre of closest x-wall
% std_dist_to_cent = std(dist_cent_wall, 0, 2);
% se_dist_to_cent = std_dist_to_cent/sqrt(no_it);
% 
% errorbar(1:n, mean_dist_to_cent, se_dist_to_cent, 'LineStyle', "none", "Marker", "o", ...
%     "Color", "k")
% xlim([0.5 5.5])
% ylabel("Distance to tank centre")
% 
% plt_count = plt_count + 1;
% figure(plt_count)
% 
% for i = 1:n
% 
%     scatter(i, dist_cent_wall(i,:), 10, "k", 'filled')
%     hold on
% 
% end
% xlim([0.5 5.5])
% ylabel("Distance to tank centre")


