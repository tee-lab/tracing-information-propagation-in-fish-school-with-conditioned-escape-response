close all
clear
clc

tic

%% fixed parameters

no_pol_edges = 21;
no_gc_edges = 31;
no_time_edges = 501;
p_edges = 41;

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

%% loading all .mat files and calculating group properties

fileInfo = dir('*.mat');
fname = fileInfo.name;

load(fname) % load data
fname = fname(1:end-4);

no_it = no_it*no_exp; % total no.of iterations (to calculate group properties)
en_start = reshape(en_start, n, no_it); % start time
en_centre = reshape(en_centre, n, no_it); % centre crossing time
en_end = reshape(en_end, n, no_it); % crossing time of 2nd border.
rank_order_atk = reshape(rank_order_atk, n, no_it); % rank of crossing
pos_t = reshape(pos_t, n, 2, n_iter, no_it); % position data
s_t = reshape(s_t, n, n_iter, no_it); % speed data
theta_t = reshape(theta_t, n, n_iter, no_it); % orientation data

% identify when all fish crossed the barrier
all_not_cross_id = isnan(en_end);
all_not_cross_id = sum(all_not_cross_id,1);
frac_not_cross = all_not_cross_id/n;
frac_not_cross = mean(frac_not_cross);
all_not_cross_id = all_not_cross_id == 0;

% choose only those simulations where all fish crossed barrier
en_start = en_start(:,all_not_cross_id);
en_centre = en_centre(:,all_not_cross_id);
en_end = en_end(:,all_not_cross_id);
rank_order_atk = rank_order_atk(:,all_not_cross_id);
pos_t = pos_t(:,:,:,all_not_cross_id);
s_t = s_t(:,:,all_not_cross_id);
theta_t = theta_t(:,:,all_not_cross_id);
no_it = sum(all_not_cross_id);

% Encounter time

% Escape phase time

esp_time_st = nan(no_it,1); % time at which fish crosses the first barrier
esp_time_end = nan(no_it,1); % time at which fish crosses the second barrier

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

% different phases
t_ini = find(t_plt < 0);
t_esp = t_plt > 0 & t_plt < 1;
t_relax = find(t_plt > 1);

% Grp polarisation

grp_pol = nan(n_iter, no_it);

for i = 1:no_it

    theta_temp = squeeze(theta_t(:,:,i));
    vel_x = mean(cos(theta_temp),1);
    vel_y = mean(sin(theta_temp),1);
    pol_tem = sqrt(vel_x.^2 + vel_y.^2);
    grp_pol(:,i) = pol_tem.';

end

pol = [t_plt, grp_pol(:)];

grp_pol_ini = pol(t_ini,2);
grp_pol_esp = pol(t_esp,2);
grp_pol_relax = pol(t_relax,2);

[p_ini_hist, p_ini_edges] = histcounts(grp_pol_ini, p_edges, 'Normalization', 'pdf');
[p_esp_hist, p_esp_edges] = histcounts(grp_pol_esp, p_edges, 'Normalization', 'pdf');
[p_relax_hist, p_relax_edges] = histcounts(grp_pol_relax, p_edges, 'Normalization', 'pdf');

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

% pol in ini, esp, relax

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(p_ini_edges(1:end-1), p_ini_hist, '-', 'Color', ini_color, 'LineWidth', lw_xline)
hold on
plot(p_esp_edges(1:end-1), p_esp_hist, '-', 'Color', esp_color, 'LineWidth', lw_xline)
hold on
plot(p_relax_edges(1:end-1), p_relax_hist, '-', 'Color', relax_color, 'LineWidth', lw_xline)

set(gca, 'LineWidth', lw_axis, 'XLim', [0,1], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

legend({'Initial', 'Escape', 'Relax'}, 'Location', 'best')
legend('Box', 'off')

xlabel('Polarisation')
ylabel('PDF')

% Group speed dynamics

median_grp_spd = nan(n_iter, no_it);
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

% storing pdfs of speed in ini, esp, and relax phase

speed_ini_all = speed_rank(t_ini,:);
speed_ini_all = speed_ini_all(:);
id_min_spd_ini = speed_ini_all > min_spd;
speed_ini_all = speed_ini_all(id_min_spd_ini);

speed_escape_all = speed_rank(t_esp,:);
speed_escape_all = speed_escape_all(:);
id_min_spd_esp = speed_escape_all > min_spd;
speed_escape_all = speed_escape_all(id_min_spd_esp);

speed_relax_all = speed_rank(t_relax,:);
speed_relax_all = speed_relax_all(:);
id_min_spd_relax = speed_relax_all > min_spd;
speed_relax_all = speed_relax_all(id_min_spd_relax);

[spd_ini_hist, spd_ini_edges] = histcounts(speed_ini_all, 'Normalization', 'pdf');
[spd_esp_hist, spd_esp_edges] = histcounts(speed_escape_all, 'Normalization', 'pdf');
[spd_relax_hist, spd_relax_edges] = histcounts(speed_relax_all, 'Normalization', 'pdf');

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(spd_ini_edges(1:end-1), spd_ini_hist, '-', 'Color', ini_color, 'LineWidth', lw_xline)
hold on
plot(spd_esp_edges(1:end-1), spd_esp_hist, '-', 'Color', esp_color, 'LineWidth', lw_xline)
hold on
plot(spd_relax_edges(1:end-1), spd_relax_hist, '-', 'Color', relax_color, 'LineWidth', lw_xline)

set(gca, 'LineWidth', lw_axis, 'XLim', [0.03,8], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Speed (cm/s)')
ylabel('PDF')

legend({'Initial', 'Escape', 'Relax'}, 'Location', 'best')
legend('Box', 'off')

% Group cohesion

grp_cohesion = nan(n_iter, no_it); % grp cohesion 2D
grp_coh_x = nan(n_iter, no_it); % gc along x axis
grp_coh_y = nan(n_iter, no_it); % gc along y axis
grp_coh_wcf = nan(n_iter, no_it); % gc without conditioned fish

for i = 1:no_it

    pos_it = pos_t(rank_order_atk(:,i),:,:,i);

    gc_x = mean(pos_it(:,1,:), 1);
    gc_y = mean(pos_it(:,2,:), 1);

    pos_gc(:,1,:) = pos_it(:,1,:) - gc_x;
    pos_gc(:,2,:) = pos_it(:,2,:) - gc_y;

    dist_gc = vecnorm(pos_gc,2,2);
    dist_gc = squeeze(dist_gc);
    dist_gc_wcf = dist_gc(2:end,:);
    mean_dist_gc_wcf = mean(dist_gc_wcf,1);
    mean_dist_gc = mean(dist_gc,1);

    grp_cohesion(:,i) = mean_dist_gc;
    grp_coh_wcf(:,i) = mean_dist_gc_wcf;

    % dispersion only along x axis
    pos_gc_x = squeeze(pos_gc(:,1,:));
    pos_gc_x = abs(pos_gc_x);
    mean_gc_x = mean(pos_gc_x, 1);
    grp_coh_x(:,i) = mean_gc_x;

    % dispersion along y axis
    pos_gc_y = squeeze(pos_gc(:,2,:));
    pos_gc_y = abs(pos_gc_y);
    mean_gc_y = mean(pos_gc_y, 1);
    grp_coh_y(:,i) = mean_gc_y;

end

grp_coh = [t_plt, grp_cohesion(:)];
grp_coh_x = [t_plt, grp_coh_x(:)];
grp_coh_y = [t_plt, grp_coh_y(:)];
grp_coh_wcf = [t_plt, grp_coh_wcf(:)];

% storing gc data
gc_ini = grp_coh(t_ini,2);
gc_esp = grp_coh(t_esp,2);
gc_relax = grp_coh(t_relax,2);

[gc_ini_hist, gc_ini_edges] = histcounts(gc_ini, 'Normalization', 'pdf');
[gc_esp_hist, gc_esp_edges] = histcounts(gc_esp, 'Normalization', 'pdf');
[gc_relax_hist, gc_relax_edges] = histcounts(gc_relax, 'Normalization', 'pdf');

% gc 2d

[~, gc_sort_id] = sort(grp_coh(:,1));
grp_coh = grp_coh(gc_sort_id, :);

t_edges = linspace(min(t_plt), max(t_plt), no_time_edges);
min_gc = min(grp_coh(:,2));
max_gc = max(grp_coh(:,2));
gc_edges = linspace(min_gc, max_gc, no_gc_edges);
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

% group cohesion along x axis

[~, gc_sort_id] = sort(grp_coh_x(:,1));
grp_coh_x = grp_coh_x(gc_sort_id, :);

t_edges = linspace(min(t_plt), max(t_plt), no_time_edges);
min_gc_x = min(grp_coh_x(:,2));
max_gc_x = max(grp_coh_x(:,2));
gc_x_edges = linspace(min_gc_x, max_gc_x, no_gc_edges);
[histcount_gc_x, ~, ~, t_bins, gc_bins] = histcounts2(grp_coh_x(:,1), grp_coh_x(:,2), t_edges, gc_x_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
gc_x_edges = gc_x_edges(1:end-1) + (gc_x_edges(2) - gc_x_edges(1))/2;

mean_gc_x = nan(1,size(histcount_gc_x,1));
std_gc_x = nan(1,size(histcount_gc_x,1));
se_gc_x = nan(1,size(histcount_gc_x,1));
median_gc_x = nan(1,size(histcount_gc_x,1));

for i = 1:size(histcount_gc_x,1)
    t_id = t_bins == i;
    gc_x_temp = grp_coh_x(t_id,2);
    mean_gc_x(i) = mean(gc_x_temp);
    std_gc_x(i) = std(gc_x_temp);
    se_gc_x(i) = std(gc_x_temp)/sqrt(length(gc_x_temp));
    median_gc_x(i) = median(gc_x_temp);
end

% gc along y axis

[~, gc_sort_id] = sort(grp_coh_y(:,1));
grp_coh_y = grp_coh_y(gc_sort_id, :);

t_edges = linspace(min(t_plt), max(t_plt), no_time_edges);
min_gc_y = min(grp_coh_y(:,2));
max_gc_y = max(grp_coh_y(:,2));
gc_y_edges = linspace(min_gc_y, max_gc_y, no_gc_edges);
[histcount_gc_y, ~, ~, t_bins, gc_bins] = histcounts2(grp_coh_y(:,1), grp_coh_y(:,2), t_edges, gc_y_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
gc_y_edges = gc_y_edges(1:end-1) + (gc_y_edges(2) - gc_y_edges(1))/2;

mean_gc_y = nan(1,size(histcount_gc_y,1));
std_gc_y = nan(1,size(histcount_gc_y,1));
se_gc_y = nan(1,size(histcount_gc_y,1));
median_gc_y = nan(1,size(histcount_gc_y,1));
for i = 1:size(histcount_gc_y,1)
    t_id = t_bins == i;
    gc_y_temp = grp_coh_y(t_id,2);
    mean_gc_y(i) = mean(gc_y_temp);
    std_gc_y(i) = std(gc_y_temp);
    se_gc_y(i) = std(gc_y_temp)/sqrt(length(gc_y_temp));
    median_gc_y(i) = median(gc_y_temp);
end

% gc without conditioned fish

[~, gc_sort_id] = sort(grp_coh_wcf(:,1));
grp_coh_wcf = grp_coh_wcf(gc_sort_id, :);

t_edges = linspace(min(t_plt), max(t_plt), no_time_edges);
min_gc_wcf = min(grp_coh_wcf(:,2));
max_gc_wcf = max(grp_coh_wcf(:,2));
gc_wcf_edges = linspace(min_gc_wcf, max_gc_wcf, no_gc_edges);
[histcount_gc_wcf, ~, ~, t_bins, gc_bins] = histcounts2(grp_coh_wcf(:,1), grp_coh_wcf(:,2), t_edges, gc_wcf_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
gc_wcf_edges = gc_wcf_edges(1:end-1) + (gc_wcf_edges(2) - gc_wcf_edges(1))/2;

mean_gc_wcf = nan(1,size(histcount_gc_wcf,1));
std_gc_wcf = nan(1,size(histcount_gc_wcf,1));
se_gc_wcf = nan(1,size(histcount_gc_wcf,1));
median_gc_wcf = nan(1,size(histcount_gc_wcf,1));

for i = 1:size(histcount_gc_wcf,1)
    t_id = t_bins == i;
    gc_wcf_temp = grp_coh_wcf(t_id,2);
    mean_gc_wcf(i) = mean(gc_wcf_temp);
    std_gc_wcf(i) = std(gc_wcf_temp);
    se_gc_wcf(i) = std(gc_wcf_temp)/sqrt(length(gc_wcf_temp));
    median_gc_wcf(i) = median(gc_wcf_temp);
end

% ploting gc_x, gc_y, gc, and gc without conditioned fish together

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(t_edges, mean_gc, 'o-', 'Color', '#A52A2A', 'LineWidth', lw_plot, ...
    'MarkerFaceColor', '#A52A2A')
hold on
plot(t_edges, mean_gc_x, 'o-', 'Color', ini_color, 'LineWidth', lw_plot, ...
    'MarkerFaceColor', ini_color)
hold on
plot(t_edges, mean_gc_y, 'o-', 'Color', relax_color, 'LineWidth', lw_plot, ...
    'MarkerFaceColor', relax_color)
hold on
plot(t_edges, mean_gc_wcf, 'o-', 'Color', 'k', 'LineWidth', lw_plot, 'MarkerFaceColor', 'k')
hold on
xline(0, '--r', 'LineWidth', lw_plot)
hold on
xline(1, '--r', 'LineWidth', lw_plot)
hold off

legend({'C', 'C_x', 'C_y', 'C_{wc}'}, 'Location', 'best')
legend('Box', 'off')

set(gca, 'XLim', [-0.5 3], 'YLim', [0, 7], 'LineWidth', lw_axis, ...
    'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('t_n', 'FontSize', label_fs)
ylabel('Dispersion (cm)', 'FontSize', label_fs)

% time crossing

[fcross_tmin, rank_id_min] = sort(en_start, 1, 'ascend'); % arranging start in ascending order.
[fcross_tcentre, rank_id_centre] = sort(en_centre, 1, 'ascend');
[fcross_tmax, rank_id_max] = sort(en_end, 1, 'ascend');

c_id = repmat((1:n)', no_it, 1);

tcross_min_diff = nan(n,no_it);
tcross_min_diff(1,:) = fcross_tmin(1,:) - t_atk;
tcross_min_diff(2:end,:) = fcross_tmin(2:end,:) - fcross_tmin(1:(end-1),:);
tmind_all = [c_id, tcross_min_diff(:)];

tcross_time_all = [c_id, (fcross_tmin(:) - t_atk)];

tcross_centre_diff = nan(n,no_it);
tcross_centre_diff(1,:) = fcross_tcentre(1,:) - t_atk;
tcross_centre_diff(2:end,:) = fcross_tcentre(2:end,:) - fcross_tcentre(1:(end-1),:);
tcd_all = [c_id, tcross_centre_diff(:)];

tcross_max_diff = nan(n,no_it);
tcross_max_diff(1,:) = fcross_tmax(1,:) - t_atk;
tcross_max_diff(2:end,:) = fcross_tmax(2:end,:) - fcross_tmax(1:(end-1),:);
tmaxd_all = [c_id, tcross_max_diff(:)];

%  Time difference plots

mean_td_min = nan(1,n);
std_td_min = nan(1,n);
se_td_min = nan(1,n);
median_td_min = nan(1,n);
en_time_diff = nan(n, no_it);

for i = 1:n
    td_id = tmind_all(:,1) == i;
    td_temp = tmind_all(td_id,2);
    en_time_diff(i,:) = td_temp;
    mean_td_min(i) = mean(td_temp);
    std_td_min(i) = std(td_temp);
    se_td_min(i) = (1.96*std(td_temp))/sqrt(length(td_temp));
    median_td_min(i) = median(td_temp);
end

mean_td_centre = nan(1,n);
std_td_centre = nan(1,n);
se_td_centre = nan(1,n);
median_td_centre = nan(1,n);

for i = 1:n
    td_id = tcd_all(:,1) == i;
    td_temp = tcd_all(td_id,2);
    mean_td_centre(i) = mean(td_temp);
    std_td_centre(i) = std(td_temp);
    se_td_centre(i) = (1.96*std(td_temp))/sqrt(length(td_temp));
    median_td_centre(i) = median(td_temp);
end

mean_td_max = nan(1,n);
std_td_max = nan(1,n);
se_td_max = nan(1,n);
median_td_max = nan(1,n);

for i = 1:n
    td_id = tmaxd_all(:,1) == i;
    td_temp = tmaxd_all(td_id,2);
    mean_td_max(i) = mean(td_temp);
    std_td_max(i) = std(td_temp);
    se_td_max(i) = (1.96*std(td_temp))/sqrt(length(td_temp));
    median_td_max(i) = median(td_temp);
end

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

st_id = 2;

errorbar((1:n)-0.1, mean_td_min*dt, se_td_min*dt, 'o', ...
        'LineWidth', lw_plot, 'Color', ini_color, 'MarkerSize', mr_size, ...
        'MarkerFaceColor', ini_color)
hold on
errorbar((1:n), mean_td_centre*dt, se_td_centre*dt, 'o', ...
        'LineWidth', lw_plot, 'Color', esp_color, 'MarkerSize', mr_size, ...
        'MarkerFaceColor', esp_color)
hold on
errorbar((1:n)+0.1, mean_td_max*dt, se_td_max*dt, 'o', ...
        'LineWidth', lw_plot, 'Color', relax_color, 'MarkerSize', mr_size, ...
        'MarkerFaceColor', relax_color)

legend({'barrier-start', 'barrier-centre', 'barrier-last'}, 'Location', 'best')
legend('Box', 'off')

set(gca, 'XLim', [st_id-0.4 n+0.2], 'XTick', 1:n, 'LineWidth', lw_axis, ...
        'Xcolor', 'k', 'YColor', 'k', ...
        'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', label_fs)
ylabel('Time since previous fish crossed (s)', 'FontSize', label_fs)

% time to cross

mean_tcross = nan(1,n);
std_tcross = nan(1,n);
se_tcross = nan(1,n);
median_tcross = nan(1,n);

st_id = 1;
for i = st_id:n
    tc_id = tcross_time_all(:,1) == i;
    td_temp = tcross_time_all(tc_id,2);
    td_temp = td_temp*dt;
    mean_tcross(i) = mean(td_temp);
    std_tcross(i) = std(td_temp);
    se_tcross(i) = (std(td_temp))/sqrt(length(td_temp));
    median_tcross(i) = median(td_temp);
end

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

% mean_tc_min_model = [(1:n)', mean_tcross', std_tcross', se_tcross', median_tcross'];
% writematrix(mean_tc_min_model, 'mean_tc_model.csv')

errorbar((st_id:n) - 0.075, mean_tcross(st_id:n), se_tcross(st_id:n), ...
    'LineStyle', "none", "Marker", "o", "Color", ini_color, ...
    'LineWidth', lw_plot, 'MarkerSize', 10, 'MarkerFaceColor', ini_color)

set(gca, 'XLim', [st_id-0.2 n+0.2], 'XTick', st_id:n, ...
    'YLim', [0,18], 'YTick', 0:3:18, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel('Time since green light on (s)', 'FontSize', 23)

% distance to conditioned fish

dist_min_wall = nan(n, no_it); % distance to the closest wall (only along x axis)
dist_cent_wall = nan(n, no_it); % distance to centre of closest wall (2-D)
dist_to_conditioned_fish_x = nan(n-1, no_it); % distance to conditioned fish (x-asis only)
dist_to_conditioned_fish = nan(n-1, no_it); % distance to conditioned fish
order_to_conditioned_fish = nan(n-1, no_it); % check the order if the one to cross the border is the closest to the conditioned fish
order_to_conditioned_fish_x = nan(n-1, no_it);

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

    [~, rank_pos_id] = sort(vecnorm(dist_to_cond,2,2), 'ascend');
    [~, rank_pos_id_x] = sort(abs(dist_to_cond(:,1)), 'ascend');
    order_to_conditioned_fish_x(:,e) = rank_pos_id_x;
    order_to_conditioned_fish(:,e) = rank_pos_id;

    % check which side of the tank fish are in

    dist_x  = mini_box_len - pos_fish_light_on(:,1); % x_min_wall - x
    dist_y = (box_width/2) - pos_fish_light_on(:,2); % distance to centre of the tank

    dist_min_wall(:,e) = dist_x;
    dist_cent_wall(:,e) = sqrt(dist_x.^2 + dist_y.^2);

end

% distribution of distance to conditioned fish

mean_dist_to_cond = mean(dist_to_conditioned_fish, 2); % mean dist to conditioned fish
std_dist_to_cond = std(dist_to_conditioned_fish, 0, 2); % sd
se_dist_to_cond = std_dist_to_cond/sqrt(no_it); % se

mean_dist_to_cond_x = mean(dist_to_conditioned_fish_x,2);
std_dist_to_cond_x = std(dist_to_conditioned_fish_x, 0, 2);
se_dist_to_cond_x = std_dist_to_cond_x/sqrt(no_it);

plt_count = plt_count + 1;
fig = figure(plt_count); fig.Position = [300, 1200, 800, 700];
figure(plt_count)

errorbar((2:n)-0.075, mean_dist_to_cond_x, se_dist_to_cond_x, ...
    'LineStyle', "none", "Marker", "s", ...
    "Color", ini_color, 'LineWidth', lw_plot, 'MarkerSize', ...
    10, 'MarkerFaceColor', ini_color)
hold on
errorbar((2:n)-0.075, mean_dist_to_cond, se_dist_to_cond, 'LineStyle', ...
    "none", "Marker", "d", ...
    "Color", ini_color, 'LineWidth', lw_plot, 'MarkerSize', 10, ...
    'MarkerFaceColor', ini_color)

legend({'x^i_{C}', 'r^i_{C}'}, 'Location', 'southwest')
legend('Box', 'off')

set(gca, 'XLim', [1.7 5+0.3], 'XTick', 1:n, 'YLim', [0, 7], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel("Distance to conditioned fish (cm)", 'FontSize', 23)

% relative orientation and viewing angle

psi_ic = nan(n-1, no_it); % viewing angle when the light was turned on
phi_ic = nan(n-1, no_it); % relative orientation when light was turned on

for e = 1:no_it

    % fish data for this given experiment
    pos_it = pos_t(:,:,:,e); % input what pos data to use.
    phi_it = theta_t(:,:,e); % orientation
    vel_it = nan(size(phi_it,1),2,size(phi_it,2));
    vel_it(:,1,:) = cos(phi_it);
    vel_it(:,2,:) = sin(phi_it);
    time_light_on = min(t_atk); % when was the light turned on
    rank_id_min = rank_order_atk(:,e); % rank of crossing

    % position when light was turned on in the order of crossing
    pos_fish_light_on = pos_it(rank_id_min,:,time_light_on);
    vel_fish_light_on = vel_it(rank_id_min,:,time_light_on);

    vel_naive = vel_fish_light_on(2:end,:);
    vel_naive = vel_naive./(vecnorm(vel_naive,2,2) + eps);
    % r_ic
    r_ic = pos_fish_light_on(1,:) - pos_fish_light_on(2:end,:);
    r_ic = r_ic./(vecnorm(r_ic,2,2) + eps);

    % calculating viewing angle
    psi_ic_temp = dot(vel_naive, r_ic, 2);
    psi_ic(:,e) = acos(psi_ic_temp);

    % vel conditioned fish
    vel_con = vel_fish_light_on(1,:);
    vel_con = vel_con./vecnorm(vel_con,2,2);
    vel_con = repmat(vel_con, n-1, 1);

    % calculating relative orientation
    phi_ic_temp = dot(vel_naive, vel_con, 2);
    phi_ic(:,e) = acos(phi_ic_temp);

end

% psi

psi_ic = psi_ic';
psi_ic = psi_ic*180/pi;

mean_psi_ic = mean(psi_ic,1);
std_psi_ic = std(psi_ic,1);
error_psi_ic = (std_psi_ic)/sqrt(no_it);

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

errorbar((2:n) - 0.05, mean_psi_ic, error_psi_ic, "o", 'MarkerEdgeColor',  ini_color, ...
    'MarkerFaceColor',  ini_color, 'color', ini_color, 'LineWidth', lw_plot, ...
    'MarkerSize', 10)

set(gca, 'XLim', [1.8 n+0.2], 'YLim', [0,120], 'XTick', 1:n, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel('Viewing angle (in degree)', 'FontSize', 23)

% phi

phi_ic = phi_ic';
phi_ic = phi_ic*180/pi;

mean_phi_ic = mean(phi_ic,1);
std_phi_ic = std(phi_ic,1);
error_phi_ic = (std_phi_ic)/sqrt(no_it);

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

errorbar((2:n)-0.05, mean_phi_ic, error_phi_ic, "o", 'MarkerEdgeColor',  ini_color, ...
    'MarkerFaceColor',  ini_color, 'Color', ini_color, 'LineWidth', lw_plot, ...
    'MarkerSize', 10)

set(gca, 'XLim', [1.8 n+0.2], 'YLim', [0,120], 'XTick', 1:n, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel('Relative orientation (in degree)', 'FontSize', 23)

toc