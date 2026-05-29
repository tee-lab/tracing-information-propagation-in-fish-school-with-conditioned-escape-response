close all
clear
clc

%% load data

all_bd_data = load("all_bd_data.mat"); % load all pos, vel data
load("all_bd_data.mat", "dt", "n", "no_exp", "tank_max_x", "tank_min_x", "tank_middle", "del_step")
time_scale = 2; % if 1 start time of escape is tcross_centre, if 2 then it is time of light on
smooth_window = 5;
min_cc_spd = 0.3; % minimum correlation strength
min_cc_vel = 0.5;
max_delay = 75; % max delay
min_ts_length = 150; % minimum frames to calculate correlation
no_treat = 3; % ini, escape, relax
t_cut_min = -0.5;
t_cut_max = 3;
min_lag = 8; % minimum lag for correlation to be significant

tcross_data = struct();

ini_color = "#33b1ff";
esp_color = "#fa4d56";
relax_color = "#198038";
speed_plt_color = ["#6929c4", "#1192e8", "#005d5d", "#9f1853", "#d2a106"];
font_size = 25;
lw_plot = 2;
lw_axis = 2;
lw_xline = 4;

%% storing all data

pol = []; % polarisation
pol_x = []; % alingnment along x axis
pol_y = []; % alignment along y axis
grp_coh = []; % cohesion
grp_coh_x = []; % cohesion along x axis
grp_coh_y = []; % cohesion along y axis
grp_coh_wo_con = [];

tcd_all = []; % time crossing difference for centre of the barrier t2 - t1
tcross_time_all = []; % time to cross after the first fish crossed
tmind_all = []; % time crossing difference for 1st barrier
tmaxd_all = []; % time crossing difference for last barrier
rank_ord_all = []; % ids of fish in the order the cross first, second and last barrier
time_barrier_all = [];

speed_all = []; % speed 
speed_ini_all = []; % store all speeds in initial phase
speed_escape_all = []; % store all speed in escape phase.
speed_relax_all = []; % store all speeds in relax phase

pol_ini_all = []; % store all pol in ini phase
pol_esp_all = []; % store all pol in esp phase
pol_relax_all = []; % store all pol in relax phase

gc_ini_all = []; % store all gc in ini phase
gc_esp_all = []; % store all gc in esp phase
gc_relax_all = []; % store all gc in relax phase

% ids of leader fish in
spd_leader_ini = []; % ini
spd_leader_esp = []; % escape
spd_leader_rlx = []; % relax

% ids of follower fish in
spd_follower_ini = []; % ini
spd_follower_esp = []; % escape
spd_follower_rlx = []; % relax

% storing lag and cc over all
cc_all = [];
cc_lag_esp = [];
cc_lag_relax = [];

% storing distance travelled by fish during ini, esp and relax phase
dist_travelled_ini = nan(n,no_exp);
dist_travelled_esp = nan(n,no_exp);
dist_travelled_relax = nan(n,no_exp);

% storing median speeds by fish during ini, esp and relax phase
median_speed_ini = nan(n,no_exp);
median_speed_esp = nan(n,no_exp);
median_speed_relax = nan(n,no_exp);

%%

for e = 1:no_exp
    
    % fish data for this given experiment
    pos_fish = all_bd_data.(strcat('pos_t_ex_', num2str(e))); % input what pos data to use. 
    vel_fish = all_bd_data.(strcat('vel_t_ex_', num2str(e))); % input what vel data to use
    speed_fish = all_bd_data.(strcat('speed_t_ex_', num2str(e))); % input what speed data to use
    time_light_on = all_bd_data.(strcat('frame_light_on_ex_', num2str(e))); % when was the light turned on
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

    time_barrier = (max(fcross_tmax) - min(fcross_tmin))*dt;
    time_barrier_all = cat(1, time_barrier_all, time_barrier);
    
    % ranks of fish as they cross the barrier. 1 - 1st fish that cross the
    % barrier, 2 - 2nd and so on

    [~, rank_id] = sort(fcross_tmin, 'ascend');
    [~, rank_id_centre] = sort(fcross_tcentre, 'ascend');
    [~, rank_id_max] = sort(fcross_tmax, 'ascend');
    rank_ord = [rank_id, rank_id_centre, rank_id_max];
    rank_ord_all = cat(1,rank_ord_all, rank_ord);
    
    % store tcross data and crossing rank
    tcross_data.(['tcross_min', num2str(e)]) = tcross_min;
    tcross_data.(['tcross_max', num2str(e)]) = tcross_max;
    tcross_data.(['rank_id_centre', num2str(e)]) = rank_id_centre;
    tcross_data.(['rank_id_min', num2str(e)]) = rank_id;
    tcross_data.(['rank_id_max', num2str(e)]) = rank_id_max;
    % tcross_data.(['time_light_on', num2str(e)]) = time_light_on;
    
    % rescale time. if 1 then t = 0 implies time at which 1st crossed the
    % center of the tank. if 2 then t = 0 is when light was turned on
    if time_scale == 1

        t_plt = (t_plt - tcross_centre_min)/(tcross_centre_max - tcross_centre_min);

    elseif time_scale == 2

        t_plt = (t_plt - time_light_on)/(tcross_centre_max - time_light_on);

    end

    % crossing time difference t2 - t1 (1st barrier, 2nd barrier and center of the tank)
    tcross_min_diff = nan(n,1);
    tcross_min_diff(1) = fcross_tmin(rank_id(1)) - (time_light_on/dt);
    if tcross_min_diff(1) < 0
        x = 1;
    end
    tcross_min_diff(2:n) = fcross_tmin(rank_id(2:end)) - fcross_tmin(rank_id(1:end-1));
    
    tcross_max_diff = nan(n,1);
    tcross_max_diff(1) = fcross_tmax(rank_id(1)) - (time_light_on/dt);
    tcross_max_diff(2:n) = fcross_tmax(rank_id(2:end)) - fcross_tmax(rank_id(1:end-1));
    
    tcross_centre_diff = nan(n,1);
    tcross_centre_diff(1) = fcross_tcentre(rank_id(1)) - (time_light_on/dt);
    tcross_centre_diff(2:n) = fcross_tcentre(rank_id(2:end)) - fcross_tcentre(rank_id(1:end-1));

    % group polarisation
    vel_norm = vel_fish./(vecnorm(vel_fish,2,2)+eps);
    vn_x = squeeze(vel_norm(:,1,:)); % vx
    vn_y = squeeze(vel_norm(:,2,:)); % vy
    mx = mean(vn_x,1); % mx
    my = mean(vn_y,1); % my
    m = sqrt(mx.^2 + my.^2); % m

    m_ini = m(1:round((time_light_on/dt)));
    m_esp = m(round((time_light_on/dt))+1:max(fcross_tmax));
    m_relax = m(max(fcross_tmax)+1:end);

    pol_ini_all = cat(1, pol_ini_all, m_ini');
    pol_esp_all = cat(1, pol_esp_all, m_esp');
    pol_relax_all = cat(1, pol_relax_all, m_relax');

    % m = smoothdata(m, 'gaussian', smooth_window);
    % mx = smoothdata(mx, 'gaussian', smooth_window);

    p_temp = m(t_st:t_skip:t_end);
    px_temp = abs(mx(t_st:t_skip:t_end));
    py_temp = abs(my(t_st:t_skip:t_end));
    p_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', p_temp(t_plt >= t_cut_min & t_plt <= t_cut_max)'];
    pol_x_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', px_temp(t_plt >= t_cut_min & t_plt <= t_cut_max)'];
    pol_y_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', py_temp(t_plt >= t_cut_min & t_plt <= t_cut_max)'];
    pol = cat(1, pol, p_temp);
    pol_x = cat(1, pol_x, pol_x_temp);
    pol_y = cat(1, pol_y, pol_y_temp);

    % cohesion (dispersion)
    gc_x = mean(pos_fish(:,1,:), 1); % group centre x
    gc_y = mean(pos_fish(:,2,:), 1); % group centre y

    % position in group centre reference frame
    pos_gc(:,1,:) = pos_fish(:,1,:) - gc_x;
    pos_gc(:,2,:) = pos_fish(:,2,:) - gc_y;
    dist_gc = vecnorm(pos_gc,2,2);
    mean_dist_gc = squeeze(mean(dist_gc,1));

    gc_ini_model = mean_dist_gc(1:round((time_light_on/dt)));
    gc_esp_model = mean_dist_gc(round((time_light_on/dt))+1:max(fcross_tmax));
    gc_relax_model = mean_dist_gc(max(fcross_tmax)+1:end);

    gc_ini_all = cat(1, gc_ini_all, gc_ini_model);
    gc_esp_all = cat(1, gc_esp_all, gc_esp_model);
    gc_relax_all = cat(1, gc_relax_all, gc_relax_model);

    % dispersion only along x axis
    pos_gc_x = squeeze(pos_gc(:,1,:));
    pos_gc_x = abs(pos_gc_x);
    mean_gc_x = mean(pos_gc_x, 1);

    % dispersion along y axis
    pos_gc_y = squeeze(pos_gc(:,2,:));
    pos_gc_y = abs(pos_gc_y);
    mean_gc_y = mean(pos_gc_y, 1);

    % mean_dist_gc = smoothdata(mean_dist_gc, 'gaussian', smooth_window);
    % mean_gc_x = smoothdata(mean_gc_x, 'gaussian', smooth_window);
    % mean_gc_y = smoothdata(mean_gc_y, 'gaussian', smooth_window);
    
    % storing gc - 2d, x and y.
    gc_temp = mean_dist_gc(t_st:t_skip:t_end);
    gcx_temp = mean_gc_x(t_st:t_skip:t_end);
    gcy_temp = mean_gc_y(t_st:t_skip:t_end);

    gc_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', gc_temp(t_plt >= t_cut_min & t_plt <= t_cut_max)];
    gcx_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', gcx_temp(t_plt >= t_cut_min & t_plt <= t_cut_max)'];
    gcy_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', gcy_temp(t_plt >= t_cut_min & t_plt <= t_cut_max)'];

    grp_coh = cat(1, grp_coh, gc_temp);
    grp_coh_x = cat(1, grp_coh_x, gcx_temp);
    grp_coh_y = cat(1, grp_coh_y, gcy_temp);

    % distance traveled by conditioned fish and naive during initial phase
    pos_fish_rank_ord_ini = pos_fish(rank_id,:,1:round((time_light_on/dt)));
    diff_pos_ro_ini = pos_fish_rank_ord_ini(:,:,2:end) - pos_fish_rank_ord_ini(:,:,1:end-1);
    dist_diff_pos_ro_ini = squeeze(vecnorm(diff_pos_ro_ini, 2, 2));
    dist_diff_pos_ro_ini = sum(dist_diff_pos_ro_ini,2);
    dist_travelled_ini(:,e) = dist_diff_pos_ro_ini;

    % distance travelled by conditioned fish and naive during escape phase
    pos_fish_rank_ord_esp = pos_fish(rank_id,:,round((time_light_on/dt))+1:max(fcross_tmax));
    diff_pos_ro_esp = pos_fish_rank_ord_esp(:,:,2:end) - pos_fish_rank_ord_esp(:,:,1:end-1);
    dist_diff_pos_ro_esp = squeeze(vecnorm(diff_pos_ro_esp, 2, 2));
    dist_diff_pos_ro_esp = sum(dist_diff_pos_ro_esp,2);
    dist_travelled_esp(:,e) = dist_diff_pos_ro_esp;

    % distance travelled by conditioned fish and naive during relax phase
    pos_fish_rank_ord_relax = pos_fish(rank_id,:,max(fcross_tmax)+1:end);
    diff_pos_ro_relax = pos_fish_rank_ord_relax(:,:,2:end) - pos_fish_rank_ord_relax(:,:,1:end-1);
    dist_diff_pos_ro_relax = squeeze(vecnorm(diff_pos_ro_relax, 2, 2));
    dist_diff_pos_ro_relax = sum(dist_diff_pos_ro_relax,2);
    dist_travelled_relax(:,e) = dist_diff_pos_ro_relax;

    % gc without conditioned fish
    pos_fish_wo_con = pos_fish(rank_id,:,:);
    pos_fish_wo_con = pos_fish_wo_con(2:end,:,:);

    pf_wc_gc_x = mean(pos_fish_wo_con(:,1,:), 1);
    pf_wc_gc_y = mean(pos_fish_wo_con(:,2,:), 1);
    pf_gc_wc(:,1,:) = pos_fish_wo_con(:,1,:) - pf_wc_gc_x;
    pf_gc_wc(:,2,:) = pos_fish_wo_con(:,2,:) - pf_wc_gc_y;
    dist_pf_gc_wc = vecnorm(pf_gc_wc,2,2);
    mean_dist_pf_gc_wc = squeeze(mean(dist_pf_gc_wc,1));

    grp_coh_wo_con_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', mean_dist_pf_gc_wc(t_plt >= t_cut_min & t_plt <= t_cut_max)];
    grp_coh_wo_con = cat(1, grp_coh_wo_con, grp_coh_wo_con_temp);

    % tcross vs fish
    fcross_tcentre = fcross_tcentre*dt; % time of each fish crossing the centre
    tcross_min_diff = tcross_min_diff*dt;
    tcross_centre_diff = tcross_centre_diff*dt; % tcr_centre(2) - tcr_centre(1)
    tcross_max_diff = tcross_max_diff*dt;

    % storing tcross data
    tmind_temp = [(1:n)', tcross_min_diff];
    tcd_temp = [(1:n)', tcross_centre_diff]; % t2 - t1
    tmaxd_temp = [(1:n)', tcross_max_diff];
    tc_temp = [(1:n)', fcross_tmin(rank_id) - (time_light_on/dt)]; % crossing time
    tmind_all = cat(1, tmind_all, tmind_temp);
    tcd_all = cat(1, tcd_all, tcd_temp);
    tmaxd_all = cat(1, tmaxd_all, tmaxd_temp);
    tcross_time_all = cat(1, tcross_time_all, tc_temp);
    
    % plotting speed time series
    % figure(fig_spd_ts)
    % nexttile
    speed_rank = speed_fish(rank_id,:);
    speed_ini = speed_rank(:,1:round((time_light_on/dt)));
    speed_esp = speed_rank(:,round((time_light_on/dt))+1:max(fcross_tmax));
    speed_relax = speed_rank(:,max(fcross_tmax)+1:end);

    % median speeds of conditioned and naive fish in ini phase
    speed_ini_all = cat(2, speed_ini_all, speed_ini);
    speed_escape_all = cat(2, speed_escape_all, speed_esp);
    speed_relax_all = cat(2, speed_relax_all, speed_relax);
    median_speed_ini(:,e) = mean(speed_ini,2);
    median_speed_esp(:,e) = mean(speed_esp,2);
    median_speed_relax(:,e) = mean(speed_relax,2);
    
    spd_all_temp = speed_rank(:,t_plt >= t_cut_min & t_plt <= t_cut_max);
    spd_all_temp = [t_plt(t_plt >= t_cut_min & t_plt <= t_cut_max)', spd_all_temp'];
    speed_all = cat(1, speed_all, spd_all_temp);

end

save('tcross_data.mat', '-struct', "tcross_data")

%% group polarisation

plt_count = 0;

no_pol_edges = 21;
no_time_edges = 51;
[~, pol_sort_id] = sort(pol(:,1));
pol = pol(pol_sort_id, :);

t_edges = linspace(t_cut_min, t_cut_max, no_time_edges);
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

% along x axis

[~, px_sort_id] = sort(pol_x(:,1));
pol_x = pol_x(px_sort_id, :);

tx_edges = linspace(t_cut_min, t_cut_max, no_time_edges);
px_edges = linspace(0, 1, no_pol_edges);
[histcount_px, ~, ~, t_bins, px_bins] = histcounts2(pol_x(:,1), pol_x(:,2), tx_edges, px_edges);

tx_edges = tx_edges(1:end-1) + (tx_edges(2) - tx_edges(1))/2;
px_edges = px_edges(1:end-1) + (px_edges(2) - px_edges(1))/2;

mean_px = nan(1,size(histcount_px,1));
std_px = nan(1,size(histcount_px,1));
se_px = nan(1,size(histcount_px,1));
median_px = nan(1,size(histcount_px,1));
for i = 1:size(histcount_px,1)
    t_id = t_bins == i;
    px_temp = pol_x(t_id,2);
    mean_px(i) = mean(px_temp);
    std_px(i) = std(px_temp);
    se_px(i) = std(px_temp)/sqrt(length(px_temp));
    median_px(i) = median(px_temp);
end

% along y axis

[~, py_sort_id] = sort(pol_y(:,1));
pol_y = pol_y(py_sort_id, :);

t_edges = linspace(t_cut_min, t_cut_max, no_time_edges);
py_edges = linspace(0, 1, no_pol_edges);
[histcount_py, ~, ~, t_bins, py_bins] = histcounts2(pol_y(:,1), pol_y(:,2), t_edges, py_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
py_edges = py_edges(1:end-1) + (py_edges(2) - py_edges(1))/2;

mean_py = nan(1,size(histcount_py,1));
std_py = nan(1,size(histcount_py,1));
se_py = nan(1,size(histcount_py,1));
median_py = nan(1,size(histcount_py,1));
for i = 1:size(histcount_py,1)
    t_id = t_bins == i;
    py_temp = pol_y(t_id,2);
    mean_py(i) = mean(py_temp);
    std_py(i) = std(py_temp);
    se_py(i) = std(py_temp)/sqrt(length(py_temp));
    median_py(i) = median(py_temp);
end

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

errorbar(t_edges, mean_pol, se_pol, 'o-', 'Color', '#A52A2A', 'LineWidth', lw_plot, ...
        'MarkerFaceColor', '#A52A2A')
hold on
errorbar(t_edges, mean_px, se_px, 'o-', 'Color', ini_color, 'LineWidth', lw_plot, ...
        'MarkerFaceColor', ini_color)
hold on
errorbar(t_edges, mean_py, se_py, 'o-', 'Color', relax_color, 'LineWidth', lw_plot, ...
        'MarkerFaceColor', relax_color)
hold on

xline(0, '--r', 'LineWidth', lw_xline)
hold on
xline(1, '--r', 'LineWidth', lw_xline)
hold off

legend({'P', 'P_x', 'P_y'}, 'Location', 'northeast', 'FontSize', 32)
legend('Box', 'off')
xlabel('t_n')
ylabel('Polarisation')
set(gca, 'XLim', [-0.5 3], 'YLim', [0,1], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', 'FontSize', font_size, ...
    'FontName', 'Helvetica')

% exportgraphics(gca, 'pol_ts_data.pdf', 'ContentType', 'vector')

%% Group cohesion 2d.

[~, gc_sort_id] = sort(grp_coh(:,1));
grp_coh = grp_coh(gc_sort_id, :);

t_edges = linspace(t_cut_min, t_cut_max, 50);
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

% gc along x-axis

[~, gcx_sort_id] = sort(grp_coh_x(:,1));
grp_coh_x = grp_coh_x(gcx_sort_id, :);

t_edges = linspace(t_cut_min, t_cut_max, 50);
min_gcx = min(grp_coh_x(:,2));
max_gcx = max(grp_coh_x(:,2));
gcx_edges = linspace(min_gcx, max_gcx, 31);
[histcount_gcx, ~, ~, t_bins, gcx_bins] = histcounts2(grp_coh_x(:,1), grp_coh_x(:,2), t_edges, gcx_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
gcx_edges = gcx_edges(1:end-1) + (gcx_edges(2) - gcx_edges(1))/2;

mean_gcx = nan(1,size(histcount_gcx,1));
std_gcx = nan(1,size(histcount_gcx,1));
se_gcx = nan(1,size(histcount_gcx,1));
median_gcx = nan(1,size(histcount_gcx,1));
for i = 1:size(histcount_gcx,1)
    t_id = t_bins == i;
    gc_temp = grp_coh_x(t_id,2);
    mean_gcx(i) = mean(gc_temp);
    std_gcx(i) = std(gc_temp);
    se_gcx(i) = std(gc_temp)/sqrt(length(gc_temp));
    median_gcx(i) = median(gc_temp);
end

% Group cohesion along y-axis.

[~, gcy_sort_id] = sort(grp_coh_y(:,1));
grp_coh_y = grp_coh_y(gcy_sort_id, :);

t_edges = linspace(t_cut_min, t_cut_max, 50);
min_gcy = min(grp_coh_y(:,2));
max_gcy = max(grp_coh_y(:,2));
gcy_edges = linspace(min_gcy, max_gcy, 31);
[histcount_gcy, ~, ~, t_bins, gcy_bins] = histcounts2(grp_coh_y(:,1), grp_coh_y(:,2), t_edges, gcy_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
gcy_edges = gcy_edges(1:end-1) + (gcy_edges(2) - gcy_edges(1))/2;

mean_gcy = nan(1,size(histcount_gcy,1));
std_gcy = nan(1,size(histcount_gcy,1));
se_gcy = nan(1,size(histcount_gcy,1));
median_gcy = nan(1,size(histcount_gcy,1));
for i = 1:size(histcount_gcy,1)
    t_id = t_bins == i;
    gc_temp = grp_coh_y(t_id,2);
    mean_gcy(i) = mean(gc_temp);
    std_gcy(i) = std(gc_temp);
    se_gcy(i) = std(gc_temp)/sqrt(length(gc_temp));
    median_gcy(i) = median(gc_temp);
end

[~, gc_wo_con_sort_id] = sort(grp_coh_wo_con(:,1));
grp_coh_wo_con = grp_coh_wo_con(gc_wo_con_sort_id, :);

t_edges = linspace(t_cut_min, t_cut_max, 50);
min_gc_wo_con = min(grp_coh_wo_con(:,2));
max_gc_wo_con = max(grp_coh_wo_con(:,2));
gc_wc_edges = linspace(min_gc_wo_con, max_gc_wo_con, 31);
[histcount_gc, ~, ~, t_bins, gc_bins] = histcounts2(grp_coh_wo_con(:,1), grp_coh_wo_con(:,2), t_edges, gc_wc_edges);

t_edges = t_edges(1:end-1) + (t_edges(2) - t_edges(1))/2;
gc_wc_edges = gc_wc_edges(1:end-1) + (gc_wc_edges(2) - gc_wc_edges(1))/2;

mean_gc_wo_con = nan(1,size(histcount_gc,1));
std_gc_wo_con = nan(1,size(histcount_gc,1));
se_gc_wo_con = nan(1,size(histcount_gc,1));
median_gc_wo_con = nan(1,size(histcount_gc,1));
for i = 1:size(histcount_gc,1)
    t_id = t_bins == i;
    gc_temp = grp_coh_wo_con(t_id,2);
    mean_gc_wo_con(i) = mean(gc_temp);
    std_gc_wo_con(i) = std(gc_temp);
    se_gc_wo_con(i) = std(gc_temp)/sqrt(length(gc_temp));
    median_gc_wo_con(i) = median(gc_temp);
end

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

errorbar(t_edges, mean_gc, se_gc, 'o-', 'Color', '#A52A2A', 'LineWidth', lw_plot, ...
        'MarkerFaceColor', '#A52A2A')
hold on
errorbar(t_edges, mean_gcx, se_gcx, 'o-', 'Color', ini_color, 'LineWidth', lw_plot, ...
        'MarkerFaceColor', ini_color)
hold on
errorbar(t_edges, mean_gcy, se_gcy, 'o-', 'Color', relax_color, 'LineWidth', lw_plot, ...
        'MarkerFaceColor', relax_color)
hold on
errorbar(t_edges, mean_gc_wo_con, se_gc_wo_con, 'o-', ...
    'Color', 'k', 'LineWidth', lw_plot, 'MarkerFaceColor', 'k')
hold on

xline(0, '--r', 'LineWidth', lw_xline)
hold on
xline(1, '--r', 'LineWidth', lw_xline)
hold off

legend({'D', 'D_x', 'D_y', 'D_{wc}'}, 'Location', 'northeast', 'FontSize', 32)
legend('Box', 'off')

xlabel('t_n')
ylabel('Dispersion (cm)')

set(gca, 'XLim', [-0.5 3], 'YLim', [0,8], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

% exportgraphics(gca, 'disp_ts_data.pdf', 'ContentType', 'vector')

%% Time difference

mean_td_min = nan(1,n);
std_td_min = nan(1,n);
se_td_min = nan(1,n);
median_td_min = nan(1,n);

for i = 1:n
    td_id = tmind_all(:,1) == i;
    td_temp = tmind_all(td_id,2);
    mean_td_min(i) = mean(td_temp);
    std_td_min(i) = std(td_temp);
    se_td_min(i) = (std(td_temp))/sqrt(length(td_temp));
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
    se_td_centre(i) = (std(td_temp))/sqrt(length(td_temp));
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
    se_td_max(i) = (std(td_temp))/sqrt(length(td_temp));
    median_td_max(i) = median(td_temp);
end

% from model

tc_min_model = readmatrix('tc_min_mual_0.1_muat_2.3_muesp_1.25_K_4_k_2_sight_2.36_gamma_0.25_omegaini_0.7.csv');
st_id = 2;
mean_id = 2;
se_id = 4;
mean_td_min_model = tc_min_model(st_id:end,mean_id)*dt;
se_td_min_model = tc_min_model(st_id:end,se_id)*dt;

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

st_id = 2; % start id - should 1 (include cond fish) or 2 (exclude cond fish)
errorbar((st_id:n) - 0.05, mean_td_min(st_id:end), se_td_min(st_id:end), 'LineStyle', "none", ...
    "Marker", "o", "Color", ini_color, 'LineWidth', lw_plot, 'MarkerSize', ...
    10, 'MarkerFaceColor', ini_color)
hold on
errorbar((st_id:n) + 0.05, mean_td_min_model, se_td_min_model, 'o', ...
    'LineWidth', lw_plot, 'Color', esp_color, 'MarkerSize', 10, ...
    'MarkerFaceColor', esp_color)

% hold on
% errorbar((st_id:n), mean_td_centre(st_id:end), se_td_centre(st_id:end), 'LineStyle', "none", ...
%     "Marker", "o", "Color", esp_color, 'LineWidth', lw_plot, 'MarkerSize', ...
%     10, 'MarkerFaceColor', esp_color)
% hold on
% errorbar((st_id:n)+0.1, mean_td_max(st_id:end), se_td_max(st_id:end), ...
%     'LineStyle', "none", "Marker", "o", "Color", relax_color, 'LineWidth', ...
%     lw_plot, 'MarkerSize', 10, 'MarkerFaceColor', relax_color)

% legend({'Experiment', 'Model'}, 'Location', 'best')
% legend('Box', 'off')

set(gca, 'XLim', [st_id-0.2 n+0.2], 'XTick', st_id:n, ...
    'YLim', [0, 4], 'YTick', 0:4,...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel('Time since previous fish crossed (s)', 'FontSize', 23)

% set(gca, 'XLim', [st_id-0.2 n+0.2], 'XTick', st_id:n, 'XTickLabel', {"0", "1", "2", "3", "4"}, ...
%     'YLim', [0, 11], ...
%     'LineWidth', 2, 'Xcolor', 'k', 'YColor', 'k', 'FontWeight', 'bold', ...
%     'FontSize', font_size)

% exportgraphics(gca, 'time_sice_prev_fish_data.pdf', 'ContentType', 'vector')

%% time to cross

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

% model 

tc_model = readmatrix('mean_tc_model.csv');
mean_tc_min_model = tc_model(:,mean_id);
se_tc_min_model = tc_model(:,se_id);

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

errorbar((st_id:n) - 0.075, mean_tcross(st_id:n), se_tcross(st_id:n), 'LineStyle', "none", ...
    "Marker", "o", "Color", ini_color, 'LineWidth', lw_plot, 'MarkerSize', ...
    10, 'MarkerFaceColor', ini_color)
hold on
errorbar((st_id:n) + 0.075, mean_tc_min_model(st_id:n), se_tc_min_model, 'LineStyle', "none", ...
    "Marker", "o", "Color", esp_color, 'LineWidth', lw_plot, 'MarkerSize', ...
    10, 'MarkerFaceColor', esp_color)

legend({'Experiment', 'Model'}, 'Location', 'south')
legend('Box', 'off', 'Orientation', 'vertical', 'EdgeColor', [0.9, 0.9, 0.9], ...
    'LineWidth', 1)

set(gca, 'XLim', [st_id-0.2 n+0.2], 'XTick', st_id:n, ...
    'YLim', [0,18], 'YTick', 0:3:18, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel('Time since green light on (s)', 'FontSize', 23)

% exportgraphics(gca, 'time_sice_green_light_on_data.pdf', 'ContentType', 'vector')

%% plotting all speed

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

min_spd = min(min(speed_all(:,2:end)));
max_spd = max(max(speed_all(:,2:end)));

for i = 2:n+1

    [~, spd_sort_id] = sort(speed_all(:,1));
    speed_all = speed_all(spd_sort_id, :);

    t_edges = linspace(t_cut_min, t_cut_max, 51);
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

    errorbar(t_edges, mean_spd_i, se_spd_i, '-o', 'LineWidth', 1.5, ...
        'Color', speed_plt_color(i-1))
    hold on

end

xline(0, '--r', 'LineWidth', lw_xline)
hold on
xline(1, '--r', 'LineWidth', lw_xline)
hold off

legend({'CR-1', 'CR-2', 'CR-3', 'CR-4', 'CR-5'}, ...
    'Location', 'northeast', 'FontSize', 32)
legend('box','off')

set(gca, 'XLim', [-0.5 3], 'YLim', [0, 5], 'YTick', 0:0.5:5, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('t_n')
ylabel('Speed (cm/s)')

% exportgraphics(gca, 'speed_ts_data.pdf', 'ContentType', 'vector')

%% distance travelled by naive and conditioned fish in ini phase

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

dist_travelled_ini = dist_travelled_ini';
mean_dist_travel_ini = mean(dist_travelled_ini,1);
std_dist_travel_ini = std(dist_travelled_ini,1);
std_error_dt_ini = std_dist_travel_ini/sqrt(no_exp);

% distance travelled by naive and conditioned fish in relax phase

dist_travelled_relax = dist_travelled_relax';
mean_dist_travel_relax = mean(dist_travelled_relax,1);
std_dist_travel_relax = std(dist_travelled_relax,1);
std_error_dt_relax = std_dist_travel_relax/sqrt(no_exp);

errorbar(1:n, mean_dist_travel_ini, std_error_dt_ini, 'o', ...
    'LineWidth', lw_plot)
hold on
errorbar(1:n, mean_dist_travel_relax, std_error_dt_relax, 'o', ...
    'LineWidth', lw_plot)
set(gca, 'XLim', [.9 5.1], 'XTick', 1:5, 'YLim', [0,100], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

legend({'Initial', 'Relax'}, 'Location', 'best')
legend('Box', 'off')

xlabel('Crossing rank')
ylabel('Distance travelled (cm)')

% exportgraphics(gca, 'dist_trav_ini_relax_data.pdf', 'ContentType', 'vector')

%% mean speed by naive and conditioned fish in ini

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

mean_spd_ini = mean(speed_ini_all,2);
std_spd_ini = std(speed_ini_all,0, 2);
std_err_spd_ini = std_spd_ini/sqrt(no_exp);

% mean speed by naive and conditioned fish in relax phase

mean_spd_relax = mean(speed_relax_all,2);
std_spd_relax = std(speed_relax_all,0,2);
std_err_spd_relax = (std_spd_relax)/sqrt(no_exp);

errorbar((1:n)-0.1, mean_spd_ini, std_err_spd_ini, 'o', ...
    'LineWidth', lw_plot)
hold on
errorbar((1:n)+0.1, mean_spd_relax, std_err_spd_relax, 'o', ...
    'LineWidth', lw_plot);

set(gca, 'XLim', [.7 5.3], 'XTick', 1:5, 'YLim', [0,2], 'YTick', 0:0.5:2, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank')
ylabel('Mean speed (cm/s)')

legend({'Initial', 'Relax'}, 'Location', 'best')
legend('Box', 'off')

% exportgraphics(gca, 'speed_ini_relax_data.pdf', 'ContentType', 'vector')

%% plotting distributions of speeds in initial, eps and relax phase
pdf_lw = 3;

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

spd_ini_model = readmatrix("spd_ini_hist_model.csv");
spd_ini_edges_model = spd_ini_model(:,1);
spd_ini_hist_model = spd_ini_model(:,2);
spd_esp_model = readmatrix("spd_esp_hist_model.csv");
spd_esp_edges_model = spd_esp_model(:,1);
spd_esp_hist_model = spd_esp_model(:,2);
spd_relax_model = readmatrix('spd_relax_hist_model.csv');
spd_relax_edges_model = spd_relax_model(:,1);
spd_relax_hist_model = spd_relax_model(:,2);

speed_ini_all = speed_ini_all(:);
speed_escape_all = speed_escape_all(:);
speed_relax_all = speed_relax_all(:);

[spd_ini_hist, spd_ini_edges] = histcounts(speed_ini_all, 'Normalization', 'pdf');
[spd_esp_hist, spd_esp_edges] = histcounts(speed_escape_all, 'Normalization', 'pdf');
[spd_relax_hist, spd_relax_edges] = histcounts(speed_relax_all, 'Normalization', 'pdf');

plot(spd_ini_edges(1:end-1), spd_ini_hist, '-', 'Color', ini_color, 'LineWidth', pdf_lw)
hold on
plot(spd_esp_edges(1:end-1), spd_esp_hist, '-', 'Color', esp_color, 'LineWidth', pdf_lw)
hold on
plot(spd_relax_edges(1:end-1), spd_relax_hist, '-', 'Color', relax_color, 'LineWidth', pdf_lw)
hold on
plot(spd_ini_edges_model, spd_ini_hist_model, '--', 'Color', ini_color, 'LineWidth', pdf_lw)
hold on
plot(spd_esp_edges_model, spd_esp_hist_model, '--', 'Color', esp_color, 'LineWidth', pdf_lw)
hold on
plot(spd_relax_edges_model, spd_relax_hist_model, '--', 'Color', relax_color, 'LineWidth', pdf_lw)
    
set(gca, 'LineWidth', lw_axis, 'XLim', [0.03,8], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Speed (cm/s)')
ylabel('PDF')

legend({'Initial', 'Escape', 'Relax'}, 'Location', 'best')
legend('Box', 'off')

% exportgraphics(gca, 'speed_dist.pdf', 'ContentType', 'vector')

%% plotting polarisation in ini, esp and relax phase

p_ini_model = readmatrix('pol_ini_hist_model.csv');
p_ini_edges_model = p_ini_model(:,1);
p_ini_hist_model = p_ini_model(:,2);
p_esp_model = readmatrix('pol_esp_hist_model.csv');
p_esp_edges_model = p_esp_model(:,1);
p_esp_hist_model = p_esp_model(:,2);
p_relax_model = readmatrix('pol_relax_hist_model.csv');
p_relax_edges_model = p_relax_model(:,1);
p_relax_hist_model = p_relax_model(:,2);

p_edges = 41;

[p_ini_hist, p_ini_edges] = histcounts(pol_ini_all, p_edges, 'Normalization', 'pdf');
[p_esp_hist, p_esp_edges] = histcounts(pol_esp_all, p_edges, 'Normalization', 'pdf');
[p_relax_hist, p_relax_edges] = histcounts(pol_relax_all, p_edges, 'Normalization', 'pdf');

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(p_ini_edges(1:end-1), p_ini_hist, '-', 'Color', ini_color, 'LineWidth', pdf_lw)
hold on
plot(p_esp_edges(1:end-1), p_esp_hist, '-', 'Color', esp_color, 'LineWidth', pdf_lw)
hold on
plot(p_relax_edges(1:end-1), p_relax_hist, '-', 'Color', relax_color, 'LineWidth', pdf_lw)
hold on
plot(p_ini_edges_model, p_ini_hist_model, '--', 'Color', ini_color, 'LineWidth', pdf_lw)
hold on
plot(p_esp_edges_model, p_esp_hist_model, '--', 'Color', esp_color, 'LineWidth', pdf_lw)
hold on
plot(p_relax_edges_model, p_relax_hist_model, '--', 'Color', relax_color, 'LineWidth', pdf_lw)

set(gca, 'LineWidth', lw_axis, 'XLim', [0,1], ...
        'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
         'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Polarisation')
ylabel('PDF')

legend({'Initial', 'Escape', 'Relax'}, 'Location', 'best')
legend('Box', 'off')

% exportgraphics(gca, 'pol_dist.pdf', 'ContentType', 'vector')

%% plotting distribution of group cohesion

gc_ini_model = readmatrix('gc_ini_hist_model.csv');
gc_ini_edges_model = gc_ini_model(:,1);
gc_ini_hist_model = gc_ini_model(:,2);
gc_esp_model = readmatrix('gc_esp_hist_model.csv');
gc_esp_edges_model = gc_esp_model(:,1);
gc_esp_hist_model = gc_esp_model(:,2);
gc_relax_model = readmatrix('gc_relax_hist_model.csv');
gc_relax_edges_model = gc_relax_model(:,1);
gc_relax_hist_model = gc_relax_model(:,2);

[gc_ini_hist, gc_ini_edges] = histcounts(gc_ini_all, 'Normalization', 'pdf');
[gc_esp_hist, gc_esp_edges] = histcounts(gc_esp_all, 'Normalization', 'pdf');
[gc_relax_hist, gc_relax_edges] = histcounts(gc_relax_all, 'Normalization', 'pdf');

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

plot(gc_ini_edges(1:end-1), gc_ini_hist, '-', 'Color', ini_color, ...
    'LineWidth', pdf_lw)
hold on
plot(gc_esp_edges(1:end-1), gc_esp_hist, '-', 'Color', esp_color, ...
    'LineWidth', pdf_lw)
hold on
plot(gc_relax_edges(1:end-1), gc_relax_hist, '-', 'Color', relax_color, ...
    'LineWidth', pdf_lw)
hold on
plot(gc_ini_edges_model, gc_ini_hist_model, '--', 'Color', ini_color, ...
    'LineWidth', pdf_lw)
hold on
plot(gc_esp_edges_model, gc_esp_hist_model, '--', 'Color', esp_color, ...
    'LineWidth', pdf_lw)
hold on
plot(gc_relax_edges_model, gc_relax_hist_model, '--', 'Color', relax_color, ...
    'LineWidth', pdf_lw)
    
set(gca, 'LineWidth', 2, 'XLim', [0,12], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Dispersion (cm)')
ylabel('PDF')

legend({'Initial', 'Escape', 'Relax'}, 'Location', 'best')
legend('Box', 'off')

% exportgraphics(gca, 'disp_dist.pdf', 'ContentType', 'vector')

%% storing time since last fish crossed data as .csv

% video = 1:no_exp;
% naive_agents = (1:n)';
% naive_agents = repmat(naive_agents, no_exp, 1);
% video = repmat(video, n, 1);
% video = video(:);
% % naive_fish_id = tmind_all(:,1) ~= 1;
% % time_last_fish_crossed_min = tmind_all(naive_fish_id,2);
% time_last_fish_crossed_min = tmind_all(:,2);
% 
% time_since_last_fish_crossed_data = [video, naive_agents, time_last_fish_crossed_min];
% writematrix(time_since_last_fish_crossed_data, 'time_since_last_fish_crossed_data.csv')