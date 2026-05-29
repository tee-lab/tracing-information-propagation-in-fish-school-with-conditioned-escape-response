close all
clear
clc

%% load data

all_bd_data = load('../main_text/all_bd_data.mat'); % load all pos, vel data
tcross_data = load('../main_text/tcross_data.mat');
load("../main_text/all_bd_data.mat", "dt", "n", "no_exp", "tank_max_x", "tank_min_x", "tank_middle")
smooth_window = 5;
tank_width = 20;

agent_col = ["red", "#00059f", "#2c2cff", "#4e91fd", "#bac2ff"];
t_cut_min = -0.5;
t_cut_max = 3;
min_lag = 0; % minimum lag for correlation to be significant.

font_size = 25;
lw_plot = 2;
lw_axis = 2;

ini_color = "#33b1ff";
esp_color = "#fa4d56";
relax_color = "#198038";

%%

dist_min_wall = nan(n,no_exp); % distance to the closest wall (only along x axis)
dist_cent_wall = nan(n,no_exp); % distance to centre of closest wall (2-D)
dist_to_conditioned_fish_x = nan(n-1, no_exp); % distance to conditioned fish (x-asis only)
dist_to_conditioned_fish = nan(n-1, no_exp); % distance to conditioned fish
order_to_conditioned_fish = nan(n-1, no_exp); % check the order if the one to cross the border is the closest to the conditioned fish
order_to_conditioned_fish_x = nan(n-1, no_exp);

for e = 1:no_exp
    
    % fish data for this given experiment
    pos_fish = all_bd_data.(strcat('pos_t_ex_', num2str(e))); % input what pos data to use.
    time_light_on = all_bd_data.(strcat('frame_light_on_ex_', num2str(e))); % when was the light turned on
    rank_id_min = tcross_data.(strcat('rank_id_min', num2str(e))); % rank of crossing
    
    % position when light was turned on in the order of crossing
    pos_fish_light_on = pos_fish(rank_id_min,:,time_light_on);
    % distance to conditioned fish
    dist_to_cond = pos_fish_light_on(2:end,:) - pos_fish_light_on(1,:);
    dist_to_conditioned_fish_x(:,e) = abs(dist_to_cond(:,1));
    dist_to_conditioned_fish(:,e) = vecnorm(dist_to_cond,2,2);
    [~, rank_pos_id] = sort(vecnorm(dist_to_cond,2,2), 'ascend');
    [~, rank_pos_id_x] = sort(abs(dist_to_cond(:,1)), 'ascend');
    order_to_conditioned_fish_x(:,e) = rank_pos_id_x;
    order_to_conditioned_fish(:,e) = rank_pos_id;
        
    % check which side of the tank fish are in
    if pos_fish_light_on(:,1) < tank_middle

        dist_x = tank_min_x - pos_fish_light_on(:,1); % x_min_wall - x

    elseif pos_fish_light_on(:,1) > tank_middle

        dist_x = pos_fish_light_on(:,1) - tank_max_x; % x - x_max_wall

    end

    dist_y = (tank_width/2) - pos_fish_light_on(:,2); % distance to centre of the tank

    dist_min_wall(:,e) = dist_x;
    dist_cent_wall(:,e) = sqrt(dist_x.^2 + dist_y.^2);

end

%% distribution of distance to conditioned fish

mean_dist_to_cond = mean(dist_to_conditioned_fish, 2); % mean dist to conditioned fish
std_dist_to_cond = std(dist_to_conditioned_fish, 0, 2); % sd
se_dist_to_cond = std_dist_to_cond/sqrt(no_exp); % se

mean_dist_to_cond_x = mean(dist_to_conditioned_fish_x,2);
std_dist_to_cond_x = std(dist_to_conditioned_fish_x, 0, 2);
se_dist_to_cond_x = std_dist_to_cond_x/sqrt(no_exp);

% model

dist_to_conditioned_fish_model = readmatrix('dist_to_cond_fish_model.csv');
dist_to_conditioned_fish_x_model = readmatrix('dist_to_cond_fish_x_model.csv');

no_it = 39*100;
mean_dist_to_cond_model = mean(dist_to_conditioned_fish_model, 2); % mean dist to conditioned fish
std_dist_to_cond_model = std(dist_to_conditioned_fish_model, 0, 2); % sd
se_dist_to_cond_model = (std_dist_to_cond_model)/sqrt(no_it); % se

mean_dist_to_cond_x_model = mean(dist_to_conditioned_fish_x_model,2);
std_dist_to_cond_x_model = std(dist_to_conditioned_fish_x_model, 0, 2);
se_dist_to_cond_x_model = (std_dist_to_cond_x_model)/sqrt(no_it);

plt_count = 1;
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
hold on
errorbar((2:n)+0.075, mean_dist_to_cond_x_model, se_dist_to_cond_x_model, ...
    'LineStyle', "none", "Marker", "s", ...
    "Color", esp_color, 'LineWidth', lw_plot, 'MarkerSize', ...
    10, 'MarkerFaceColor', esp_color)
hold on
errorbar((2:n)+0.075, mean_dist_to_cond_model, se_dist_to_cond_model, ...
    'LineStyle', "none", "Marker", "d", ...
    "Color", esp_color, 'LineWidth', lw_plot, 'MarkerSize', 10, ...
    'MarkerFaceColor', esp_color)

legend({'x^i_{C}', 'r^i_{C}'}, 'Location', 'southwest')
legend('Box', 'off')
set(gca, 'XLim', [1.7 5+0.3], 'XTick', 1:n, 'YLim', [0, 7], ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel("Distance to conditioned fish (cm)", 'FontSize', 23)

% exportgraphics(gca, 'dcf_data.pdf', 'ContentType', 'vector')

%% saving it to csv

% video = 1:no_exp;
% naive_agents = (2:n)';
% naive_agents = repmat(naive_agents, no_exp, 1);
% video = repmat(video, (n-1), 1);
% video = video(:);
% dist_to_cf = dist_to_conditioned_fish(:);
% 
% dist_to_cf_data = [video, naive_agents, dist_to_cf];
% writematrix(dist_to_cf_data, 'dist_to_cf_data.csv')

%% saving distance to tank wall to .csv

% video = 1:no_exp;
% naive_agents = (1:n)';
% naive_agents = repmat(naive_agents, no_exp, 1);
% video = repmat(video, n, 1);
% video = video(:);
% dist_to_tw = dist_cent_wall(:);
% 
% dist_to_tw_data = [video, naive_agents, dist_to_tw];
% writematrix(dist_to_tw_data, 'dist_to_tw_data.csv')

