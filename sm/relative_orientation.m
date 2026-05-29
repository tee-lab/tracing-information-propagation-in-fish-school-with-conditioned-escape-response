% code to calculate \psi_{iC} and \phi_{iC}

close all
clear
clc

%% load data

all_bd_data = load('../main_text/all_bd_data.mat'); % load all pos, vel data
tcross_data = load('../main_text/tcross_data.mat');
load("../main_text/all_bd_data.mat", "dt", "n", "no_exp", "tank_max_x", "tank_min_x", "tank_middle")
smooth_window = 5;
tank_width = 20;

t_cut_min = -0.5;
t_cut_max = 3;

font_size = 25;
lw_plot = 2;
lw_axis = 2;

ini_color = "#33b1ff";
esp_color = "#fa4d56";
relax_color = "#198038";

%% calculate relative orientation and viewing angle

psi_ic = nan(n-1, no_exp); % viewing angle when the light was turned on
phi_ic = nan(n-1, no_exp); % relative orientation when light was turned on

for e = 1:no_exp

    % fish data for this given experiment
    pos_fish = all_bd_data.(strcat('pos_t_ex_', num2str(e))); % input what pos data to use.
    vel_fish = all_bd_data.(strcat('vel_t_ex_', num2str(e))); % input velocity data. 
    time_light_on = all_bd_data.(strcat('frame_light_on_ex_', num2str(e))); % when was the light turned on
    rank_id_min = tcross_data.(strcat('rank_id_min', num2str(e))); % rank of crossing
    
    % position when light was turned on in the order of crossing
    pos_fish_light_on = pos_fish(rank_id_min,:,time_light_on);
    vel_fish_light_on = vel_fish(rank_id_min,:,time_light_on);
    vel_naive = vel_fish_light_on(2:end,:);
    vel_naive = vel_naive./(vecnorm(vel_naive,2,2) + eps);
    % r_ic
    r_ic = pos_fish_light_on(1,:) - pos_fish_light_on(2:end,:);
    r_ic = r_ic./(vecnorm(r_ic,2,2) + eps);
    
    % calculating viewing angle
    psi_ic_temp = dot(vel_naive, r_ic, 2);
    psi_ic(:,e) = acos(psi_ic_temp);
    % psi_ic_sign = cross([vel_naive, zeros(n-1,1)], [r_ic, zeros(n-1,1)], 2);
    % psi_ic_temp = psi_ic_temp.*sign(psi_ic_sign(:,3));

    % vel conditioned fish
    vel_con = vel_fish_light_on(1,:);
    vel_con = vel_con./vecnorm(vel_con,2,2);
    vel_con = repmat(vel_con, n-1, 1);
    
    % calculating relative orientation
    phi_ic_temp = dot(vel_naive, vel_con, 2);
    phi_ic(:,e) = acos(phi_ic_temp);
    
end

%% plotting \psi

psi_ic = psi_ic';
psi_ic = psi_ic*180/pi;

mean_psi_ic = mean(psi_ic,1);
std_psi_ic = std(psi_ic,1);
error_psi_ic = (std_psi_ic)/sqrt(no_exp);

% model
no_it = 39*100;
psi_ic_model = readmatrix('psi_ic_model.csv');
psi_ic_model = (180*psi_ic_model)/pi;
mean_psi_ic_model = mean(psi_ic_model, 2);
std_psi_ic_model = std(psi_ic_model, 0, 2);
se_psi_ic_model = (std_psi_ic_model)/sqrt(no_it);

plt_count = 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

errorbar((2:n) - 0.05, mean_psi_ic, error_psi_ic, "o", 'MarkerEdgeColor',  ini_color, ...
    'MarkerFaceColor',  ini_color, 'color', ini_color, 'LineWidth', lw_plot, ...
    'MarkerSize', 10)
hold on
errorbar((2:n) + 0.05, mean_psi_ic_model, se_psi_ic_model, ...
    "o", 'MarkerEdgeColor',  esp_color, ...
    'MarkerFaceColor',  esp_color, 'color', esp_color, 'LineWidth', lw_plot, ...
    'MarkerSize', 10)

legend({'Experiment', 'Model'}, 'Location', 'southwest')
legend('Box', 'off')

set(gca, 'XLim', [1.8 n+0.2], 'YLim', [0,120], 'XTick', 1:n, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel('Viewing angle (in degree)', 'FontSize', 23)

% exportgraphics(gca, 'psi_data.pdf', 'ContentType', 'vector')

%% plotting relative orientation

phi_ic = phi_ic';
phi_ic = phi_ic*180/pi;

mean_phi_ic = mean(phi_ic,1);
std_phi_ic = std(phi_ic,1);
error_phi_ic = (std_phi_ic)/sqrt(no_exp);

% model
phi_ic_model = readmatrix('phi_ic_model.csv');
phi_ic_model = (180*phi_ic_model)/pi;
mean_phi_ic_model = mean(phi_ic_model, 2);
std_phi_ic_model = std(phi_ic_model, 0, 2);
se_phi_ic_model = (std_phi_ic_model)/sqrt(no_it);

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];

errorbar((2:n)-0.05, mean_phi_ic, error_phi_ic, "o", 'MarkerEdgeColor',  ini_color, ...
    'MarkerFaceColor',  ini_color, 'Color', ini_color, 'LineWidth', lw_plot, ...
    'MarkerSize', 10)
hold on
errorbar((2:n)+0.05, mean_phi_ic_model, se_phi_ic_model, "o", 'MarkerEdgeColor',  esp_color, ...
    'MarkerFaceColor',  esp_color, 'Color', esp_color, 'LineWidth', lw_plot, ...
    'MarkerSize', 10)

% legend({'Experiment', 'Model'}, 'Location', 'northwest')
% legend('Box', 'off')

set(gca, 'XLim', [1.8 n+0.2], 'YLim', [0,120], 'XTick', 1:n, ...
    'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', font_size, 'FontName', 'Helvetica')

xlabel('Crossing rank', 'FontSize', 23)
ylabel('Relative orientation (in degree)', 'FontSize', 23)

% exportgraphics(gca, 'phi_data.pdf', 'ContentType', 'vector')
   
%% saving it to csv

% video = 1:no_exp;
% naive_agents = (2:n)';
% naive_agents = repmat(naive_agents, no_exp, 1);
% video = repmat(video, (n-1), 1);
% video = video(:);
% psi_ic = psi_ic(:);
% phi_ic = phi_ic(:);
% 
% psi_data = [video, naive_agents, psi_ic];
% writematrix(psi_data, 'psi_data.csv')
% 
% phi_data = [video, naive_agents, phi_ic];
% writematrix(phi_data, 'phi_data.csv')