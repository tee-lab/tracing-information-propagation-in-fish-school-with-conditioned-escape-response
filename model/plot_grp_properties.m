close all
clear
clc

%% Measuring group properties like polarisation, cohesion, etc.

tic

n = 5;

t_edge_id = 1;
mean_id = 2;
std_id = 3;
se_id = 4;
median_id = 5;
fn_st = 4;

font_size = 25;
lw_axis = 2;
lw_plot = 2;
lw_xline = 3;
mr_size = 10;
label_fs = 25;

fl = dir('pd*.csv');
fl = struct2table(fl);
dt = 0.04;

ini_color = "#33b1ff";
esp_color = "#fa4d56";
relax_color = "#198038";
speed_plt_color = ["#6929c4", "#1192e8", "#005d5d", "#9f1853", "#d2a106"];

spd_ini_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/spd_ini_emp.csv');
spd_esp_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/spd_esp_emp.csv');
spd_relax_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/spd_relax_emp.csv');

pol_ini_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/pol_ini_emp.csv');
pol_esp_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/pol_esp_emp.csv');
pol_relax_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/pol_relax_emp.csv');

gc_ini_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/gc_ini_emp.csv');
gc_esp_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/gc_esp_emp.csv');
gc_relax_emp = readmatrix('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/bee_dance_tracks/new_tracks_12_2_25/gc_relax_emp.csv');

for f = 93:1:height(fl)

    fname = fl.name{f};
    fname = fname(fn_st:(end-4));
    % fname = "mual_0.7_muat_1.2_muesp_0.6_k_4_gamma_0.4_omegaini_0.6";
    disp(fname)

    plt_count = 0;

    %% group polarisation

    pol_data = readmatrix(strcat('pd_', fname, '.csv'));
    t_edges = pol_data(:,t_edge_id);
    mean_pol = pol_data(:,mean_id);

    plt_count = plt_count + 1;
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

    exportgraphics(gca, 'pol_ts_model.pdf', 'ContentType', 'vector')

    %% Group cohesion

    grp_coh = readmatrix(strcat('gc_', fname, '.csv'));
    grp_coh_x = readmatrix(strcat('gc_x_', fname, '.csv'));
    grp_coh_y = readmatrix(strcat('gc_y_', fname, '.csv'));
    grp_coh_wcf = readmatrix(strcat('gc_wcf_', fname, '.csv'));

    plt_count = plt_count + 1;
    fig = figure(plt_count); 
    fig.Position = [300, 1200, 800, 700];

    t_edges = grp_coh(:,t_edge_id);
    mean_gc = grp_coh(:,mean_id);
    mean_gc_x = grp_coh_x(:,mean_id);
    mean_gc_y = grp_coh_y(:,mean_id);
    mean_gc_wcf = grp_coh_wcf(:,mean_id);

    se_gc = grp_coh(:,se_id);
    se_gc_x = grp_coh_x(:,se_id);
    se_gc_y = grp_coh_y(:,se_id);
    se_gc_wcf = grp_coh_wcf(:,se_id);

    yyaxis left
    plot(t_edges, mean_pol, 'o-', 'Color', ini_color, 'LineWidth', lw_plot, ...
        'MarkerFaceColor', ini_color)
    set(gca, 'YLim', [0 1], 'LineWidth', lw_axis, ...
        'Xcolor', 'k', 'YColor', 'k', ...
        'FontSize', font_size, 'FontName', 'Helvetica')
    ylabel('Polarisation', 'FontSize', label_fs)

    yyaxis right
    plot(t_edges, mean_gc, 'o-', 'Color', '#A52A2A', 'LineWidth', lw_plot, ...
        'MarkerFaceColor', '#A52A2A')
    hold on
    % plot(t_edges, mean_gc_x, 'o-', 'Color', ini_color, 'LineWidth', lw_plot, ...
    %     'MarkerFaceColor', ini_color)
    % hold on
    % plot(t_edges, mean_gc_y, 'o-', 'Color', relax_color, 'LineWidth', lw_plot, ...
    %     'MarkerFaceColor', relax_color)
    % hold on
    % plot(t_edges, mean_gc_wcf, 'o-', 'Color', 'k', 'LineWidth', lw_plot, 'MarkerFaceColor', 'k')
    hold on
    xline(0, '--r', 'LineWidth', lw_plot)
    hold on
    xline(1, '--r', 'LineWidth', lw_plot)
    hold off

    % errorbar(t_edges, mean_gc, se_gc, 'o-', 'Color', 'r')
    % hold on
    % errorbar(t_edges, mean_gc_x, se_gc_x, 'o-', 'Color', 'b')
    % hold on
    % errorbar(t_edges, mean_gc_y, se_gc_y, 'o-', 'Color', 'g')
    % hold on
    % errorbar(t_edges, mean_gc_wcf, se_gc_wcf, 'o-', 'Color', 'k')
    % hold on
    % xline(0, '--r', 'LineWidth', 1.5)
    % hold on
    % xline(1, '--r', 'LineWidth', 1.5)
    % hold off

    % legend({'C', 'C_x', 'C_y', 'C_{wc}'}, 'Location', 'best')
    % legend('Box', 'off')

    legend({'Polarisation', 'Dispersion'}, 'Location', 'northeast')
    legend('Box', 'off')

    set(gca, 'XLim', [-0.5 3], 'YLim', [0, 7], 'LineWidth', lw_axis, ...
        'Xcolor', 'k', 'YColor', 'k', ...
        'FontSize', font_size, 'FontName', 'Helvetica')

    xlabel('t_n', 'FontSize', label_fs)
    ylabel('Dispersion (cm)', 'FontSize', label_fs)

    exportgraphics(gca, 'pol_disp_ts_model.pdf', 'ContentType', 'vector')

    %% Speed time series

    t_edges = readmatrix(strcat('spd_tedges_data_', fname, '.csv'));
    mean_spd = readmatrix(strcat('mean_spd_data_', fname, '.csv'));
    se_spd = readmatrix(strcat('se_spd_data_', fname, '.csv'));

    plt_count = plt_count + 1;
    fig = figure(plt_count); 
    fig.Position = [300, 1200, 800, 700];

    for i = 1:n

        plot(t_edges, mean_spd(:,i), '-o', 'LineWidth', lw_plot, ...
            'Color', speed_plt_color(i))
        % errorbar(t_edges, mean_spd_i, se_spd_i, 'o')
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

    exportgraphics(gca, 'spd_ts_model_omega_07.pdf', 'ContentType', 'vector')

    %% Crossing time
    
    tc_min = readmatrix(strcat('tc_min_', fname, '.csv'));
    tc_cent = readmatrix(strcat('tc_cent_', fname, '.csv'));
    tc_max = readmatrix(strcat('tc_max_', fname, '.csv'));
    
    st_id = 2;
    mean_td_min = tc_min(st_id:end,mean_id)*dt;
    mean_td_centre = tc_cent(st_id:end,mean_id)*dt;
    mean_td_max = tc_max(st_id:end,mean_id)*dt;

    se_td_min = tc_min(st_id:end,se_id)*dt;
    se_td_centre = tc_cent(st_id:end,se_id)*dt;
    se_td_max = tc_max(st_id:end,se_id)*dt;

    plt_count = plt_count + 1;
    fig = figure(plt_count); 
    fig.Position = [300, 1200, 800, 700];

    errorbar((st_id:n), mean_td_min, se_td_min, 'o', ...
        'LineWidth', lw_plot, 'Color', ini_color, 'MarkerSize', mr_size, ...
        'MarkerFaceColor', ini_color)
    % hold on
    % errorbar((st_id:n), mean_td_centre, se_td_centre, 'o', ...
    %     'LineWidth', lw_plot, 'Color', esp_color, 'MarkerSize', mr_size, ...
    %     'MarkerFaceColor', esp_color)
    % hold on
    % errorbar((st_id:n)+0.1, mean_td_max, se_td_max, 'o', ...
    %     'LineWidth', lw_plot, 'Color', relax_color, 'MarkerSize', mr_size, ...
    %     'MarkerFaceColor', relax_color)
    % 
    % legend({'barrier-start', 'barrier-centre', 'barrier-last'}, 'Location', 'best')
    % legend('Box', 'off')

    ymax_lim = ceil(max(mean_td_min + se_td_min));
    set(gca, 'XLim', [st_id-0.4 n+0.2], 'XTick', 1:n, 'YLim', [0,ymax_lim], ...
        'YTick', 0:1:ymax_lim, 'LineWidth', lw_axis, ...
        'Xcolor', 'k', 'YColor', 'k', ...
        'FontSize', font_size, 'FontName', 'Helvetica')
    
    xlabel('Crossing rank', 'FontSize', label_fs)
    ylabel('Time since previous fish crossed (s)', 'FontSize', label_fs)
    exportgraphics(gca, 'tspfc_model.pdf', 'ContentType', 'vector')

    tslfc_min_data = readmatrix(strcat('tc_min_all_data', fname, '.csv'));

    no_it = round(size(tslfc_min_data,1)/5);
    video = 1:no_it;
    video = repmat(video, n+1-st_id, 1);
    video = video(:);

    naive_agents = (st_id:n)';
    naive_agents = repmat(naive_agents, no_it, 1);
    naive_agents_id = tslfc_min_data(:,1) ~= 1;

    lmer_tslfc_data = [video, naive_agents, tslfc_min_data(naive_agents_id,2)*dt];
    % lmer_tslfc_data = [video, naive_agents, tslfc_min_data(:,2)*dt];
    tslfc_fn = strcat('lmer_tslfc_', fname, '.csv');
    writematrix(lmer_tslfc_data, tslfc_fn);

    %% distribution of distance to conditioned fish

    % reading files

    dist_min_wall = readmatrix(strcat('dist_min_wall_', fname, '.csv'));
    dist_cent_wall = readmatrix( strcat('dist_cent_wall_', fname, '.csv'));
    dist_to_conditioned_fish = readmatrix(strcat('dist_to_cond_fish_', fname, '.csv'));
    dist_to_conditioned_fish_x = readmatrix(strcat('dist_to_cond_fish_x_', fname, '.csv'));
    % order_to_conditioned_fish = readmatrix(strcat('order_to_cond_fish_', fname, '.csv'));
    % order_to_conditioned_fish_x = readmatrix(strcat('order_to_cond_fish_x_', fname, '.csv'));

    no_it = size(dist_min_wall, 2);
    video = 1:no_it;
    video = repmat(video, n, 1);
    video = video(:);
    naive_agents = (1:n)';
    naive_agents = repmat(naive_agents, no_it, 1);

    lmer_dttw_data = [video, naive_agents, dist_cent_wall(:)];
    dttw_fn = strcat('lmer_dttw_', fname, '.csv');
    writematrix(lmer_dttw_data, dttw_fn);

    psi_ic = readmatrix(strcat('psi_ic', fname, '.csv'));
    phi_ic = readmatrix(strcat('phi_ic', fname, '.csv'));

    no_it = size(dist_min_wall, 2);
    video = 1:no_it;
    video = repmat(video, n-1, 1);
    video = video(:);
    naive_agents = (2:n)';
    naive_agents = repmat(naive_agents, no_it, 1);
    
    lmer_dtcf_data = [video, naive_agents, dist_to_conditioned_fish(:)];
    dtcf_fn = strcat('lmer_dtcf_', fname, '.csv');
    writematrix(lmer_dtcf_data, dtcf_fn);

    lmer_psi_ic_data = [video, naive_agents, psi_ic(:)];
    psi_ic_fn = strcat('lmer_psi_', fname, '.csv');
    writematrix(lmer_psi_ic_data, psi_ic_fn)

    lmer_phi_ic_data = [video, naive_agents, phi_ic(:)];
    phi_ic_fn = strcat('lmer_phi_', fname, '.csv');
    writematrix(lmer_phi_ic_data, phi_ic_fn)

    mean_dist_to_cond = mean(dist_to_conditioned_fish, 2); % mean dist to conditioned fish
    std_dist_to_cond = std(dist_to_conditioned_fish, 0, 2); % sd
    se_dist_to_cond = (std_dist_to_cond)/sqrt(no_it); % se

    mean_dist_to_cond_x = mean(dist_to_conditioned_fish_x,2);
    std_dist_to_cond_x = std(dist_to_conditioned_fish_x, 0, 2);
    se_dist_to_cond_x = (std_dist_to_cond_x)/sqrt(no_it);
    
    psi_ic = (180*psi_ic)/pi;
    mean_psi_ic = mean(psi_ic, 2);
    std_psi_ic = std(psi_ic, 0, 2);
    se_psi_ic = (std_psi_ic)/sqrt(no_it);

    phi_ic = (180*phi_ic)/pi;
    mean_phi_ic = mean(phi_ic, 2);
    std_phi_ic = std(phi_ic, 0, 2);
    se_phi_ic = (std_phi_ic)/sqrt(no_it);

    plt_count = plt_count + 1;
    fig = figure(plt_count);
    fig.Position = [300, 1200, 800, 700];

    errorbar(2:n, mean_dist_to_cond_x, se_dist_to_cond_x, 'LineStyle', "none", ...
        "Marker", "o", "Color", ini_color, 'LineWidth', lw_plot, ...
        'MarkerSize', mr_size, 'MarkerFaceColor', ini_color)
    hold on
    errorbar(2:n, mean_dist_to_cond, se_dist_to_cond, 'LineStyle', "none", ...
        "Marker", "o", "Color", esp_color, 'LineWidth', lw_plot, ...
        'MarkerSize', mr_size, 'MarkerFaceColor', esp_color)
        
    legend({'x^i_{C}', 'r^i_{C}'}, 'Location', 'southeast')
    legend('Box', 'off')
    
    ylim_max = ceil(max(mean_dist_to_cond + se_dist_to_cond));
    set(gca, 'XLim', [1.9 5.1], 'XTick', 2:n, 'YLim', [0, ylim_max], ...
        'YTick', 0:1:ylim_max, 'LineWidth', lw_axis, ...
        'Xcolor', 'k', 'YColor', 'k', ...
        'FontSize', font_size, 'FontName', 'Helvetica')

    xlabel("Crossing rank", 'FontSize', label_fs)
    ylabel("Distance to conditioned fish (cm)", 'FontSize', label_fs)

    exportgraphics(gca, 'dcf_model.pdf', 'ContentType', 'vector')
    % set(gca, 'XLim', [1.9 5.1], 'XTick', 2:n, 'LineWidth', lw_axis, ...
    %     'Xcolor', 'k', 'YColor', 'k', ...
    %     'FontSize', font_size, 'FontWeight', 'bold')

    % relative orientation

    plt_count = plt_count + 1;
    fig = figure(plt_count);
    fig.Position = [300, 1200, 800, 700];
    errorbar(2:n, mean_psi_ic, se_psi_ic, 'LineStyle', "none", ...
        "Marker", "o", "Color", ini_color, 'LineWidth', lw_plot, ...
        'MarkerSize', mr_size, 'MarkerFaceColor', ini_color)

    ylim_min = floor(min(mean_psi_ic - se_psi_ic)) - 10;
    ylim_max = ceil(max(mean_psi_ic + se_psi_ic)) + 10;
    set(gca, 'XLim', [1.9 5.1], 'XTick', 2:n, 'YLim', [0, 120], ...
        'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
        'FontSize', font_size, 'FontName', 'Helvetica')

    xlabel("Crossing rank", 'FontSize', label_fs)
    ylabel("Viewing angle (in degree)", 'FontSize', label_fs)
    
    exportgraphics(gca, 'psi_model.pdf', 'ContentType', 'vector')

    plt_count = plt_count + 1;
    fig = figure(plt_count);
    fig.Position = [300, 1200, 800, 700];
    errorbar(2:n, mean_phi_ic, se_phi_ic, 'LineStyle', "none", ...
        "Marker", "o", "Color", ini_color, 'LineWidth', lw_plot, ...
        'MarkerSize', mr_size, 'MarkerFaceColor', ini_color)
    
    ylim_min = floor(min(mean_phi_ic - se_phi_ic)) - 10;
    ylim_max = ceil(max(mean_phi_ic + se_phi_ic)) + 10;
    set(gca, 'XLim', [1.9 5.1], 'XTick', 2:n, 'YLim', [0, 100],...
        'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
        'FontSize', font_size, 'FontName', 'Helvetica')

    xlabel("Crossing rank", 'FontSize', label_fs)
    ylabel("Relative orientation (in degree)", 'FontSize', label_fs)

    exportgraphics(gca, 'phi_model.pdf', 'ContentType', 'vector')

    % plt_count = plt_count + 1;
    % fig = figure(plt_count);
    % fig.Position = [300, 1200, 800, 700];
    % errorbar(2:n, mean_dist_to_cond_x, se_dist_to_cond_x, 'LineStyle', "none", "Marker", "o", ...
    %     "Color", "k")
    % xlim([1.5 5.5])
    % ylabel("Distance to conditioned fish (x)")


    %% distribution of distance to min tank x

    % plt_count = plt_count + 1;
    % fig = figure(plt_count); 
    % fig.Position = [300, 1200, 800, 700];

    mean_dist_min_wall_x = mean(dist_min_wall, 2); % mean distance to closest x-wall
    std_dist_min_wall_x = std(dist_min_wall, 0, 2);
    se_dist_min_wall_x = std_dist_min_wall_x/sqrt(no_it);

    % for i = 1:n
    %
    %     scatter(i, dist_min_wall(i,:), sz, "k", 'filled')
    %     hold on
    %
    % end
    % ylabel("Distance to tank wall")
    % xlim([0.5 5.5])

    % plt_count = plt_count + 1;
    % fig = figure(plt_count); 
    % fig.Position = [300, 1200, 800, 700];
    % 
    % errorbar(1:n, mean_dist_min_wall_x, se_dist_min_wall_x, 'LineStyle', "none", "Marker", "o", ...
    %     "Color", "k")
    % xlim([0.5 5.5])
    % ylabel("Distance to tank wall (x)")

    %% distribution of distance to tank centre

    plt_count = plt_count + 1;
    fig = figure(plt_count);
    fig.Position = [300, 1200, 800, 700];

    mean_dist_to_cent = mean(dist_cent_wall, 2); % distance to centre of closest x-wall
    std_dist_to_cent = std(dist_cent_wall, 0, 2);
    se_dist_to_cent = std_dist_to_cent/sqrt(no_it);

    errorbar(2:n, mean_dist_to_cond, se_dist_to_cond, 'LineStyle', "none", ...
        "Marker", "o", "Color", ini_color, 'LineWidth', lw_plot, ...
        'MarkerSize', mr_size, 'MarkerFaceColor', ini_color)
    hold on
    errorbar(1:n, mean_dist_to_cent, se_dist_to_cent, 'LineStyle', "none", ...
        "Marker", "o", "Color", esp_color, 'LineWidth', lw_plot, ...
        'MarkerSize', mr_size, 'MarkerFaceColor', esp_color)

    legend({'DCF', 'DTW'}, 'Location', 'best')
    legend('Box','off')

    set(gca, 'XLim', [0.5 5.5], 'XTick', 1:n, 'YLim', [0,12],...
        'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
        'FontSize', font_size, 'FontName', 'Helvetica')
    
    xlabel("Crossing rank", 'FontSize', label_fs)
    ylabel("Distance at t = 0 (cm)", 'FontSize', label_fs)

    exportgraphics(gca, 'dcf_tw_model.pdf', 'ContentType', 'vector')

    % plotting distributions of speed

    % speed_ini_all = readmatrix(strcat('s_ini_', fname, '.csv'));
    % speed_escape_all = readmatrix(strcat('s_esp_', fname, '.csv'));
    % speed_relax_all = readmatrix(strcat('s_relax_', fname, '.csv'));
    % 
    % [spd_ini_hist, spd_ini_edges] = histcounts(speed_ini_all, 'Normalization', 'pdf');
    % [spd_esp_hist, spd_esp_edges] = histcounts(speed_escape_all, 'Normalization', 'pdf');
    % [spd_relax_hist, spd_relax_edges] = histcounts(speed_relax_all, 'Normalization', 'pdf');

    spd_ini = readmatrix(strcat('spd_ini_hist_', fname, '.csv'));
    spd_ini_edges = spd_ini(:,1);
    spd_ini_hist = spd_ini(:,2);
    spd_esp = readmatrix(strcat('spd_esp_hist_', fname, '.csv'));
    spd_esp_edges = spd_esp(:,1);
    spd_esp_hist = spd_esp(:,2);
    spd_relax = readmatrix(strcat('spd_relax_hist_', fname, '.csv'));
    spd_relax_edges = spd_relax(:,1);
    spd_relax_hist = spd_relax(:,2);

    [spd_ini_hist_emp, spd_ini_edges_emp] = histcounts(spd_ini_emp, 'Normalization', 'pdf');
    [spd_esp_hist_emp, spd_esp_edges_emp] = histcounts(spd_esp_emp, 'Normalization', 'pdf');
    [spd_relax_hist_emp, spd_relax_edges_emp] = histcounts(spd_relax_emp, 'Normalization', 'pdf');

    plt_count = plt_count + 1;
    fig = figure(plt_count);
    fig.Position = [300, 1200, 800, 700];

    plot(spd_ini_edges_emp(1:end-1), spd_ini_hist_emp,  '-', ...
        'Color', ini_color, 'LineWidth', lw_xline)
    hold on
    plot(spd_esp_edges_emp(1:end-1), spd_esp_hist_emp, '-', ...
        'Color', esp_color, 'LineWidth', lw_xline)
    hold on
    plot(spd_relax_edges_emp(1:end-1), spd_relax_hist_emp, '-', ...
        'Color', relax_color, 'LineWidth', lw_xline)
    hold on
    plot(spd_ini_edges, spd_ini_hist, '--', 'Color', ini_color, 'LineWidth', lw_xline)
    hold on
    plot(spd_esp_edges, spd_esp_hist, '--', 'Color', esp_color, 'LineWidth', lw_xline)
    hold on
    plot(spd_relax_edges, spd_relax_hist, '--', 'Color', relax_color, 'LineWidth', lw_xline)
    
    % xline(median(speed_ini_all), '--k', 'LineWidth', 2)
    % hold on
    % xline(median(speed_escape_all), '--r', 'LineWidth', 2)
    % hold on
    % xline(median(speed_relax_all), '--g', 'LineWidth', 2)

    set(gca, 'LineWidth', lw_axis, 'XLim', [0.03,8], ...
         'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
         'FontSize', font_size, 'FontName', 'Helvetica')

    xlabel('Speed (cm/s)')
    ylabel('PDF')

    % legend({'Initial', 'Escape', 'Relax'}, 'Location', 'best')
    % legend('Box', 'off')

    exportgraphics(gca, 'spd_pdf_model.pdf', 'ContentType', 'vector')

    % plotting pdf of polarisation

    % p_ini = readmatrix(strcat('p_ini_', fname, '.csv'));
    % p_esp = readmatrix(strcat('p_esp_', fname, '.csv'));
    % p_relax = readmatrix(strcat('p_relax_', fname, '.csv'));

    p_ini = readmatrix(strcat('pol_ini_hist_', fname, '.csv'));
    p_ini_edges = p_ini(:,1);
    p_ini_hist = p_ini(:,2);
    p_esp = readmatrix(strcat('pol_esp_hist_', fname, '.csv'));
    p_esp_edges = p_esp(:,1);
    p_esp_hist = p_esp(:,2);
    p_relax = readmatrix(strcat('pol_relax_hist_', fname, '.csv'));
    p_relax_edges = p_relax(:,1);
    p_relax_hist = p_relax(:,2);

    p_edges = 41;

    % [p_ini_hist, p_ini_edges] = histcounts(p_ini, p_edges, 'Normalization', 'pdf');
    % [p_esp_hist, p_esp_edges] = histcounts(p_esp, p_edges, 'Normalization', 'pdf');
    % [p_relax_hist, p_relax_edges] = histcounts(p_relax, p_edges, 'Normalization', 'pdf');

    [p_ini_hist_emp, p_ini_edges_emp] = histcounts(pol_ini_emp, p_edges, 'Normalization', 'pdf');
    [p_esp_hist_emp, p_esp_edges_emp] = histcounts(pol_esp_emp, p_edges, 'Normalization', 'pdf');
    [p_relax_hist_emp, p_relax_edges_emp] = histcounts(pol_relax_emp, p_edges, 'Normalization', 'pdf');
    
    plt_count = plt_count + 1;
    fig = figure(plt_count);
    fig.Position = [300, 1200, 800, 700];

    plot(p_ini_edges_emp(1:end-1), p_ini_hist_emp, '-', ...
        'Color', ini_color, 'LineWidth', lw_xline)
    hold on
    plot(p_esp_edges_emp(1:end-1), p_esp_hist_emp, '-', ...
        'Color', esp_color, 'LineWidth', lw_xline)
    hold on
    plot(p_relax_edges_emp(1:end-1), p_relax_hist_emp, '-', ...
        'Color', relax_color, 'LineWidth', lw_xline)
    hold on
    plot(p_ini_edges, p_ini_hist, '--', 'Color', ini_color, 'LineWidth', lw_xline)
    hold on
    plot(p_esp_edges, p_esp_hist, '--', 'Color', esp_color, 'LineWidth', lw_xline)
    hold on
    plot(p_relax_edges, p_relax_hist, '--', 'Color', relax_color, 'LineWidth', lw_xline)
    
    % hold on
    % xline(median(p_ini), '--k', 'LineWidth', 2)
    % hold on
    % xline(median(p_esp), '--r', 'LineWidth', 2)
    % hold on
    % xline(median(p_relax), '--g', 'LineWidth', 2)

    set(gca, 'LineWidth', lw_axis, 'XLim', [0,1], ...
        'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
         'FontSize', font_size, 'FontName', 'Helvetica')

    xlabel('Polarisation')
    ylabel('PDF')

    % legend({'Initial', 'Escape', 'Relax'}, 'Location', 'best')
    % legend('Box', 'off')

    exportgraphics(gca, 'pol_pdf_model.pdf', 'ContentType', 'vector')

    % plotting distribution of group cohesion

    % gc_ini = readmatrix(strcat('gc_ini_', fname, '.csv'));
    % gc_esp = readmatrix(strcat('gc_esp_', fname, '.csv'));
    % gc_relax = readmatrix(strcat('gc_relax_', fname, '.csv'));

    gc_ini = readmatrix(strcat('gc_ini_hist_', fname, '.csv'));
    gc_ini_edges = gc_ini(:,1);
    gc_ini_hist = gc_ini(:,2);
    gc_esp = readmatrix(strcat('gc_esp_hist_', fname, '.csv'));
    gc_esp_edges = gc_esp(:,1);
    gc_esp_hist = gc_esp(:,2);
    gc_relax = readmatrix(strcat('gc_relax_hist_', fname, '.csv'));
    gc_relax_edges = gc_relax(:,1);
    gc_relax_hist = gc_relax(:,2);

    % [gc_ini_hist, gc_ini_edges] = histcounts(gc_ini, 'Normalization', 'pdf');
    % [gc_esp_hist, gc_esp_edges] = histcounts(gc_esp, 'Normalization', 'pdf');
    % [gc_relax_hist, gc_relax_edges] = histcounts(gc_relax, 'Normalization', 'pdf');

    [gc_ini_hist_emp, gc_ini_edges_emp] = histcounts(gc_ini_emp, 'Normalization', 'pdf');
    [gc_esp_hist_emp, gc_esp_edges_emp] = histcounts(gc_esp_emp, 'Normalization', 'pdf');
    [gc_relax_hist_emp, gc_relax_edges_emp] = histcounts(gc_relax_emp, 'Normalization', 'pdf');
    
    plt_count = plt_count + 1;
    fig = figure(plt_count);
    fig.Position = [300, 1200, 800, 700];

    plot(gc_ini_edges_emp(1:end-1), gc_ini_hist_emp, '-', ...
        'Color', ini_color, 'LineWidth', lw_xline)
    hold on
    plot(gc_esp_edges_emp(1:end-1), gc_esp_hist_emp, '-', ...
        'Color', esp_color, 'LineWidth', lw_xline)
    hold on
    plot(gc_relax_edges_emp(1:end-1), gc_relax_hist_emp, '-', ...
        'Color', relax_color, 'LineWidth', lw_xline)
    hold on
    plot(gc_ini_edges, gc_ini_hist, '--', 'Color', ini_color, 'LineWidth', lw_xline)
    hold on
    plot(gc_esp_edges, gc_esp_hist, '--', 'Color', esp_color, 'LineWidth', lw_xline)
    hold on
    plot(gc_relax_edges, gc_relax_hist, '--', 'Color', relax_color, 'LineWidth', lw_xline)
    
    % hold on
    % xline(median(gc_ini), '--k', 'LineWidth', 2)
    % hold on
    % xline(median(gc_esp), '--r', 'LineWidth', 2)
    % hold on
    % xline(median(gc_relax), '--g', 'LineWidth', 2)

    set(gca, 'LineWidth', lw_axis, 'XLim', [0,12], ...
        'LineWidth', lw_axis, 'Xcolor', 'k', 'YColor', 'k', ...
         'FontSize', font_size, 'FontName', 'Helvetica')

    xlabel('Dispersion (cm)')
    ylabel('PDF')

    legend({'Initial', 'Escape', 'Relax'}, 'Location', 'north')
    legend('Box', 'off')

    exportgraphics(gca, 'disp_pdf_model.pdf', 'ContentType', 'vector')

    close all
    % plt_count = plt_count + 1;
    % figure(plt_count)
    %
    % for i = 1:n
    %
    %     scatter(i, dist_cent_wall(i,:), sz, "k", 'filled')
    %     hold on
    %
    % end
    % xlim([0.5 5.5])
    % ylabel("Distance to tank centre")

    %% rank order of crossing and rank order of closest neighbors

    % plt_count = plt_count + 1;
    % figure(plt_count)
    % 
    % corr_order_cross_neigh = nan(n-1,n-1);
    % no_exp = size(order_to_conditioned_fish,2);
    % 
    % for j = 1:(n-1)
    % 
    %     for i = 1:(n-1)
    % 
    %         corr_order_cross_neigh(i,j) = (sum(order_to_conditioned_fish(j,:) == i))/no_exp;
    % 
    %     end
    % 
    % end
    % 
    % for i = 1:(n-1)
    % 
    %     subplot(2,2,i)
    %     scatter(1:(n-1), corr_order_cross_neigh(:,i), "Color", 'k')
    %     title(num2str(i))
    %     xlabel('Rank pos order')
    %     ylabel('Fraction of time crossing (2d)')
    % 
    % end
    % 
    % plt_count = plt_count + 1;
    % figure(plt_count)
    % corr_order_cross_neigh_x = nan(n-1,n-1);
    % 
    % for j = 1:(n-1)
    % 
    %     for i = 1:(n-1)
    % 
    %         corr_order_cross_neigh_x(i,j) = (sum(order_to_conditioned_fish_x(j,:) == i))/no_exp;
    % 
    %     end
    % 
    % end
    % 
    % for i = 1:(n-1)
    % 
    %     subplot(2,2,i)
    %     scatter(1:(n-1), corr_order_cross_neigh_x(:,i), "Color", 'k')
    %     title(num2str(i))
    %     xlabel('Rank pos order')
    %     ylabel('Fraction of time crossing (x)')
    % end

end