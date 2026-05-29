close all
clear
clc

%% load files

fname = "sim";
method = "_pm_"; % stat or pm
phase = 2; % 1, 2, 3
fnet_name = strcat('full_net', method, 'phase_', num2str(phase), '_', fname, '.csv');
dnet_name = strcat('diff_net', method, 'phase_', num2str(phase), '_', fname, '.csv');

full_net = readmatrix(fnet_name);
% diff_net = readmatrix(dnet_name);
 
row_leader = 1;
row_follower = 2;
row_avg_cc = 3;
row_avg_edges = 4;

min_edges = 0.05;

plt_count = 1;

ini_color = "#0096FF";
esp_color = "#EE4B2B";
relax_color = "#097969";
font_size = 25;

%% Check if the time series are stationary

if method == "_stat_"

    is_stat_fname = strcat('is_stat_phase', num2str(1), '_', fname, '.csv');
    is_stat = readmatrix(is_stat_fname);

    figure(plt_count)
    plt_count = plt_count + 1;

    histogram(is_stat, 'Normalization', 'pdf')
    h = gca;
    h.XTick = [0, 1];
    h.XTickLabel = {'Non-Stat', 'Stat'};
    ylabel('PDF')

end

%% constructing full network

if ~isempty(full_net)

    fn_avg_edges = full_net(row_avg_edges, :);
    fn_sig_edge_id = fn_avg_edges > min_edges;
    fn_avg_edges = fn_avg_edges(fn_sig_edge_id);

    fn_leader = full_net(row_leader, fn_sig_edge_id);
    fn_follower = full_net(row_follower, fn_sig_edge_id);
    fn_avg_cc = full_net(row_avg_cc, fn_sig_edge_id);


    full_int_graph = digraph(fn_follower, fn_leader, round(fn_avg_edges,3));

    figure(plt_count)
    plt_count = plt_count + 1;
    plot(full_int_graph, 'LineWidth', (full_int_graph.Edges.Weight)*4, ...
        'EdgeLabel', round(full_int_graph.Edges.Weight,2),...
        'EdgeFontSize', 10+full_int_graph.Edges.Weight, ...
        'Layout', 'layered', 'EdgeFontWeight', 'bold', ...
        'ArrowSize', 20, 'EdgeColor', esp_color, 'NodeColor', esp_color)

end

% exportgraphics(gca, 'pert_test_net_model.pdf', 'ContentType', 'vector')

%% constructing diff network

% if ~isempty(diff_net)
% 
%     dn_avg_edges = diff_net(row_avg_edges, :);
%     dn_sig_edge_id = dn_avg_edges > min_edges;
%     dn_avg_edges = dn_avg_edges(dn_sig_edge_id);
% 
%     dn_leader = diff_net(row_leader, dn_sig_edge_id);
%     dn_follower = diff_net(row_follower, dn_sig_edge_id);
%     dn_avg_cc = diff_net(row_avg_cc, dn_sig_edge_id);
% 
% 
%     diff_int_graph = digraph(dn_follower, dn_leader, round(dn_avg_edges,3));
% 
%     figure(plt_count)
%     plt_count = plt_count + 1;
%     plot(diff_int_graph, 'LineWidth', diff_int_graph.Edges.Weight + 1, ...
%         'EdgeLabel', diff_int_graph.Edges.Weight, 'Layout', 'layered')
% 
% end

%%

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [300, 1200, 800, 700];
omega_ini = 0.7:0.1:1;
int_strength = [.67, 0.23, 0.21, 0.07];

plot(omega_ini, int_strength, 'o-', 'Color', '#A52A2A', 'LineWidth', 4, ...
        'MarkerFaceColor', '#A52A2A')

set(gca, 'XLim', [0.69 1.01], 'YLim', [0, 1], 'YTick', 0:0.2:1, ...
    'LineWidth', 2, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', 25, 'FontName', 'Helvetica')

xlabel('\omega(0)', 'FontSize', 25)
ylabel('Leadership consistency (1 -> 2)', 'FontSize', 25)

% exportgraphics(gca, 'leadership_consistency.pdf', 'ContentType', 'vector')

