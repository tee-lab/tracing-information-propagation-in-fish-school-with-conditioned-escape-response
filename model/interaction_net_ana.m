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

%% constructing full network

if ~isempty(full_net)

    fn_avg_edges = full_net(row_avg_edges, :);
    fn_sig_edge_id = fn_avg_edges > min_edges;
    fn_avg_edges = fn_avg_edges(fn_sig_edge_id);

    fn_leader = full_net(row_leader, fn_sig_edge_id);
    fn_follower = full_net(row_follower, fn_sig_edge_id);
    fn_avg_cc = full_net(row_avg_cc, fn_sig_edge_id);


    full_int_graph = digraph(fn_follower, fn_leader, round(fn_avg_edges,3));

    fig = figure(plt_count);
    fig.Position = [300, 1200, 800, 700];
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


