%% Code to calculate speed correlation in null-model (without detrending)

close all
clear
clc

%% load data

% tic

fname = dir('*.mat');
fname = fname.name;
load(fname) % load data
fname = fname(1:end-4);

add_enc_time_before = 0;
add_enc_time_after = 0;
max_lag = 75;
min_lag = 0; % min lag for cc to be considered meaningful

% total no.of phase = 3; % 1 = initial, 2 = escape and 3 = relax
no_cal_corr_event = 3; 

%% identifying leadership 

for cal_corr_event = 1:no_cal_corr_event

    % check if the time-series is stationary.
    is_stationary = nan(no_exp*no_it*n, 1);
    is_stat_count = 1;

    avg_exp_net = []; % average interaction network from different sets of simulation
    % leader; follower; cc; no.of.edges

    for e = 1:no_exp

        disp(e)

        % for a given experiment, array of leaders and followers.
        spd_leader_esp = []; % leader id
        spd_follower_esp = []; % follower id
        spd_corr = []; % cc

        for iter = 1:no_it

            % load speed data for given iter and exp
            s_iter = s_t(rank_order_atk(:,iter,e),:,iter,e);

            if cal_corr_event == 1 % initial phase
                en_st_t = 1; % perturbation start time.
                en_end_t = min(t_atk) - 1;
                % en_end_t = min(en_start(:,iter,e)) - 1; % escape time.
                % en_end_t = en_end_t + add_enc_time_after;
            elseif cal_corr_event == 2 % escape phase
                % en_st_t = min(en_start(:,iter,e)); % perturbation start time.
                en_st_t = min(t_atk);
                en_end_t = max(en_end(:,iter,e)); % escape time.
                en_end_t = en_end_t + add_enc_time_after;
                % en_start = max(1, en_start - add_enc_time_before); % some time before start
            elseif cal_corr_event == 3 % relax phase
                en_st_t = max(en_end(:,iter,e)) + 1; % escape time.
                en_st_t = en_st_t + add_enc_time_after;
                en_end_t = size(s_iter,2);
            end

            s_iter = s_iter(:,en_st_t:en_end_t); % speed data during escape.

            for i = 1:n

                spd_i = s_iter(i,:); % speed of fish i
                is_stationary(is_stat_count) = adftest(spd_i);
                is_stat_count = is_stat_count + 1;
                spd_i = spd_i - mean(spd_i, 'omitmissing'); % s - \bar{s}

                for j = 1:n

                    if j ~= i

                        % similarly as above for fish, j
                        spd_j = s_iter(j,:);
                        spd_j = spd_j - mean(spd_j, 'omitmissing');

                        % calculate correlations.
                        [cor_u, cor_ut] = xcorr(spd_i, spd_j, max_lag, 'normalized');
                        tau_lag_spd = cor_ut(cor_u == max(cor_u)); % lag at max correlation

                        % statistical significance - 2/sqrt(T) (95%)
                        c_min = 2/sqrt(length(spd_i));
                        if max(cor_u) > c_min && tau_lag_spd < 0 && abs(tau_lag_spd) > min_lag
                            % store leaders and followers
                            spd_leader_esp = [spd_leader_esp i];
                            spd_follower_esp = [spd_follower_esp j];
                            spd_corr = [spd_corr, max(cor_u)];
                        end

                    end

                end

            end

        end

        % avg network for a given experiment

        avg_esp_leader = []; % store leader
        avg_esp_follower = []; % store follower
        avg_esp_edges = []; % store no.of edges
        avg_esp_cc = []; % store cc

        for l = 1:n

            for f = 1:n

                if l ~= f

                    % id where leader is l and follower is f
                    no_l_lead = find(spd_leader_esp == l & spd_follower_esp == f);
                    avg_corr_l = spd_corr(no_l_lead);
                    no_l_lead = length(no_l_lead);
                    avg_edges_lf = (no_l_lead)/no_it;
                    avg_corr_lf = mean(avg_corr_l);

                    if avg_edges_lf > 0
                        avg_esp_leader = [avg_esp_leader, l];
                        avg_esp_follower = [avg_esp_follower, f];
                        avg_esp_edges = [avg_esp_edges, avg_edges_lf];
                        avg_esp_cc = [avg_esp_cc, avg_corr_lf];
                    end

                end

            end

        end

        avg_exp_net_temp = [avg_esp_leader; avg_esp_follower; avg_esp_cc; avg_esp_edges];
        avg_exp_net = [avg_exp_net, avg_exp_net_temp];

        % struct_graph = digraph(avg_esp_leader, avg_esp_follower, avg_esp_edges, n);
        % plot(struct_graph, 'Marker', 'o', 'NodeColor', '#2c7fb8', 'MarkerSize', ...
        %     1, 'LineWidth', 2, 'NodeLabel', 1:n, 'EdgeLabel', struct_graph.Edges.Weight, ...
        %    'EdgeColor', '#bcbddc', ...
        %    'ArrowSize', 8, 'NodeFontSize', 12, 'NodeFontName', 'Arial')

    end

    % avg network from all experiment

    avg_esp_leader = [];
    avg_esp_follower = [];
    avg_esp_cc = [];
    avg_esp_edges = [];
    no_lead_cutoff = 0.05;

    for l = 1:n

        for f = 1:n

            if l ~= f

                % id where leader is l and follower is f
                no_l_lead = avg_exp_net(1,:) == l & avg_exp_net(2,:) == f;
                avg_l_cc = avg_exp_net(3,no_l_lead); % mean cc
                no_l_lead = sum(avg_exp_net(4,no_l_lead)); % normalised edges
                no_l_lead = no_l_lead/no_exp; % normalised edges

                if no_l_lead > no_lead_cutoff
                    avg_esp_leader = [avg_esp_leader, l];
                    avg_esp_follower = [avg_esp_follower, f];
                    avg_esp_edges = [avg_esp_edges, no_l_lead];
                    avg_esp_cc = [avg_esp_cc, mean(avg_l_cc)];
                end

            end

        end

    end
    
    % figure(1)

    avg_esp_graph = digraph(avg_esp_follower, avg_esp_leader, round(avg_esp_edges,3));
    % plot(avg_esp_graph, 'LineWidth', avg_esp_graph.Edges.Weight + 1, ...
    %     'EdgeLabel', avg_esp_graph.Edges.Weight, 'Layout', 'layered')

    full_net = [avg_esp_leader; avg_esp_follower; avg_esp_cc; avg_esp_edges];
    net_fname = strcat('full_net_stat_phase_', num2str(cal_corr_event), '_', fname, '.csv');
    writematrix(full_net, net_fname)

    % difference plot

    % figure(2)

    net_esp_leader = []; % leader id
    net_esp_follower = []; % follower id
    net_esp_cc = []; % cc
    net_esp_edges = []; % avg no.of edges from follower to leader.
    no_lead_cutoff = 0.0; % draw edge only if no.of edges are greater than this

    for l = 1:n

        for f = 1:n

            if l < f

                % ids where leader is l and follower is f
                no_l_lead = find(avg_exp_net(1,:) == l & avg_exp_net(2,:) == f);
                avg_cc_l = avg_exp_net(3,no_l_lead);
                no_l_lead = sum(avg_exp_net(4,no_l_lead)); % normalised edges;
                % ids where leader is 'f' and follower is 'l'
                no_f_lead = find(avg_exp_net(1,:) == f & avg_exp_net(2,:) == l);
                avg_cc_f = avg_exp_net(3,no_f_lead);
                no_f_lead = sum(avg_exp_net(4,no_f_lead)); % normalised edges

                avg_cc = mean([avg_cc_l, avg_cc_f]); % average cc
                avg_edges_lf = (no_l_lead - no_f_lead)/no_exp; % normalised avg difference edges

                if avg_edges_lf > 0 && abs(avg_edges_lf) > no_lead_cutoff
                    net_esp_leader = [net_esp_leader, l];
                    net_esp_follower = [net_esp_follower, f];
                    net_esp_edges = [net_esp_edges, avg_edges_lf];
                    net_esp_cc = [net_esp_cc, avg_cc];
                elseif avg_edges_lf < 0 && abs(avg_edges_lf) > no_lead_cutoff
                    net_esp_leader = [net_esp_leader, f];
                    net_esp_follower = [net_esp_follower, l];
                    net_esp_edges = [net_esp_edges, abs(avg_edges_lf)];
                    net_esp_cc = [net_esp_cc, avg_cc];
                end

            end

        end

    end

    net_esp_graph = digraph(net_esp_follower, net_esp_leader, round(net_esp_edges,3));
    % plot(net_esp_graph, 'LineWidth', net_esp_graph.Edges.Weight + 1, ...
    %     'EdgeLabel', net_esp_graph.Edges.Weight, 'Layout', 'layered')

    diff_net = [net_esp_leader; net_esp_follower; net_esp_cc; net_esp_edges];
    diff_net_fname = strcat('diff_net_stat_phase_', num2str(cal_corr_event), '_', fname, '.csv');
    writematrix(diff_net, diff_net_fname)

    is_stat_fname = strcat('is_stat_phase', num2str(cal_corr_event), '_', fname, '.csv');
    writematrix(is_stationary, is_stat_fname)

end

end_statement_stat = sprintf("cal_corr_event = %i, e = %i", cal_corr_event, e);
writelines(end_statement_stat, strcat('end_statement_stat_', fname, '.txt'))

% toc