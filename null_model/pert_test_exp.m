%% Code to compute Permute-match tests for null-model

close all
clear
clc

%%

% tic

fname = dir('*.mat');
fname = fname.name;
load(fname) % load data
fname = fname(1:end-4);

max_lag = 75; % max lag to calculate correlation
min_lag = 0; % min lag for cc to be considered meaningful
p_val = 0.05; % p value cutoff
add_enc_time_before = 0;
add_enc_time_after = 0;
t_atk = min(t_atk);

quant_cutoff = 0.95;
min_edges = 0.05;

% construct network for which page? 
% total no.of phase = 3; % 1 = initial, 2 = escape and 3 = relax
no_cal_corr_event = 3; 

no_permutations = 10000; % no.of permutations to calculate p_value
% check if the time-series is stationary
is_stationary = nan(no_exp*no_it*n, 1);
is_stat_count = 1;

%% permute test

for cal_corr_event = 2:2

    avg_exp_net = []; % average interaction network from different sets of simulation
    avg_combi_net = [];
    % leader; follower; cc; p-value

    for e = 1:no_exp

        disp(e)
        % store all non-nan correlations between i and j from same trail and from other trails
        cc_pt_trl_nnan_all = [];

        for iter = 1:no_it

            s_iter = s_t(rank_order_atk(:,iter,e),:,iter,e); % speed from a given iteration

            % selecting the phase to construct network
            if cal_corr_event == 1 % initial phase
                en_st_t = 1; % perturbation start time.
                % en_end_t = min(en_start(:,iter,e)) - 1; % escape time.
                en_end_t = t_atk-1;
                en_end_t = en_end_t + add_enc_time_after;
            elseif cal_corr_event == 2 % escape phase
                % en_st_t = min(en_start(:,iter,e)); % perturbation start time.
                en_st_t = t_atk;
                en_st_t = en_st_t - add_enc_time_before;
                en_end_t = max(en_end(:,iter,e)); % escape time.
                en_end_t = en_end_t + add_enc_time_after;
            elseif cal_corr_event == 3 % relax phase
                en_st_t = max(en_end(:,iter,e)) + 1; % escape time.
                en_st_t = en_st_t + add_enc_time_after;
                en_end_t = size(s_iter,2);
            end

            s_iter = s_iter(:,en_st_t:en_end_t); % speed between fist and last encounter

            cc_pt_trail = []; % permute tests for a given trail

            for i = 1:n

                % permute test for a given individual
                cc_pt = nan(n, no_it+3);
                cc_pt(:,1) = i; % first column gives individual id

                s_i = s_iter(i,:); % speed of a given individual
                f_count = 5; % from where to save permute cc
                is_stationary(is_stat_count) = adftest(s_i);
                is_stat_count = is_stat_count + 1;
                % 5 because, col 1 is ind id i, col 2 is j from within the trail, col
                % 3 is lag, 4 is cc between i and j. col 5 onwards, permute match test

                % permute test for all other trails except trail 'iter'
                for f_pt = 1:no_it

                    if f_pt ~= iter

                        % speed of all inds from different trail
                        s_it_per = s_t(rank_order_atk(:,f_pt,e),:,f_pt,e);

                        if cal_corr_event == 1 % initial phase
                            en_srt_per = 1; % perturbation start time.
                            % en_end_per = min(en_start(:,f_pt,e)) - 1; % escape time.
                            en_end_per = t_atk-1;
                            en_end_per = en_end_per + add_enc_time_after;
                        elseif cal_corr_event == 2 % escape phase
                            % en_srt_per = min(en_start(:,f_pt,e)); % perturbation start time.
                            en_srt_per = t_atk;
                            en_srt_per = en_srt_per - add_enc_time_before;
                            en_end_per = max(en_end(:,f_pt,e)); % escape time.
                            en_end_per = en_end_per + add_enc_time_after;
                        elseif cal_corr_event == 3 % relax phase
                            en_srt_per = max(en_end(:,f_pt,e)) + 1; % escape time.
                            en_srt_per = en_srt_per + add_enc_time_after;
                            en_end_per = size(s_it_per,2);
                        end

                        s_it_per = s_it_per(:,en_srt_per:en_end_per);

                        % take the minimum time between the 2 escape times
                        pert_time = min(length(s_i), size(s_it_per,2));
                        s_1 = s_i(1:pert_time);
                        s_it_per = s_it_per(:,1:pert_time);

                        % calculate and store cross-correlation between ind i and
                        % ind j' from different trail
                        cc_pt_temp = nan(1,n);
                        for ind_p = 1:n
                            s_i_per = s_it_per(ind_p,:);
                            [pert_corr,~] = crosscorr(s_1, s_i_per, "NumLags", min(max_lag, length(s_i_per)-1));
                            cc_pt_temp(ind_p) = max(pert_corr);
                        end

                        % take the max of correlations between ind i from a given trail
                        % and ind j'from other trail
                        % cc_pt(:,f_count) = max(cc_pt_temp);
                        % cc_pt(:,f_count) = cc_pt_temp(1 + randperm(n-1,1));
                        cc_pt(:,f_count) = cc_pt_temp;
                        % if i == 1
                        %     cc_pt(:,f_count) = cc_pt_temp(1 + randperm(n-1,1));
                        % else
                        %     cc_pt(1,f_count) = cc_pt_temp(1);
                        %     cc_pt(2:n,f_count) = cc_pt_temp(1 + randperm(n-1,1));
                        % end
                        f_count = f_count + 1;
                        % cc_pt(:,f_count:(f_count+n-1)) = cc_pt_temp.*ones(n);
                        % f_count = f_count + n;

                    end

                end

                % calculate correlation between ind i and j from the same trail

                for j = 1:n

                    if  i ~= j

                        s_j = s_iter(j,:);
                        [cc, ilag] = crosscorr(s_i, s_j, "NumLags", min(max_lag, length(s_1)-1));
                        if ilag(cc == max(cc)) > 0 && abs(ilag(cc == max(cc))) > min_lag
                            cc_pt(j,2) = j;
                            cc_pt(j,3) = ilag(cc == max(cc));
                            cc_pt(j,4) = max(cc);
                        end

                    end

                end

                % concatanate accross all inds
                cc_pt_trail = cat(1,cc_pt_trail, cc_pt);

            end

            nan_id = ~isnan(cc_pt_trail(:,4)); % non-nan cc
            cc_pt_trl_nnan = cc_pt_trail(nan_id,:);
            cc_pt_trl_nnan_all = cat(1, cc_pt_trl_nnan_all, cc_pt_trl_nnan);

        end

        % calculating pvalue

        cc_nan_id = abs(cc_pt_trl_nnan_all(:,3)) > min_lag;
        cc_pt_trl_nnan_all = cc_pt_trl_nnan_all(cc_nan_id,:);

        surr_quantile = cc_pt_trl_nnan_all(:,5:end);
        surr_quantile = quantile(surr_quantile, quant_cutoff, 2);
        p_cc_pt = cc_pt_trl_nnan_all(:,4) - surr_quantile;

        % p_cc_pt = cc_pt_trl_nnan_all(:,5:end) - cc_pt_trl_nnan_all(:,4);
        % p_cc_pt(p_cc_pt < 0) = 0;
        % p_cc_pt(p_cc_pt > 0) = 1;
        % pval_cc_pt = mean(p_cc_pt, 2);

        % pval_id = find(pval_cc_pt < p_val);
        pval_id = find(p_cc_pt > 0);
        pval_cc_net_all = cc_pt_trl_nnan_all(pval_id,:);

        pval_cc_net = pval_cc_net_all(:,1:4);

        % avg network for a given experiment
        l_id = []; % store leader
        f_id = []; % store follower
        avg_edges = []; % store no.of edges
        avg_cc = []; % store cc

        for l = 1:n

            for f = 1:n

                if l ~= f

                    no_l_lead = find(pval_cc_net(:,1) == l & pval_cc_net(:,2) == f);
                    avg_corr_l = mean(pval_cc_net(no_l_lead,4));
                    no_l_lead = length(no_l_lead);
                    avg_edges_lf = (no_l_lead)/no_it;

                    if avg_edges_lf > min_edges
                        l_id = [l_id, l];
                        f_id = [f_id, f];
                        avg_edges = [avg_edges, avg_edges_lf];
                        avg_cc = [avg_cc, avg_corr_l];
                    end

                end

            end

        end

        avg_net_temp = [l_id; f_id; avg_edges; avg_cc];
        avg_combi_net = [avg_combi_net, avg_net_temp];

        if ~isempty(l_id)
            struct_graph = digraph(f_id, l_id, avg_edges, n);
            plot(struct_graph, 'Marker', 'o', 'NodeColor', '#2c7fb8', 'MarkerSize', ...
                1, 'LineWidth', 2, 'NodeLabel', 1:n, 'EdgeLabel', struct_graph.Edges.Weight, ...
                'EdgeColor', '#bcbddc', ...
                'ArrowSize', 8, 'NodeFontSize', 12, 'NodeFontName', 'Arial')
        end

        % sum(pval_cc_pt < 0.05)/length(pval_cc_pt)
        % disp(avg_net_temp)
        % sum(sum(isnan(cc_pt_trl_nnan_all)))

        % Constructing avg network using Permutation test and monte carlo

        avg_esp_leader = []; % leader ids
        avg_esp_follower = []; % follower ids
        avg_esp_pval = []; % pval
        avg_esp_cc = []; % max correlation
        col_cc_pair = 4; % column with corr btw i and j from same traile

        for i = 1:n

            for j = 1:n

                if j ~= i

                    pair_id = find(cc_pt_trl_nnan_all(:,1) == i & cc_pt_trl_nnan_all(:,2) == j);

                    if isempty(pair_id)
                        continue
                    end

                    cc_pair_ij = cc_pt_trl_nnan_all(pair_id,:);
                    avg_cc_ij = mean(cc_pair_ij(:,col_cc_pair));
                    pert_cc = nan(1,no_permutations);

                    for p = 1:no_permutations

                        pert_ord = randperm(no_it, size(cc_pair_ij,1)) + (col_cc_pair-1);

                        pert_cc_temp = nan(1,size(cc_pair_ij,1));
                        for k = 1:size(cc_pair_ij,1)
                            pert_cc_temp(k) = cc_pair_ij(k,pert_ord(k));
                        end

                        pert_cc(p) = mean(pert_cc_temp);

                    end

                    pert_cc_diff = pert_cc - avg_cc_ij;
                    pert_cc_diff(pert_cc_diff < 0) = 0;
                    pert_cc_diff(pert_cc_diff > 0) = 1;
                    pval = mean(pert_cc_diff);

                    avg_esp_leader = [avg_esp_leader, i];
                    avg_esp_follower = [avg_esp_follower, j];
                    avg_esp_cc = [avg_esp_cc, avg_cc_ij];
                    avg_esp_pval = [avg_esp_pval, pval];

                end

            end

        end

        id_sig = avg_esp_pval < p_val;
        avg_esp_leader = avg_esp_leader(id_sig);
        avg_esp_follower = avg_esp_follower(id_sig);
        avg_esp_cc = avg_esp_cc(id_sig);
        avg_esp_pval = avg_esp_pval(id_sig);

        avg_exp_net_temp = [avg_esp_leader; avg_esp_follower; avg_esp_cc; avg_esp_pval];
        avg_exp_net = [avg_exp_net, avg_exp_net_temp];

        % avg_pm_graph = digraph(avg_esp_follower, avg_esp_leader, n);
        % p = plot(avg_pm_graph, 'LineWidth', 3, 'ArrowSize', 12, 'Layout', 'layered');
        % labeledge(p, 1:numedges(avg_pm_graph), round(avg_esp_cc,3))

    end

    % Constructing the average network across iterations and experiments

    avg_comb_leader = []; % leader id
    avg_comb_follower = []; % follower id
    avg_comb_cc = []; % cc
    avg_comb_edges = []; % avg no.of edges from follower to leader.
    no_lead_cutoff = 0.0; % draw edge only if no.of edges are greater than this

    for l = 1:n

        for f = 1:n

            if l ~= f

                no_l_lead = find(avg_combi_net(1,:) == l & avg_combi_net(2,:) == f);
                avg_cc = mean(avg_combi_net(4,no_l_lead));
                no_l_lead = sum(avg_combi_net(3,no_l_lead));
                no_l_lead = no_l_lead/no_exp;

                if no_l_lead > no_lead_cutoff
                    avg_comb_leader = [avg_comb_leader, l];
                    avg_comb_follower = [avg_comb_follower, f];
                    avg_comb_edges = [avg_comb_edges, no_l_lead];
                    avg_comb_cc = [avg_comb_cc, avg_cc];
                end

            end

        end

    end

    avg_comb_graph = digraph(avg_comb_follower, avg_comb_leader, round(avg_comb_edges,4));
    plot(avg_comb_graph, 'LineWidth', avg_comb_graph.Edges.Weight, ...
        'EdgeLabel', avg_comb_graph.Edges.Weight, 'Layout', 'layered')

    % figure(1)

    avg_esp_leader = []; % leader id
    avg_esp_follower = []; % follower id
    avg_esp_cc = []; % cc
    avg_esp_edges = []; % avg no.of edges from follower to leader.
    no_lead_cutoff = 0.0; % draw edge only if no.of edges are greater than this

    for l = 1:n

        for f = 1:n

            if l ~= f

                no_l_lead = find(avg_exp_net(1,:) == l & avg_exp_net(2,:) == f);
                avg_cc = mean(avg_exp_net(3,no_l_lead));
                no_l_lead = length(no_l_lead);
                no_l_lead = no_l_lead/no_exp;

                if no_l_lead > no_lead_cutoff
                    avg_esp_leader = [avg_esp_leader, l];
                    avg_esp_follower = [avg_esp_follower, f];
                    avg_esp_edges = [avg_esp_edges, no_l_lead];
                    avg_esp_cc = [avg_esp_cc, avg_cc];
                end

            end

        end

    end

    avg_esp_graph = digraph(avg_esp_follower, avg_esp_leader, round(avg_esp_edges,3));
    % plot(avg_esp_graph, 'LineWidth', avg_esp_graph.Edges.Weight, ...
    %     'EdgeLabel', avg_esp_graph.Edges.Weight, 'Layout', 'layered')

    full_net = [avg_esp_leader; avg_esp_follower; avg_esp_cc; avg_esp_edges];
    net_fname = strcat('full_net_pm_phase_', num2str(cal_corr_event), '_', fname, '.csv');
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
                no_l_lead = length(no_l_lead);
                % ids where leader is 'f' and follower is 'l'
                no_f_lead = find(avg_exp_net(1,:) == f & avg_exp_net(2,:) == l);
                avg_cc_f = avg_exp_net(3,no_f_lead);
                no_f_lead = length(no_f_lead);

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
    % plot(net_esp_graph, 'LineWidth', net_esp_graph.Edges.Weight, ...
    %     'EdgeLabel', net_esp_graph.Edges.Weight, 'Layout', 'layered')

    diff_net = [net_esp_leader; net_esp_follower; net_esp_cc; net_esp_edges];
    diff_net_fname = strcat('diff_net_pm_phase_', num2str(cal_corr_event), '_', fname, '.csv');
    writematrix(diff_net, diff_net_fname)

end


% toc