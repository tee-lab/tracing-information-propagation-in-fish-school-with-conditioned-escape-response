close all
clear
clc

%% 

all_bd_data = load("all_bd_data.mat");
tcross_data = load('tcross_data.mat');
load("all_bd_data.mat", "dt", "n", "no_exp", "tank_max_x", "tank_min_x", "tank_middle")
smooth_window = 5;
max_delay = 75; % max delay
min_ts_length = 10;
no_treat = 3;
tdiff = 1;
min_lag = 5;
p_val = 0.05;
quant_cutoff = 0.95;
min_edges = 0.05;

ini_color = "#0096FF";
esp_color = "#EE4B2B";
relax_color = "#097969";

%%

cc_pt_trl_nnan_all = []; % store all non-nan correlations between i and j from same trail and from other trails
sig_cc_net_all = []; % significant correlationss - where correlation within trail is greater than between trails
pval_cc_net_all = [];

tmin = nan(1,no_exp);

for e = 1:no_exp
    
    speed_fish = all_bd_data.(strcat('speed_t_ex_', num2str(e))); % input what speed data to use
    tcross_min = tcross_data.(strcat('tcross_min', num2str(e)));
    tcross_max = tcross_data.(strcat('tcross_max', num2str(e)));
    rank_id_centre = tcross_data.(strcat('rank_id_centre', num2str(e)));
    time_light_on = all_bd_data.(strcat('frame_light_on_ex_', num2str(e)));
    tmin(e) = tcross_min;

    % diff_time = min(min_ts_length,tcross_min - 1);
    % diff_time = ceil(diff_time/4);
    diff_time = min(tdiff, tcross_min-1);
    last_time = tcross_max + diff_time + min_ts_length;
    t_esp = time_light_on; % here change to tcross_min or time_light_on

    speed_trtment = speed_fish(rank_id_centre,1:(t_esp - diff_time));
    speed_escape = speed_fish(rank_id_centre,t_esp:tcross_max);
    speed_pescape = speed_fish(rank_id_centre,(tcross_max + diff_time):last_time);

    % network for speed
    for treat_count = 3:3

        % network during escape
        agent_i_spd = []; % leaders
        agent_j_spd = []; % followers

        if treat_count == 1
            speed_cc = speed_trtment;
        elseif treat_count == 2
            speed_cc = speed_escape;
        elseif treat_count == 3
            speed_cc = speed_pescape;
        end

        if size(speed_cc,2) < min_ts_length
            continue
        end

        cc_pt_trail = []; % permute tests for a given trail

        for i = 1:n
            
            % cc_pt = nan(n,n*(no_exp-1)+4); % permute test for a given individual
            cc_pt = nan(n, no_exp+3);
            cc_pt(:,1) = i; % first column gives individual id
            f_count = 5; % from where to save permute cc
            % 4 because, col 1 is ind id i, col 2 is j from within the trail, col
            % 3 is cc between i and j. col 4 onwards, permute match test

            spd_i = speed_cc(i,:); % speed of fish i
            spd_i = smoothdata(spd_i, 'gaussian', smooth_window); % smooth the data
            spd_i = spd_i - mean(spd_i, 'omitmissing'); % s - \bar{s}

            % permute test for all other trails except trail iter
            for e_pt = 1:no_exp

                if e_pt ~= e

                    speed_fish_pt = all_bd_data.(strcat('speed_t_ex_', num2str(e_pt))); % input what speed data to use
                    tcross_min_pt = tcross_data.(strcat('tcross_min', num2str(e_pt)));
                    tcross_max_pt = tcross_data.(strcat('tcross_max', num2str(e_pt)));
                    rank_id_centre_pt = tcross_data.(strcat('rank_id_centre', num2str(e_pt)));
                    time_light_on_pt = all_bd_data.(strcat('frame_light_on_ex_', num2str(e_pt)));

                    % diff_time_pt = min(min_ts_length,tcross_min_pt - 1);
                    % diff_time_pt = ceil(diff_time_pt/4);
                    diff_time_pt = min(tdiff, tcross_min_pt - 1);
                    last_time_pt = tcross_max_pt + min_ts_length + diff_time_pt;
                    t_esp_pt = time_light_on_pt; % here change to tcross_min or time_light_on

                    if treat_count == 1
                        speed_cc_pt = speed_fish_pt(rank_id_centre_pt,1:(t_esp_pt - diff_time_pt));
                    elseif treat_count == 2
                        speed_cc_pt = speed_fish_pt(rank_id_centre_pt,t_esp_pt:tcross_max_pt);
                    elseif treat_count == 3
                        speed_cc_pt = speed_fish_pt(rank_id_centre_pt,(tcross_max_pt + diff_time_pt):last_time_pt);
                    end
        
                    if size(speed_cc_pt,2) < min_ts_length
                        continue
                    end
                    % take the minimum time between the 2 escape times
                    pert_time = min(length(spd_i), size(speed_cc_pt,2));
                    s_1 = spd_i(1:pert_time);
                    speed_cc_pt = speed_cc_pt(:,1:pert_time);

                    % calculate and store cross-correlation between ind i and
                    % ind j' from different trail
                    cc_pt_temp = nan(1,n);
                    for ind_p = 1:n
                        s_i_per = speed_cc_pt(ind_p,:);
                        s_i_per = smoothdata(s_i_per, 'gaussian', smooth_window); % smooth the data
                        s_i_per = s_i_per - mean(s_i_per);
                        [pert_corr, pt_lag] = crosscorr(s_1, s_i_per, "NumLags", min(max_delay, pert_time-1));
                        [pert_xcorr, pt_xlag] = xcorr(s_1, s_i_per, min(max_delay, pert_time-1), 'normalized');
                        cc_pt_temp(ind_p) = max(pert_corr);
                    end

                    % take the max of correlations between ind i from a given trail
                    % and ind j'from other trail
                    % cc_pt(:,f_count) = max(cc_pt_temp, [], 'omitmissing');
                    % cc_pt(:,f_count) = cc_pt_temp(1 + randperm(n-1,1));
                    cc_pt(:,f_count) = cc_pt_temp;
                    % if i == 1
                    %     cc_pt(:,f_count) = cc_pt_temp(1 + randperm(n-1,1));
                    % else
                    %     cc_pt(1,f_count) = cc_pt_temp(1);
                    %     cc_pt(2:n,f_count) = cc_pt_temp(1 + randperm(n-1,1));
                    % end
                    % cc_pt(:,f_count:(f_count+n-1)) = cc_pt_temp.*ones(n);
                    f_count = f_count + 1;
                    % f_count = f_count + n;

                end

            end


            for j = 1:n

                if i ~= j
                    
                    % similarly as above for fish, j
                    spd_j = speed_cc(j,:);
                    spd_j = smoothdata(spd_j, 'gaussian', smooth_window);
                    spd_j = spd_j - mean(spd_j, 'omitmissing');

                    [cc, ilag] = crosscorr(spd_i, spd_j, "NumLags", min(max_delay, length(spd_i)-1));
                    [cc_x, xlag] = xcorr(spd_i, spd_j, min(max_delay, length(spd_i)-1), 'normalized');
                    if ilag(cc == max(cc)) > 0 && abs(ilag(cc == max(cc))) >= min_lag 
                        cc_pt(j,2) = j;
                        cc_pt(j,3) = ilag(cc == max(cc));
                        cc_pt(j,4) = max(cc);
                    end

                end

            end

             % concatanate accross all inds
             cc_pt_trail = cat(1,cc_pt_trail, cc_pt);

        end

    end
    
    if f_count < no_exp+3
        cc_pt_trail = cc_pt_trail(:,1:f_count-1);
    end
        
    nan_id = ~isnan(cc_pt_trail(:,4)); % non-nan cc
    cc_pt_trl_nnan = cc_pt_trail(nan_id,:);
    cc_pt_trl_nnan_all = cat(1, cc_pt_trl_nnan_all, cc_pt_trl_nnan);
    
    surr_quantile = cc_pt_trl_nnan(:,5:end);
    surr_quantile = quantile(surr_quantile, quant_cutoff, 2);
    % p_cc_pt = cc_pt_trl_nnan(:,5:end) - cc_pt_trl_nnan(:,4);
    p_cc_pt = cc_pt_trl_nnan(:,4) - surr_quantile;
    % p_cc_pt(p_cc_pt < 0) = 0;
    % p_cc_pt(p_cc_pt > 0) = 1;
    
    % only those cc where cc between ind i and j from same trail is greater
    % than cc from permute tests
    sig_cc_id = nan(1,size(cc_pt_trl_nnan,1));
    % pval_cc_pt = mean(p_cc_pt,2);

    % pval_id = find(pval_cc_pt < p_val);
    pval_id = find(p_cc_pt > 0);
    pval_cc_net = cc_pt_trl_nnan(pval_id,:);
    pval_cc_net_all = cat(1,pval_cc_net_all, pval_cc_net);

    for id_c = 1:length(sig_cc_id)
        cc_pt_temp = cc_pt_trl_nnan(id_c,4:end);
        cc_max_id = find(cc_pt_temp == max(cc_pt_temp));
        sig_cc_id(id_c) = cc_max_id;
    end

    sig_cc_id = find(sig_cc_id == 1);
    sig_cc_net = cc_pt_trl_nnan(sig_cc_id,:);
    sig_cc_net_all = cat(1, sig_cc_net_all, sig_cc_net);

end

%% 

pval_cc_net = pval_cc_net_all(:,1:4);
int_net = pval_cc_net;

for i = 1:size(int_net,1)
    if pval_cc_net(i,3) < 0
        int_net(i,1) = pval_cc_net(i,2);
        int_net(i,2) = pval_cc_net(i,1);
        int_net(i,3) = -pval_cc_net(i,3);
    end
end

% avg network for a given experiment
l_id = []; % store leader
f_id = []; % store follower
avg_edges = []; % store no.of edges
avg_cc = []; % store cc

for l = 1:n

    for f = 1:n

        if l ~= f

            no_l_lead = find(int_net(:,1) == l & int_net(:,2) == f);
            avg_corr_l = mean(int_net(no_l_lead,4));
            no_l_lead = length(no_l_lead);
            avg_edges_lf = (no_l_lead)/no_exp;

            if avg_edges_lf > min_edges
                l_id = [l_id, l];
                f_id = [f_id, f];
                avg_edges = [avg_edges, avg_edges_lf];
                avg_cc = [avg_cc, avg_corr_l];
            end

        end

    end

end

if ~isempty(l_id)
    struct_graph = digraph(f_id, l_id, avg_edges, n);
    plot(struct_graph, 'Marker', 'o', 'NodeColor', '#2c7fb8', 'MarkerSize', ...
        6, 'LineWidth', 3, 'NodeLabel', 1:n, 'EdgeLabel', struct_graph.Edges.Weight, ...
        'EdgeColor', '#bcbddc', ...
        'ArrowSize', 12, 'NodeFontSize', 12, 'NodeFontName', 'Arial')
end

size(pval_cc_net,1)/size(cc_pt_trl_nnan_all,1)
sum(sum(isnan(cc_pt_trl_nnan_all)))

%% Constructing avg network using Permutation test and monte carlo

avg_esp_leader = [];
avg_esp_follower = [];
avg_esp_pval = [];
avg_esp_cc = [];
col_cc_pair = 4;
no_permutations = 10000;
avg_esp_cc_lag = [];

for i = 1:n
    
    for j = 1:n
        
        if j ~= i
            
            pair_id = find(cc_pt_trl_nnan_all(:,1) == i & cc_pt_trl_nnan_all(:,2) == j);
            cc_pair_ij = cc_pt_trl_nnan_all(pair_id,:);
            avg_cc_ij = mean(cc_pair_ij(:,col_cc_pair));
            avg_cc_ig_lag = mean(cc_pair_ij(:,col_cc_pair-1));
            
            pert_cc = nan(1,no_permutations);
            for p = 1:no_permutations

                pert_ord = randperm(no_exp, size(cc_pair_ij,1)) + (col_cc_pair-1);

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
            avg_esp_cc_lag = [avg_esp_cc_lag, avg_cc_ig_lag];

        end

    end

end

id_sig = avg_esp_pval <= p_val;
avg_esp_leader = avg_esp_leader(id_sig);
avg_esp_follower = avg_esp_follower(id_sig);
avg_esp_cc = avg_esp_cc(id_sig);
avg_esp_pval = avg_esp_pval(id_sig);
avg_esp_cc_lag = avg_esp_cc_lag(id_sig)*dt;

[avg_esp_leader; avg_esp_follower; avg_esp_cc; avg_esp_cc_lag; avg_esp_pval]

fig_2 = figure(2);
fig_2.Position = [300, 300, 1000, 1000];
avg_pm_graph = digraph(avg_esp_follower, avg_esp_leader);
p = plot(avg_pm_graph, 'LineWidth', 5, 'ArrowSize', 20, 'Layout', 'layered', ...
    'EdgeColor', esp_color, 'NodeColor', esp_color);
% labeledge(p, 1:numedges(avg_pm_graph), round(avg_esp_cc_lag,3))
% labeledge(p, 1:numedges(avg_pm_graph), 1)

exportgraphics(gca, 'pert_test_net.pdf', 'ContentType', 'vector')