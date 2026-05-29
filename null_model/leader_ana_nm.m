%% Code to calculate speed correlation in null-model (with and without detrending)

close all
clear
clc

%% load data
load('vs_model_nm.mat')
fig_pos = [300 300 1600 1600];
add_enc_time = 100;
t_skip = 10;
smooth_window = 8;
max_lag = 90;
n_temp = n;
% c_min = 0.13;

%%

fig_spd_ts = figure('Position', fig_pos);
t = tiledlayout(ceil(sqrt(no_it)), ceil(sqrt(no_it)));
t.Padding = "compact"; 
t.TileSpacing = "compact";
xlabel(t, 'Time')
ylabel(t, 'Speed')

fig_acf_spd = figure('Position', fig_pos);
t = tiledlayout(ceil(sqrt(no_it)), ceil(sqrt(no_it)));
t.Padding = "compact"; 
t.TileSpacing = "compact";
xlabel(t, '\tau')
ylabel(t, 'AC')

for i = 1:no_it

    s_iter = s_t(rank_order_atk(:,i),:,i);
    en_st_t = min(en_start(:,i));
    en_st_t = max(1, en_st_t - add_enc_time);
    % en_start = 1;
    en_end_t = max(en_end(:,i));
    % en_end_t = en_end_t + add_enc_time;
    en_end_t = min(en_end_t, size(s_iter,2));
    s_iter = s_iter(:,en_st_t:en_end_t);
    s_gaus = smoothdata(s_iter, 2, 'gaussian', smooth_window);
    s_res = s_iter - s_gaus;
    tot_time = size(s_iter,2);

    figure(fig_spd_ts)
    nexttile 
    for ind = 1:n_temp
        
        s_plt = s_iter(ind, 1:tot_time);
        % s_plt = s_res(ind, 1:tot_time);
        s_plt = smoothdata(s_plt, "gaussian", smooth_window);
        plot(1:t_skip:tot_time, s_plt(1:t_skip:tot_time))
        hold on

    end
    xlim([1 tot_time])

    figure(fig_acf_spd)
    nexttile
    for ind = 1:ceil(n)

        s_acf = s_iter(ind, 1:tot_time);
        % s_acf = s_res(ind, 1:tot_time);
        s_acf = s_acf(2:end) - s_acf(1:end-1);
        s_acf = s_acf - mean(s_acf);
        [acf, lags] = xcorr(s_acf, max_lag, 'normalized');
        plot(lags, acf)
        hold on

    end

end

%% identifying leadership 

% addpath('/Users/vivek/Library/CloudStorage/OneDrive-IndianInstituteofScience/IISc/phd/phd_thesis/1c_4n_project/tracked_data');

fig_ccf_spd = figure('Position', fig_pos);
% t = tiledlayout(ceil(sqrt(no_it)), ceil(sqrt(no_it)));
% t.Padding = "compact"; 
% t.TileSpacing = "compact";
% xlabel(t, '\tau')
% ylabel(t, 'C(\tau)')

% fig_ccf_diff = figure('Position', fig_pos);
% t = tiledlayout(ceil(sqrt(no_it)), ceil(sqrt(no_it)));
% t.Padding = "compact"; 
% t.TileSpacing = "compact";
% xlabel(t, '\tau')
% ylabel(t, 'C(\tau)')

% fig_ccf_res = figure('Position', fig_pos);
% t = tiledlayout(ceil(sqrt(no_it)), ceil(sqrt(no_it)));
% t.Padding = "compact"; 
% t.TileSpacing = "compact";
% xlabel(t, '\tau')
% ylabel(t, 'C(\tau)')

spd_leader_esp = [];
spd_follower_esp = [];

for iter = 1:no_it

    s_iter = s_t(rank_order_atk(:,iter),:,iter);
    en_st_t = min(en_start(:,iter));
    % en_st_t = max(1, en_st_t - add_enc_time);
    en_end_t = max(en_end(:,iter));
    % en_end_t = en_end_t + add_enc_time;
    s_iter = s_iter(:,en_st_t:en_end_t);
    s_gaus = smoothdata(s_iter, 2, 'gaussian', smooth_window);
    s_res = s_iter - s_gaus;
    % s_iter = s_res;
    tot_time = size(s_iter,2);

    sig_cc = 0;
    for i = 1:n_temp

        spd_i = s_iter(i,:); % speed of fish i
        spd_i_res = s_res(i,:);
        % spd_i = smoothdata(spd_i, 'gaussian', smooth_window); % smooth the data
        spd_i = spd_i - mean(spd_i, 'omitmissing'); % s - \bar{s}
        spd_i_res = spd_i_res - mean(spd_i_res, 'omitmissing');
        del_spd_i = spd_i(2:end) - spd_i(1:(end-1)); % s(t+1) - s(t)
        del_spd_i = del_spd_i - mean(del_spd_i);

        for j = 1:n_temp

            if j ~= i

                % similarly as above for fish, j
                spd_j = s_iter(j,:);
                spd_j_res = s_res(j,:);
                % spd_j = smoothdata(spd_j, 'gaussian', smooth_window);
                spd_j = spd_j - mean(spd_j, 'omitmissing');
                spd_j_res = spd_j_res - mean(spd_j_res, 'omitmissing');
                del_spd_j = spd_j(2:end) - spd_j(1:(end-1));
                del_spd_j = del_spd_j - mean(del_spd_j);

                % calculate speed cross-correlation for fish i and j.
                [del_spd_cc, tlag] = xcorr(del_spd_i, del_spd_j, max_lag, 'normalized');
                del_tau_lag = tlag(del_spd_cc == max(del_spd_cc));

                % [cor_ut,cor_u,tau_lag,sig_ci]=getCor_scalar(spd_i,spd_j,max_lag);
                [cor_u, cor_ut] = xcorr(spd_i,spd_j,max_lag, 'normalized');
                tau_lag_spd = cor_ut(cor_u == max(cor_u));
                % disp(max(del_spd_cc), del_tau_lag)

                [cor_res, tlag_res] = xcorr(spd_i_res, spd_j_res, max_lag, 'normalized');
                tau_lag_res = tlag_res(cor_res == max(cor_res));
                
                figure(fig_ccf_spd)
                subplot(ceil(sqrt(no_it)), ceil(sqrt(no_it)), iter)
                    
                c_min = 2/sqrt(length(spd_i));
                if max(cor_res) > c_min && tau_lag_res < 0 
                    plot(tlag_res*dt, cor_res)
                    hold on
                    yline(c_min)
                    hold on
                    yline(-c_min)
                    hold on
                    spd_leader_esp = [spd_leader_esp i];
                    spd_follower_esp = [spd_follower_esp j];
                end
                % if max(cor_u) > c_min && tau_lag_spd < 0
                %     plot(cor_ut*dt, cor_u)
                %     hold on
                %     yline(c_min)
                %     hold on
                %     yline(-c_min)
                %     hold on
                %     spd_leader_esp = [spd_leader_esp i];
                %     spd_follower_esp = [spd_follower_esp j];
                % end
                
            end

        end

    end

    % disp(sig_cc)

end

%% avg network

figure(6)

avg_esp_leader = [];
avg_esp_follower = [];
avg_esp_edges = [];

for l = 1:n

    for f = 1:n

        if l < f
            
            no_l_lead = find(spd_leader_esp == l & spd_follower_esp == f);
            no_l_lead = length(no_l_lead);
            no_f_lead = find(spd_leader_esp == f & spd_follower_esp == l);
            no_f_lead = length(no_f_lead);
            avg_edges_lf = no_l_lead - no_f_lead;
            if avg_edges_lf > 0
                avg_esp_leader = [avg_esp_leader, l];
                avg_esp_follower = [avg_esp_follower, f];
                avg_esp_edges = [avg_esp_edges, avg_edges_lf];
            elseif avg_edges_lf < 0
                avg_esp_leader = [avg_esp_leader, f];
                avg_esp_follower = [avg_esp_follower, l];
                avg_esp_edges = [avg_esp_edges, abs(avg_edges_lf)];
            end

        end

    end

end

avg_esp_graph = digraph(avg_esp_follower, avg_esp_leader, avg_esp_edges);
% plot(avg_esp_graph, 'LineWidth', avg_esp_graph.Edges.Weight)
plot(avg_esp_graph)

%% t2 - t1

% fig_entime = figure('Position', fig_pos);
% t = tiledlayout(ceil(sqrt(no_it)), ceil(sqrt(no_it)));
% t.Padding = "compact"; 
% t.TileSpacing = "compact";
% xlabel(t, 'Ind')
% ylabel(t, 't_2 - t_1')
% 
% figure(fig_entime)
% 
% for i = 1:no_it
% 
%     ent_temp = en_start(:,i);
%     ent_temp = sort(ent_temp);
%     diff_ent_temp = ent_temp(2:end) - ent_temp(1:end-1);
%     nexttile
%     plot(2:n, diff_ent_temp)
% 
% end