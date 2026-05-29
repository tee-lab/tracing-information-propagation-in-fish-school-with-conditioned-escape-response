close all
clear
clc

%%

load('dir_walk_n5_atk.mat') % load position, speed and heading data
vel_t = zeros(n,2,n_iter,no_it);
vel_t(:,1,:,:) = cos(theta_t);
vel_t(:,2,:,:) = sin(theta_t);

tm_delay = 20;
cmin = 0.5;

%%

fig_net = 1;
fig_ccf = 1;
node_indeg = nan(no_it, n);
ccf_lag = [];

for i = 1:no_it

    en_start = min(en_time(:,i));
    en_start = en_start - 50;
    en_end = max(en_time(:,i));
    en_end = en_end + 800;

    vx = squeeze(vel_t(:,1,en_start:en_end,i));
    vy = squeeze(vel_t(:,2,en_start:en_end,i));

    agent_i = [];
    agent_j = [];
    wghts = [];

    ccf_lag_temp = nan(n, n);

    for ind = 1:n

        vx_i = vx(ind,:);
        vy_i = vy(ind,:);

        for j = 1:n

            if j ~= ind

                vx_j = vx(j,:);
                vy_j = vy(j,:);

                [ccf_vix_vjx_temp, ~] = xcorr(vx_i, vx_j, tm_delay, 'unbiased');
                [ccf_viy_vjy_temp, tlag] = xcorr(vy_i, vy_j, tm_delay, 'unbiased');

                ccf_vel_temp = (ccf_vix_vjx_temp + ccf_viy_vjy_temp);
                ccf_vel_temp_abs = abs(ccf_vel_temp);
                tlag_id = find(ccf_vel_temp_abs == max(ccf_vel_temp_abs));

                if max(ccf_vel_temp_abs) >= cmin && tlag(tlag_id) < 0

                    ccf_lag_temp(ind,j) = tlag(tlag_id);
                    ccf_lag_temp(j,ind) = -tlag(tlag_id);
                    agent_i = [agent_i ind];
                    agent_j = [agent_j j];
                    wghts = [wghts abs(tlag(tlag_id))];

%                     set(0, 'CurrentFigure', figure(fig_ccf))
%                     subplot(ceil(n/2),2,ind)
%                     plot(tlag, ccf_vel_temp, '-', 'LineWidth', 0.5)
%                     hold on

                end

            end

        end

        % axis([-max(tlag) max(tlag) -1 1])
        % hold off

    end

    if isempty(agent_i) == 0
        set(0, 'CurrentFigure', figure(fig_net))
        struct_graph = digraph(agent_j, agent_i, wghts, n);
        plot(struct_graph, 'Layout', 'layered', 'Marker', 's', 'NodeColor', 'r', 'MarkerSize', 10, ...
        'EdgeColor', [0.5 0.5 0.5], 'ArrowSize', 8, 'NodeFontSize', 13, 'LineWidth', 0.75,...
        'EdgeLabel', wghts)
        node_indeg(i,:) = indegree(struct_graph)';
    end

    fig_net = fig_net + 1;
    fig_ccf = fig_ccf + 1;

    ccf_lag = cat(3, ccf_lag, ccf_lag_temp);

end