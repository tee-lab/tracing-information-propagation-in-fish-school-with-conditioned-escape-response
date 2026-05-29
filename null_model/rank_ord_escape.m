close all
clear
clc

%% Load data

load('vs_model_nm.mat')
no_of_cp = 2; % no.of change points - escape and relax

rank_order_escape = nan(size(rank_order_atk));
esp_time = nan(2, no_it, no_exp);
smooth_window = 10;

%%

for e = 1:no_exp

    for iter = 1:no_it

        s_iter = s_t(:,:,iter,e);
        rank_esp_time = nan(n,2);

        for i = 1:n

            s_i = s_iter(i,:);
            s_i = smoothdata(s_i, "gaussian", smooth_window);
            cpt_temp = findchangepts(s_i, MaxNumChanges = no_of_cp, ...
                Statistic = "mean");
            if ~isempty(cpt_temp)
                rank_esp_time(i,:) = cpt_temp;
            end

            subplot(5,1,i)
            plot(1:length(s_i), s_i, 'k')
            hold on
            xline(rank_esp_time(i,1), 'r')
            hold on
            xline(rank_esp_time(i,2), 'r')
            hold on
            xline(t_atk, 'blue')
            hold on
            xline(max(en_end(:,iter,e)), 'blue')
            axis([1 4000 0 5])
            hold off

        end
        
        [~, rank_esp_id] = sort(rank_esp_time(:,1), 'ascend');
        rank_order_escape(:,iter,e) = rank_esp_id;
        esp_time(1,iter,e) = min(rank_esp_time(:,1));
        esp_time(2,iter,e) = max(rank_esp_time(:,2));

    end

end

esp_time_interval = esp_time(2,:,:) - esp_time(1,:,:);
esp_time_interval = squeeze(esp_time_interval);
esp_time_interval = esp_time_interval(:);
rank_ord_diff = rank_order_atk == rank_order_escape;
rank_ord_diff = sum(rank_ord_diff,1);
rank_ord_diff = rank_ord_diff(:);

plt_no = 1;
figure(plt_no)
histogram(esp_time_interval, 'Normalization', 'pdf')

plt_no = plt_no + 1;
figure(plt_no)
histogram(rank_ord_diff, 5, 'Normalization', 'pdf')

order_data = struct('rank_order_escape', rank_order_escape, 'esp_time', esp_time);
save('rank_order_data.mat', '-struct', 'order_data')