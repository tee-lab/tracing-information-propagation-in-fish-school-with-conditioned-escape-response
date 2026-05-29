close all
clear
clc

%
tic

load('vs_nm.mat')

min(min(min(pos_t(:,1,:,:))))
max(max(max(pos_t(:,1,:,:))))
max(max(max(pos_t(:,2,:,:))))
min(min(min(pos_t(:,2,:,:))))

% plotting speed trajectories

iter = round(no_it*rand()); 
exp = 1;
disp(iter)

figure(1)

pos_t = pos_t(:,:,:,iter);
theta_t = theta_t(:,:,iter);
s_t = s_t(:,:,iter);
% rank_order = rank_order_deter(:,iter);
rank_order = rank_order_atk(:,iter);

for i = 1:n
    
    subplot(3,2,i)
    plot(1:n_iter, s_t(i,:))
    hold on

end

figure(2)

for i = 1:n

    subplot(3,2,i)
    plot(squeeze(pos_t(i,1,:)), squeeze(pos_t(i,2,:)))
    hold on
    hold on

    xline(0, 'Color', 'k')
    hold on
    xline(box_len, 'Color', 'k')
    hold on
    xline(mini_box_len, 'Color', 'r')
    hold on
    xline(box_len - mini_box_len, 'Color', 'r')
    hold on
    yline(0, 'Color', 'k')
    hold on
    yline(box_width, 'Color', 'k')

end

%%

figure(3)

mo = VideoWriter('null_model_icts', 'MPEG-4');
mo.FrameRate = 10;
mo.Quality = 100;
open(mo);

t_st = 2400;
t_et = 3300;

agent_col = ["red", "cyan", "#A2142F", "#0072BD",  "#7E2F8E", "k", "magenta"]; % red, blue, brown, dark blue, purple.

for t = t_st:3:t_et

    pos_x = pos_t(:,1,t);
    pos_y = pos_t(:,2,t);
    vel_x = cos(theta_t(:,t));
    vel_y = sin(theta_t(:,t));

    quiver(pos_x, pos_y, vel_x, vel_y, 0.1, 'LineWidth', 1, 'ShowArrowHead','on',...
        'Color', 'k')

    hold on

    for i = 1:n
        
        % col_id = find(rank_order == i);
        plot(pos_x(i), pos_y(i), '.', 'Color', 'k', 'MarkerSize', 25);
        hold on

    end

    xline(0, 'Color', 'k')
    hold on
    xline(box_len, 'Color', 'k')
    hold on
    xline(mini_box_len, 'Color', 'r')
    hold on
    xline(box_len - mini_box_len, 'Color', 'r')
    hold on
    yline(0, 'Color', 'k')
    hold on
    yline(box_width, 'Color', 'k')
    hold off

    tim_disp = strcat('T =  ', num2str(t));
    title(tim_disp)
    
    axis([0, box_len, 0, box_width+5])
    axis('equal')
    drawnow('limitrate')

    image = getframe(figure(1));
    writeVideo(mo, image);

end

close(mo)
toc