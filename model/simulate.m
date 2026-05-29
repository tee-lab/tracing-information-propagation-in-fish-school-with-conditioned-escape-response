close all
clear
clc

%
tic

% fname = dir('*.mat');
% fname = fname.name;
fname = 'test_sim.mat';
load(fname) % load data

%rank_order_atk

%

min(min(min(pos_t(:,1,:,:))))
max(max(max(pos_t(:,1,:,:))))
max(max(max(pos_t(:,2,:,:))))
min(min(min(pos_t(:,2,:,:))))

% plotting speed trajectories

iter = round(no_it*rand());
disp(iter)
iter = 31;
exp = 1;

figure(1)

pos_t = pos_t(:,:,:,iter,exp);
theta_t = theta_t(:,:,iter,exp);
s_t = s_t(:,:,iter,exp);

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

% mo = VideoWriter('n_5_cats', 'MPEG-4');
% mo.FrameRate = 10;
% mo.Quality = 100;
% open(mo);

t_st = 1;
t_et = size(pos_t,3);

agent_col = {'red', 'green', 'black', 'magenta', 'purple'};

for t = 2400:5:(t_et-1300)

    pos_x = pos_t(:,1,t);
    pos_y = pos_t(:,2,t);
    vel_x = cos(theta_t(:,t));
    vel_y = sin(theta_t(:,t));

    quiver(pos_x, pos_y, vel_x, vel_y, 0.1, 'LineWidth', 1, 'ShowArrowHead','on',...
        'Color', 'k')

    hold on

    plot(pos_x(1), pos_y(1), '.', 'Color', 'red', 'MarkerSize', 25);
    hold on 
    plot(pos_x(2), pos_y(2), '.', 'Color', 'magenta', 'MarkerSize', 25);
    hold on
    plot(pos_x(3:n), pos_y(3:n), '.', 'Color', 'k', 'MarkerSize', 25);

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

    % plot(pert_pos(1,1,t), pert_pos(1,2,t), '.', 'Color', 'r', 'MarkerSize', 30);
    % 
    % hold on
    % 
    % viscircles(pert_pos(1,:,t), atk_dist, 'Color', 'r');

    hold off

    tim_disp = strcat('T =  ', num2str(t));
    title(tim_disp)
    
    axis([-2, box_len+2, -2, box_width+2])
    axis('equal')
    drawnow('limitrate')

    % image = getframe(figure(1));
    % writeVideo(mo, image);

end

% close(mo)
toc