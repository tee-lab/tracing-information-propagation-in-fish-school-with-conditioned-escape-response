close all
clear
clc

%
tic

fname = dir('*.mat'); % file name here
fname = fname.name;
load(fname) % load data

%rank_order_atk

min(min(min(pos_t(:,1,:,:))))
max(max(max(pos_t(:,1,:,:))))
max(max(max(pos_t(:,2,:,:))))
min(min(min(pos_t(:,2,:,:))))

% plotting speed trajectories

iter = ceil(no_it*rand());
exp = 1;
disp(iter)
en_end_escape = max(en_end(:,iter));

plt_count = 0;
plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [50, 100, 1400, 700];

pos_t = pos_t(:,:,:,iter,exp);
theta_t = theta_t(:,:,iter,exp);
s_t = s_t(:,:,iter,exp);

for i = 1:n

    subplot(3,2,i)
    plot(1:n_iter, s_t(i,:))
    hold on

end

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [50, 100, 1400, 700];

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

plt_count = plt_count + 1;
fig = figure(plt_count);
fig.Position = [50, 100, 1400, 700];

% mo = VideoWriter('sm_video_3', 'MPEG-4');
% mo.FrameRate = 10;
% mo.Quality = 100;
% open(mo);

t_atk = min(t_atk);
t_st = 1;
t_et = size(pos_t,3);

agent_col = {'red', 'green', 'black', 'magenta', 'purple'};

for t = 800:2:2200 %(t_et)

    pos_x = pos_t(:,1,t);
    pos_y = pos_t(:,2,t);
    vel_x = cos(theta_t(:,t));
    vel_y = sin(theta_t(:,t));

    quiver(pos_x, pos_y, vel_x, vel_y, 0.1, 'LineWidth', 1, 'ShowArrowHead','on',...
        'Color', 'k')

    hold on

    plot(pos_x(1), pos_y(1), '.', 'Color', 'red', 'MarkerSize', 25);
    hold on 
    plot(pos_x(2), pos_y(2), '.', 'Color', 'k', 'MarkerSize', 25);
    hold on
    plot(pos_x(3:n), pos_y(3:n), '.', 'Color', 'k', 'MarkerSize', 25);

    hold on

    xline(0, 'Color', 'k', 'LineWidth', 3)
    hold on
    xline(box_len, 'Color', 'k', 'LineWidth', 3)
    hold on
    xline(mini_box_len, 'Color', 'r', 'LineWidth', 3)
    hold on
    xline(box_len - mini_box_len, 'Color', 'r', 'LineWidth', 3)
    hold on
    yline(0, 'Color', 'k', 'LineWidth', 3)
    hold on
    yline(box_width, 'Color', 'k', 'LineWidth', 3)
    hold on
    text(box_width, box_width/2, 'Hurdle start', ...
     'Rotation', 90, 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', 'FontSize', 20, 'FontName', 'Helvetica')
    hold on
    text(box_len - box_width, box_width/2, 'Hurdle end', ...
     'Rotation', 90, 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom', 'FontSize', 20, 'FontName', 'Helvetica')
    hold on
    rectangle('Position', [box_width 0 10 box_width],...
    'FaceColor', "#BCBCBC", 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3)

    % plot(pert_pos(1,1,t), pert_pos(1,2,t), '.', 'Color', 'r', 'MarkerSize', 30);
    % 
    % hold on
    % 
    % viscircles(pert_pos(1,:,t), atk_dist, 'Color', 'r');

    hold off

    if t>=t_atk && t < t_atk+300
        rectangle('Position', [0 0 20 20], 'FaceColor', "#dcffdb", 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3)
    end
    
    if t < t_atk 
        title("Initial phase", 'FontSize', 30)
    elseif t >= t_atk && t <= en_end_escape
        title("Escape phase", 'FontSize', 30)
    else
        title("Relaxation phase", 'FontSize', 30)
    end

    set(gca, 'XLim', [0 box_len], 'YLim', [0, box_width], ...
    'LineWidth', 2, 'Xcolor', 'k', 'YColor', 'k', ...
    'FontSize', 30, 'FontName', 'Helvetica', 'Box', 'on')
    
    xlabel('X (cm)')
    ylabel('Y (cm)')
    % axis('equal')
    drawnow('limitrate')

    % image = getframe(figure(plt_count));
    % writeVideo(mo, image);

end

% close(mo)
toc