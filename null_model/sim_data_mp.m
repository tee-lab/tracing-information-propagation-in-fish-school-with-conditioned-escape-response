close all 
clear
clc

%%

% all measurements are in cm and seconds

% tic

n = 5; % no.of individuals
dt = 0.04; % integration time
T = 250; % total simulation time
n_iter = round(T/dt);  
zor = 5; % zone of repulsion
l_wall = 3;
s_be = 1.2; % speed before detecting threat

g_wall = 30;
g_esp = 1.25;

box_len = 50;
box_width = 20;
mini_box_len = 20;

beta = 0.3;
D_phi = 0.3; % angular noise
D_s = 0.3; % velocity noise
alpha = 1; % turning friction
mu_d = 2.3; 
m_d = 2.5; 

gamma = 0.25; % rate at which escape direction info is lost.

no_it = 39; % no.of iteration (same as no.of trails)
no_exp = 5; % no.of sets of trails

stb_time = 500;

theta_f = nan(n,n_iter-stb_time,no_it,no_exp); % store theta
pos_f = nan(n,2,n_iter-stb_time,no_it,no_exp); % store pos
s_f = nan(n,n_iter-stb_time,no_it,no_exp); % store speed
en_start_f = nan(n,no_it,no_exp); % encounter time
en_centre_f = nan(n,no_it,no_exp); % time crossing of centre of box
en_end_f = nan(n,no_it,no_exp); % time at crossing
rank_order_atk_f = nan(n,no_it,no_exp); % id of inds in the order of detecting threats
rank_order_deter_f = nan(n,no_it,no_exp);

% attack time for all individual.
t_atk = ones(n,1)*round(n_iter/2);
% t_atk = round(n_iter/2);
% t_atk = t_atk + atk_diff*(0:(n-1));
[~, rank_ord_deterministic] = sort(t_atk);

for e = 1:no_exp

    disp(e)

    parfor i = 1:no_it

        [phi_t, pos_t, s_t, en_start, en_end, en_centre, rank_order_atk] = bb_traj_mp(n, ...
            dt, n_iter, zor, s_be, beta, D_s, D_phi, box_len, mini_box_len, ...
            mu_d, m_d, box_width, l_wall, g_wall, g_esp, alpha, t_atk, gamma)

        theta(:,:,i) = phi_t(:,stb_time+1:end);
        pos(:,:,:,i) = pos_t(:,:,stb_time+1:end);
        s(:,:,i) = s_t(:,stb_time+1:end);
        en_start_t(:,i) = en_start - stb_time;
        en_centre_t(:,i) = en_centre - stb_time;
        en_end_t(:,i) = en_end - stb_time;
        rank_order_atk_t(:,i) = rank_order_atk;
        rank_ord_deterministic_t(:,i) = rank_ord_deterministic;

    end

    theta_f(:,:,:,e) = theta;
    pos_f(:,:,:,:,e) = pos;
    s_f(:,:,:,e) = s;
    en_start_f(:,:,e) = en_start_t;
    en_end_f(:,:,e) = en_end_t;
    en_centre_f(:,:,e) = en_centre_t;
    rank_order_atk_f(:,:,e) = rank_order_atk_t;
    rank_order_deter_f(:,:,e) = rank_ord_deterministic_t;

end

n_iter = n_iter - stb_time;
t_atk = t_atk - stb_time;

sdata = struct('pos_t', pos_f, 'theta_t', theta_f, 'n', n, 's_t', s_f, 'dt', dt, ...
    'n_iter', n_iter,'s_be', s_be, 'zor', zor, 'l_wall', l_wall,...
    'no_it', no_it, 'beta', beta, 'D_s', D_s, 'D_phi', D_phi, 'box_len', box_len,...
    'mini_box_len', mini_box_len, 'box_width', box_width, ...
    'g_esp', g_esp, 'g_wall', g_wall, 'en_start', en_start_f, 'alpha', alpha, ...
    'en_end', en_end_f, 'en_centre', en_centre_f, 'rank_order_atk', rank_order_atk_f, ...
    'rank_order_deter', rank_order_deter_f, 'no_exp', no_exp, 't_atk', t_atk, 'gamma', ...
    gamma);

save('vs_nm.mat', '-struct', 'sdata')

% toc
disp("Simulation Complete")