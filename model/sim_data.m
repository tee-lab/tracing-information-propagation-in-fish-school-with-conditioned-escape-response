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
l_wall = 3; % distance below which repulsive force of wall is active
s_be = 1.2; % speed before detecting threat

mu_d = 2.3; % repulsion and attraction strength
md = 2.5; % 
mu_alg = 0.1; % alignment strength
mu_wall = 30; % strength of repulstion from wall
mu_esp = 1.25; % strength of escape

box_len = 50; % length of box in cm
box_width = 20; % width of box in cm
mini_box_len = 20; % length of smaller box

beta = 0.3; % speed relaxation coefficient
D_phi = 0.3; % diffusion coefficient (angular noise)
D_s = 0.3; % diffusion coefficient (velocity noise)
alpha = 1; % turning friction
sight = 3*pi/4; % visible range (sight as defined in Couzin et al, 2002)

K = 4; % no.of individuals they can perceive
k_alg = 2; % no of individuals to align to
k_atr = 2; % no of individuals to attract to
no_inf_neighbours = 1; % no.of neighbours that know that they have to escape
gamma = 0.25; % rate at which escape direction info is lost.
omega_ini = 0.7; % initial social interaction (1 - omega)

% attack time for the trained individual.
t_atk = round(n_iter/2);

no_it = 36; % no.of iteration (same as no.of trails)
no_exp = 5; % no.of sets of trails

stb_time = 500;

theta_f = nan(n,n_iter-stb_time,no_it,no_exp); % store theta
pos_f = nan(n,2,n_iter-stb_time,no_it,no_exp); % store pos
s_f = nan(n,n_iter-stb_time,no_it,no_exp); % store speed
en_start_f = nan(n,no_it,no_exp); % encounter time
en_centre_f = nan(n,no_it,no_exp); % time crossing of centre of box
en_end_f = nan(n,no_it,no_exp); % time at crossing
rank_order_atk_f = nan(n,no_it,no_exp); % id of inds in the order of detecting threats

for e = 1:no_exp

    disp(e)

    parfor i = 1:no_it

        [phi_t, pos_t, s_t, en_start, en_centre, en_end, rank_order_atk] = traj_sim(n, ...
            dt, n_iter, zor, sight, l_wall, box_len, mini_box_len, box_width, s_be, beta, ...
            alpha, D_s, D_phi, k_alg, k_atr, K, mu_alg, mu_d, md, mu_wall, mu_esp, ...
            no_inf_neighbours, gamma, omega_ini, t_atk)

        theta(:,:,i) = phi_t(:,stb_time+1:end);
        pos(:,:,:,i) = pos_t(:,:,stb_time+1:end);
        s(:,:,i) = s_t(:,stb_time+1:end);
        en_start_t(:,i) = en_start - stb_time;
        en_centre_t(:,i) = en_centre - stb_time;
        en_end_t(:,i) = en_end - stb_time;
        rank_order_atk_t(:,i) = rank_order_atk;

    end

    theta_f(:,:,:,e) = theta;
    pos_f(:,:,:,:,e) = pos;
    s_f(:,:,:,e) = s;
    en_start_f(:,:,e) = en_start_t;
    en_centre_f(:,:,e) = en_centre_t;
    en_end_f(:,:,e) = en_end_t;
    rank_order_atk_f(:,:,e) = rank_order_atk_t;

end

n_iter = n_iter - stb_time;
t_atk = t_atk - stb_time;

sdata = struct('pos_t', pos_f, 'theta_t', theta_f, 's_t', s_f, 'en_start', en_start_f,...
    'rank_order_atk', rank_order_atk_f, 'en_centre', en_centre_f, 'en_end', en_end_f, ...
    'n', n,  'dt', dt, 'n_iter', n_iter, 'zor', zor, 'l_wall', l_wall, 'box_len', box_len, ...
    'mini_box_len', mini_box_len, 'box_width', box_width, 's_be', s_be, ...
    'beta', beta, 'alpha', alpha, 'D_s', D_s, 'D_phi', D_phi, 'k_alg', k_alg, ...
    'k_atr', k_atr, 'K', K, 'mu_alg', mu_alg, 'mu_d', mu_d, 'md', md, 'mu_wall', mu_wall, ...
    'mu_esp', mu_esp, 'no_inf_neighbours', no_inf_neighbours, 'gamma', gamma,...
    'omega_ini', omega_ini, 'no_it', no_it, 'no_exp', no_exp, 't_atk', t_atk, ...
    'sight', sight);

% fname = strcat('mual_', num2str(mu_alg), '_muat_', num2str(mu_d), ...
%     '_muesp_', num2str(mu_esp), '_K_', num2str(K), '_k_', num2str(k_alg), ...
%     '_omega_ini_', num2str(omega_ini), '.mat');
% save(fname, '-struct', 'sdata')

save('sim.mat', '-struct', 'sdata')

disp("Simulation Complete")
% toc