function [phi_t, pos_t, s_t, en_start, en_centre, en_end, rank_order_atk] = traj_sim(n, ...
    dt, n_iter, zor, sight, l_wall, box_len, mini_box_len, box_width, s_be, beta, ...
    alpha, D_s, D_phi, k_alg, k_atr, K, mu_alg, mu_d, md, mu_wall, mu_esp, ...
    no_inf_neighbours, gamma, omega_ini, t_atk)

% initial condition. putting all agents on the left side of the dabba.
pos(:,1) = zor*rand(n,1);
pos(:,2) = zor*rand(n,1);

% position of centre of box of the right box. after perturbatrion
% trained individuals tend to move in this direction.
pos_cen_box = [(2*box_len - mini_box_len)/2, box_width/2];

% initial randon orientation before the start of perturbation
phi = 2*pi*rand(n,1); % initial heading angle

% initial speed
s = s_be + 0.1*randn(n,1); % initial speed
s(s < 0) = s_be; % ensuring that speed >= 0

pos_t_1 = pos; % pos at time, t-1.
phi_t_1 = phi; % orientation at time, t-1.
s_t_1 = s; % speed at time, t-1.

% stores position, speed and orientation over all time steps
pos_t = zeros(n,2,n_iter);
s_t = zeros(n,n_iter);
phi_t = zeros(n,n_iter);

% storing t = 1
pos_t(:,:,1) = pos;
phi_t(:,1) = phi;
s_t(:,1) = s;

% time of crossing the boarders
en_start = nan(n,1); % time at which agents cross first barrier
en_centre = nan(n,1); % time at which agents cross centre barrier
en_end = nan(n,1); % time at which agents cross the last barrier

rank_order_atk = nan(1,n); % order at which they cross the first barrier
rank_count = 1;

omega = zeros(n,1); % weight of escape force. omega_ini for conditioned fish
omega(1:no_inf_neighbours) = omega_ini;

for t = 2:n_iter

    % if rem(t, 1e4) == 0
    %     disp(t*dt);
    % end

    for i = 1:n

        e_v = [cos(phi_t_1(i)) sin(phi_t_1(i))]; % current velocity heading
        e_phi = [-sin(phi_t_1(i)) cos(phi_t_1(i))]; % perpendicular to velocity heading

        % speed noise
        spd_noise = sqrt(2*D_s)*randn(1);

        % angular noise
        angular_noise = sqrt(2*D_phi)*randn(1);

        % find the nearest k neighbour and 
        % repulsion with neighbours
        dis_vect = pos_t_1 - repmat(pos_t_1(i,:), size(pos_t_1,1), 1); % rij
        mag_rij = sqrt(dis_vect(:,1).^2 + dis_vect(:,2).^2); % |rij|
        rij_hat = dis_vect./(mag_rij + eps);
        
        phi_ij = dot(repmat(e_v,n,1),rij_hat,2);
        phi_ij = acos(phi_ij); % angle between rij and e_v
        phi_ij(i) = pi + 0.1;
        vis_neigh_id = find(phi_ij < sight); % identify neighbours within visual zone.
        [~, near_neigh_id] = sort(mag_rij(vis_neigh_id)); % sort visible neighbours based on their distance. 
        near_neigh_id = near_neigh_id(1:min(K,length(near_neigh_id)));
        % chose k randomly from K neighbours. 
        alg_neigh_id = near_neigh_id(randperm(length(near_neigh_id), min(k_alg, length(near_neigh_id)))); % randomly choosing k out of K neighbours
        alg_neigh_id = vis_neigh_id(alg_neigh_id); 
        % atr_neigh_id = near_neigh_id(randperm(K,k_atr));
        atr_neigh_id  = alg_neigh_id;

        % force from wall
        f_wall = rep_wall_force(pos_t_1, phi_t_1, i, no_inf_neighbours, ...
            box_len, mini_box_len, box_width, l_wall, mu_wall, t, t_atk);

        if t < t_atk % before the light is turned on
            
            % alignment and attraction force
            f_alg = alignment_force(s_t_1, phi_t_1, i, alg_neigh_id, mu_alg);
            f_atr = attraction_force(pos_t_1, i, atr_neigh_id, mu_d, md, zor);
            
            f_total = f_wall + f_alg + f_atr;

            % total along heading direction
            f_v = dot(f_total, e_v, 2);

            % total along turning direction
            f_phi = dot(f_total, e_phi, 2);

            % equations of motion

            % speed
            s(i) = s(i) + beta*(s_be - s(i))*dt + f_v*dt + spd_noise*sqrt(dt);

            % heading
            phi(i) = phi_t_1(i) + (1/(s_t_1(i) + alpha))*(f_phi*dt + angular_noise*sqrt(dt));

        % when perturbation has started and for conditioned individuals
        elseif t >= t_atk && i <= no_inf_neighbours && pos_t_1(i,1) < box_len - mini_box_len
            
            % alignment and attraction force
            f_alg = alignment_force(s_t_1, phi_t_1, i, alg_neigh_id, mu_alg);
            f_atr = attraction_force(pos_t_1, i, atr_neigh_id, mu_d, md, zor);
            f_social = f_alg + f_atr;
            
            % escape force
            f_esp = pos_cen_box - pos_t_1(i,:);
            f_esp = mu_esp*(f_esp/(vecnorm(f_esp)+eps));
            
            % total force weighted by omega (preference for social interaction)
            f_total = f_wall + (1 - omega(i))*f_social + omega(i)*f_esp;

            % total along heading direction
            f_v = dot(f_total, e_v, 2);

            % total along turning direction
            f_phi = dot(f_total, e_phi, 2);
            
            % equations of motion
            
            % speed
            spd_noise = sqrt(2*D_s)*randn(1);
            s(i) = s(i) + beta*(s_be - s(i))*dt + f_v*dt + spd_noise*sqrt(dt);

            % orientation
            phi(i) = phi_t_1(i) + (1/(s_t_1(i) + alpha))*(f_phi*dt + angular_noise*sqrt(dt));
             
            if ismember(i, rank_order_atk) == 0 && pos_t_1(i,1) > mini_box_len
                rank_order_atk(rank_count) = i;
                rank_count = rank_count + 1;
                en_start(i) = t;
            end

            if isnan(en_centre(i)) && pos_t_1(i,1) > (box_len/2)
                en_centre(i) = t;
            end
        
        % if light is turned on and the conditioned fish has reached the
        % other side
        elseif t > t_atk && i <= no_inf_neighbours && pos_t_1(i,1) > box_len - mini_box_len  

            if isnan(en_end(i)) && pos_t_1(i,1) > box_len - mini_box_len
                en_end(i) = t;
            end
            
            % alignment and attraction force
            f_alg = alignment_force(s_t_1, phi_t_1, i, alg_neigh_id, mu_alg);
            f_atr = attraction_force(pos_t_1, i, atr_neigh_id, mu_d, md, zor);
            f_social = f_alg + f_atr;
            
            % escape force
            f_esp = pos_cen_box - pos_t_1(i,:);
            f_esp = mu_esp*(f_esp/(vecnorm(f_esp)+eps));
            
            % total force weighted by omega
            f_total = f_wall + (1 - omega(i))*f_social + omega(i)*f_esp;
            % rate of decrease of omega
            omega(i) = omega(i)*(1 - gamma*dt);

            % total force along speed
            f_v = dot(f_total, e_v, 2);

            % total force along heading direction
            f_phi = dot(f_total, e_phi, 2);

            % equations of motion

            % speed
            s(i) = s(i) + beta*(s_be - s(i))*dt + f_v*dt + spd_noise*sqrt(dt);

            % heading
            phi(i) = phi_t_1(i) + (1/(s_t_1(i) + alpha))*(f_phi*dt + angular_noise*sqrt(dt));
        
        % if light is turned on and for naive individuals
        else

            if ismember(i, rank_order_atk) == 0 && pos_t_1(i,1) > mini_box_len
                rank_order_atk(rank_count) = i;
                rank_count = rank_count + 1;
                en_start(i) = t;
            end

            if isnan(en_centre(i)) && pos_t_1(i,1) > (box_len/2)
                en_centre(i) = t;
            end

            if isnan(en_end(i)) && pos_t_1(i,1) > box_len - mini_box_len
                en_end(i) = t;
            end
            
            % alignment and attraction force
            f_alg = alignment_force(s_t_1, phi_t_1, i, alg_neigh_id, mu_alg);
            f_atr = attraction_force(pos_t_1, i, atr_neigh_id, mu_d, md, zor);
            f_total = f_alg + f_atr + f_wall;

            % total force along speed
            f_v = dot(f_total, e_v, 2);

            % total force along heading direction
            f_phi = dot(f_total, e_phi, 2);

            % equations of motion

            % speed
            s(i) = s(i) + beta*(s_be - s(i))*dt + f_v*dt + spd_noise*sqrt(dt);

            % heading
            phi(i) = phi_t_1(i) + (1/(s_t_1(i) + alpha))*(f_phi*dt + angular_noise*sqrt(dt));

        end

        % if isempty(nn) == 0
        %     rij = mean(dis_vect(nn,:),1);
        %     rij = -rij/(sqrt(rij(1,1)^2 + rij(1,2)^2)+eps);
        %     theta_d_rand = atan2(sin(theta_d_rand) + rij(1,2), cos(theta_d_rand) + rij(1,1));
        % end

        if s(i) < 0
            s(i) = 0;
        end
        
        if phi(i) < 0 
            phi(i) = phi(i) + 2*pi;
        elseif phi(i) > 2*pi
            phi(i) = phi(i) - 2*pi;
        end
        pos(i,:) = pos_t_1(i,:) + dt*s(i)*[cos(phi(i)) sin(phi(i))];

    end

    phi_t_1 = phi;
    pos_t_1 = pos;
    s_t_1 = s;

    pos_t(:,:,t) = pos;
    phi_t(:,t) = phi;
    s_t(:,t) = s;

end