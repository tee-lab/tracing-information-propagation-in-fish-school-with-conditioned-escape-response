function [phi_t, pos_t, s_t, en_start, en_end, en_centre, rank_order_atk] = bb_traj_mp(n, ...
    dt, n_iter, zor, s_be, beta, D_s, D_phi, box_len, mini_box_len, ...
    mu_d, m_d, box_width, l_wall, g_wall, g_esp, alpha, t_atk, gamma)

% initial condition. putting all agents on the left side of the dabba.
pos(:,1) = mini_box_len*rand(n,1);
pos(:,2) = box_width*rand(n,1);

% position of centre of box of the right box. after perturbatrion
% individuals tend to move in this direction.
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

pos_t(:,:,1) = pos;
phi_t(:,1) = phi;
s_t(:,1) = s;

% angle of escape
en_start = nan(n,1);
en_centre = nan(n,1);
en_end = nan(n,1);

rank_order_atk = nan(1,n);
rank_count = 1;

omega = ones(n,1); % weight of escape force. 

for t = 2:n_iter

    % if rem(t, 1e4) == 0
    %     disp(t*dt);
    % end

    for i = 1:n

        e_v = [cos(phi_t_1(i)) sin(phi_t_1(i))]; % current velocity heading
        e_phi = [-sin(phi_t_1(i)) cos(phi_t_1(i))];

        % speed noise
        spd_noise = sqrt(2*D_s)*randn(1);

        % angular noise
        angular_noise = sqrt(2*D_phi)*randn(1);

        % repulsion with neighbours and correlated random walk
        dis_vect = pos_t_1 - repmat(pos_t_1(i,:), size(pos_t_1,1), 1); % rij
        mag_rij = sqrt(dis_vect(:,1).^2 + dis_vect(:,2).^2); % |rij|
        mag_rij(i) = 2*zor;

        nn = find(mag_rij < zor); % find agents within zor

        if ~isempty(nn)
            r_rep = dis_vect(nn,:)./mag_rij(nn);
            mag_rij = mag_rij(nn);
            f_rep = mu_d*tanh(m_d*(mag_rij - zor)).*(r_rep);
            f_rep = mean(f_rep,1);
        else
            f_rep = [0, 0];
        end

        f_wall = rep_wall_force(pos_t_1, phi_t_1, i, box_len, mini_box_len, ...
            box_width, l_wall, g_wall, t, t_atk);

        f_total = f_rep + f_wall;

        if t < t_atk(i)

            % wall force along heading direction
            f_v = dot(f_total, e_v, 2);

            % wall force along turning direction
            f_phi = dot(f_total, e_phi, 2);

            % equations of motion

            % speed
            s(i) = s(i) + beta*(s_be - s(i))*dt + f_v*dt + spd_noise*sqrt(dt);

            % heading
            phi(i) = phi_t_1(i) + (1/(s_t_1(i) + alpha))*(f_phi*dt + angular_noise*sqrt(dt));

        % check if perturbation has started and if agents can see it and close to it.
        elseif t >= t_atk(i) && pos_t_1(i,1) < box_len - mini_box_len

            f_esp = pos_cen_box - pos_t_1(i,:);
            f_esp = g_esp*(f_esp/(vecnorm(f_esp)+eps));

            f_total = f_total + f_esp;

            % wall force along heading direction
            f_v = dot(f_total, e_v, 2);

            % wall force along turning direction
            f_phi = dot(f_total, e_phi, 2);
            
            % equations of motion
            
            % speed
            spd_noise = sqrt(2*D_s)*randn(1);
            s(i) = s(i) + beta*(s_be - s(i))*dt + f_v*dt +spd_noise*sqrt(dt);
            % if pos_t_1(i,1) > box_width
            %     s(i) = s(i) + beta*(s_ae - s(i))*dt + f_v*dt + spd_noise*sqrt(dt);
            % else
            %     s(i) = s(i) + beta*(s_be - s(i))*dt + f_v*dt + spd_noise*sqrt(dt);
            % end

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

        elseif t > t_atk(i) && pos_t_1(i,1) > box_len - mini_box_len  % if perturbation has not yet started

            if isnan(en_end(i))
                en_end(i) = t;
            end

            f_esp = pos_cen_box - pos_t_1(i,:);
            f_esp = g_esp*(f_esp/(vecnorm(f_esp)+eps));

            f_total = f_total + (omega(i))*f_esp;
            omega(i) = omega(i)*(1 - gamma*dt);

            % wall force along speed
            f_v = dot(f_total, e_v, 2);

            % wall force along heading direction
            f_phi = dot(f_total, e_phi, 2);

            % equations of motion

            % speed
            s(i) = s(i) + beta*(s_be - s(i))*dt + f_v*dt + spd_noise*sqrt(dt);

            % heading
            phi(i) = phi_t_1(i) + (1/(s_t_1(i) + alpha))*(f_phi*dt + angular_noise*sqrt(dt));

        end

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