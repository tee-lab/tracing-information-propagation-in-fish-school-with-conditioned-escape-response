function f_wall = rep_wall_force(pos_t_1, phi_t_1, i, no_inf_neighbours, ...
    box_len, mini_box_len, box_width, l_wall, mu_wall, t, t_atk)

all_e_desired = [1 0; 0 1];
e_i = [cos(phi_t_1(i)) sin(phi_t_1(i))];

% disp(t)
if t < t_atk || (t >= t_atk && i > no_inf_neighbours && max(pos_t_1(1:no_inf_neighbours,1)) < box_width && pos_t_1(i,1) < box_len - mini_box_len) % before perturbation

    % check if there is any repulsion from wall.
    % [x-0, l-x, y, w-l]
    dist_wall = [pos_t_1(i,1), (mini_box_len - pos_t_1(i,1)), ...
        pos_t_1(i,2), (box_width - pos_t_1(i,2))];
    
    [min_dist_wall, min_dist_wall_id] = min(dist_wall);

    if min_dist_wall > l_wall

        f_wall = [0, 0];

    else

        if min_dist_wall_id == 1 && e_i(1) < 0

            e_d = all_e_desired(2,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        elseif min_dist_wall_id == 2 && e_i(1) > 0

            e_d = all_e_desired(2,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        elseif min_dist_wall_id == 3 && e_i(2) < 0

            e_d = all_e_desired(1,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        elseif min_dist_wall_id == 4 && e_i(2) > 0

            e_d = all_e_desired(1,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        else
            f_wall = [0, 0];
        end

    end

elseif t >= t_atk && (i <= no_inf_neighbours || max(pos_t_1(1:no_inf_neighbours,1)) > box_width) && pos_t_1(i,1) < box_len - mini_box_len % start of perturbation
    
    % no side wall repulsion in this case. 
    dist_wall = [pos_t_1(i,1), l_wall + 1, pos_t_1(i,2), (box_width - pos_t_1(i,2))];

    [min_dist_wall, min_dist_wall_id] = min(dist_wall);

    if min_dist_wall > l_wall

        f_wall = [0, 0];

    else

        if min_dist_wall_id == 1 && e_i(1) < 0

            e_d = all_e_desired(2,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

           f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        elseif min_dist_wall_id == 3 && e_i(2) < 0

            e_d = all_e_desired(1,:);
            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        elseif min_dist_wall_id == 4 && e_i(2) > 0

            e_d = all_e_desired(1,:);
            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        else
            f_wall = [0, 0];
        end

    end

elseif t > t_atk && pos_t_1(i,1) > box_len - mini_box_len  % relax

    dist_wall = [(pos_t_1(i,1) - (box_len - mini_box_len)), (box_len - pos_t_1(i,1)), ...
        pos_t_1(i,2), (box_width - pos_t_1(i,2))];

    [min_dist_wall, min_dist_wall_id] = min(dist_wall);

    if min_dist_wall > l_wall

        f_wall = [0, 0];

    else

        if min_dist_wall_id == 1 && e_i(1) < 0

            e_d = all_e_desired(2,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        elseif min_dist_wall_id == 2 && e_i(1) > 0

            e_d = all_e_desired(2,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        elseif min_dist_wall_id == 3 && e_i(2) < 0

            e_d = all_e_desired(1,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        elseif min_dist_wall_id == 4 && e_i(2) > 0

            e_d = all_e_desired(1,:);
            if dot(e_i, e_d) < 0
                e_d = -e_d;
            end

            f_wall = mu_wall*exp(-(min_dist_wall/l_wall))*(e_d - e_i);

        else
            f_wall = [0, 0];
        end

    end

end