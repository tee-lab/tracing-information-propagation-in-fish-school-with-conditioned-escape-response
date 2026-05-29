function f_alg = alignment_force(s_t_1, phi_t_1, i, alg_neigh_id, mu_alg)

if isempty(alg_neigh_id)
    f_alg = [0,0];
else
    vel = s_t_1.*[cos(phi_t_1) sin(phi_t_1)];
    vel_j = vel(alg_neigh_id,:);
    vel_ji = mu_alg*(vel_j - vel(i,:));
    f_alg = mean(vel_ji,1);
end