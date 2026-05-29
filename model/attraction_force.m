function f_atr = attraction_force(pos_t_1, i, atr_neigh_id, mu_d, md, zor)

if isempty(atr_neigh_id)
    f_atr = [0,0];
else
    rji = pos_t_1(atr_neigh_id,:) - pos_t_1(i,:);
    mag_rji = vecnorm(rji,2,2);
    unit_rji = rji./mag_rji;
    force_strength = mu_d*tanh(md*(mag_rji - zor)).*unit_rji;
    f_atr = mean(force_strength,1);
end