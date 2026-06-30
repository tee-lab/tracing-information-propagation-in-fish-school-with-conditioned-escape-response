# Tracing information propagation in a fish school with a conditioned escape response

This repository contains the code for the manuscript:

The codes are tested to run on MATLAB Version: 25.2.0.3177638 (R2025b).

## Raw data

Tracked trajectory data used in figure 1e are available in the folder `/main_text/`. Tracks of 
each fish are labelled as `Clip0036_1C_4N_E4_4_bee1.csv`, `Clip0036_1C_4N_E4_4_bee2.csv`, ..., 
`Clip0036_1C_4N_E4_4_bee5.csv`. The first, second, and third columns in the `.csv` files are 
Frame number, position along the x (cm) and y(cm) axis, respectively. 

Tracked data from all the 53 trials are in `/main_text/all_bd_data.mat`. `all_bd_data.mat` contains
position (`pos_t_ex_x`), speed (`speed_t_ex_x`), velocity (`vel_t_ex_x`), and
frame at which green light was turned on (`frame_light_on_ex_x`) for all the experiments, where 
`x` is the trial number between 1 to 53. 

Data file `all_bd_data.mat` should be in `/main_text` folder to reproduce all the figures 
in the main text and supplementary text. We have already added the file in the folder.

## Code for figure 1e and 2a,b,c

**Trajectories of individual fish, speed, polarisation and dispersion from the representative 
trial of collective escape response:** Run `/main_text/figure_1e_2ace.m` to generate figure 1e and
figure 2a-c.

## Code for figure 2 d-f, figure 3, and figure 4 a-b

Run `main_text/col_beh_ana.m`: this code calculates average polarisation, dispersion and speed
as a function of normalised time *t*. It further calculates probability density functions of 
polarisation, dispersion and speed across initial, escape and relax phases. It finally calculates 
normalised crossing time since green light is turned and normalised time since previous fish crossed
the hurdle as a function of crossing rank. 

## Code for figure 4c

**Leader-follower interaction network during collective collective escape:** 
run `main_text/fig_4c.m` to generate the leader-follower interaction for experimental data.

## Codes for results from model of collective escape dynamics

1. **Simulate the model:** run `/model/sim_data.m` to simulate the model. 
2. **Calculate group properties:** run `model/grp_properties.m` to generate figure 2gh, 
figure 3, figure 4ab, electronic supplementary material figure S7, S8. 
3. **Leader-follower interaction network:** run `/model/pert_test_exp.m` to construct 
leader-follower network described in Section 3b of the main text. `model/pert_test_exp.m` 
generates a `.csv` file with leader-follower data. Run `/model/interaction_net_ana.m` to generate
figure 4c. 
4. **Simulation video:** run `/model/simulate.m` to visualise the collective escape dynamics obtained
from the model. 

## Codes for results from null model of collective escape dynamics

1. **Simulate the model:** run `/null_model/sim_data_mp.m` to simulate the Null model described
in electronic supplementary material section S3.
2. **Calculate group properties:** run `/null_model/grp_properties.m` to generate electronic
supplementary material figure S4. This code also generates mean group speed for null model 
that can be compared with figure 2g of collective escape model. 
3. **Leader-follower interaction network:** run `/null_model/pert_test_exp.m` to construct 
leader-follower network described in Section 3b of the main text. `null_model/pert_test_exp.m` 
generates a `.csv` file with leader-follower data. Run `/model/interaction_net_ana.m` to generate
the network. For null model, we do not find any leader-follower pair, as there are no interactions
between agents. 
4. **Simulation video:** run `/null_model/simulate.m` to visualise the group dynamics obtained
from the null model. 

## Supplementary figures

1. `/main_text/col_beh_ana.m` generates electronic supplementary material figure S2, S3
2. `/sm/fig_s5a.m` generates electronic supplementary material figure S5a
3. `/sm/fig_s5bc.m` generates electronic supplementary material figure S5b-c. 





